#include "BAMoutput.h"
#include <sys/stat.h>
#include "GlobalVariables.h"
#include <pthread.h>
#include "serviceFuns.cpp"
#include "ThreadControl.h"
#include "streamFuns.h"
#include "ErrorWarning.h"

BAMoutput::BAMoutput (int iChunk, string tmpDir, Parameters &Pin) : P(Pin){//allocate bam array

    nBins=P.outBAMcoordNbins;
    binSize=P.chunkOutBAMsizeBytes/nBins;
    bamArraySize=binSize*nBins;
    bamArray = new char [bamArraySize];

    bamDir=tmpDir+to_string((uint) iChunk);//local directory for this thread (iChunk)

    mkdir(bamDir.c_str(),P.runDirPerm);
    binStart=new char* [nBins];
    binBytes=new uint64 [nBins];
    binStream=new ofstream* [nBins];
    binTotalN=new uint [nBins];
    binTotalBytes=new uint [nBins];
    for (uint ii=0;ii<nBins;ii++) {
        binStart[ii]=bamArray+bamArraySize/nBins*ii;
        binBytes[ii]=0;
        binStream[ii]=&ofstrOpen((bamDir +"/"+to_string(ii)).c_str(), ERROR_OUT, P);    //open temporary files
        binTotalN[ii]=0;
        binTotalBytes[ii]=0;
    };

    binSize1=binStart[nBins-1]-binStart[0];
    nBins=1;//start with one bin to estimate genomic bin sizes
};

BAMoutput::BAMoutput (BGZF *bgzfBAMin, Parameters &Pin) : P(Pin){//allocate BAM array with one bin, streamed directly into bgzf file

    bamArraySize=P.chunkOutBAMsizeBytes;
    bamArray = new char [bamArraySize];
    binBytes1=0;
    bgzfBAM=bgzfBAMin;

    //compress in the mapping threads unless the file is uncompressed or STAR_SERIAL_BGZF=1
    const char *serialEnv=getenv("STAR_SERIAL_BGZF");
    bgzfParallel = bgzfBAM->is_compressed && !(serialEnv!=NULL && serialEnv[0]=='1');
    bgzfOut=NULL;
    bgzfOutSize=0;
    bgzfLevel=bgzfBAM->compress_level;
    if (bgzfParallel) {
        bgzfZs.zalloc=Z_NULL; bgzfZs.zfree=Z_NULL; bgzfZs.opaque=Z_NULL;
        //same parameters as htslib bgzf_compress(): raw deflate (-15), memLevel 8, default strategy
        if (deflateInit2(&bgzfZs, bgzfLevel, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
            bgzfParallel=false;
    };
    //not used
    binSize=0;
    binStream=NULL;
    binStart=NULL;
    binBytes=NULL;
    binTotalBytes=NULL;
    binTotalN=NULL;
    nBins=0;
};

void BAMoutput::unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2) {//record one alignment to the buffer, write buffer if needed

    if (bamSize==0) return; //no output, could happen if one of the mates is not mapped

    if (binBytes1+bamSize2 > bamArraySize) {//write out this buffer
        unsortedWrite();
    };

    memcpy(bamArray+binBytes1, bamIn, bamSize);
    binBytes1 += bamSize;

};

void BAMoutput::unsortedFlush () {//flush all alignments
    unsortedWrite();
};

static inline void bgzfPack16(unsigned char *b, uint32 v) {b[0]=v; b[1]=v>>8;};
static inline void bgzfPack32(unsigned char *b, uint32 v) {b[0]=v; b[1]=v>>8; b[2]=v>>16; b[3]=v>>24;};

void BAMoutput::unsortedWrite () {//write bamArray[0..binBytes1) to the BAM file, rewind the buffer
    //The buffer holds whole BAM records, so the blocks made from it end on a record boundary and can be
    //appended as a unit: the file is a valid BGZF stream of the same records. The BGZF write buffer of
    //bgzfBAM is empty here (the header is flushed by outBAMwriteHeader, and all later writes go through
    //this function), so raw appends cannot interleave with buffered data.
    //Block boundaries differ from bgzf_write (which packs blocks across buffers), so the .bam bytes differ;
    //the decoded records do not. With several threads, record order already depended on thread timing.

    if (!bgzfParallel) {//old path
        if (g_threadChunks.threadBool) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
        bgzf_write(bgzfBAM,bamArray,binBytes1);
        if (g_threadChunks.threadBool) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
        binBytes1=0;
        return;
    };

    if (binBytes1==0) return;

    const uint64 hLen=18, fLen=8; //BGZF block header and footer, as in htslib bgzf.c
    static const unsigned char bgzfMagic[18]={037,0213,010,4,0,0,0,0,0,0377,6,0,0102,0103,2,0,0,0};

    uint64 nBlocks=(binBytes1+BGZF_BLOCK_SIZE-1)/BGZF_BLOCK_SIZE;
    if (bgzfOutSize < nBlocks*BGZF_MAX_BLOCK_SIZE) {//room for every block at its maximum size
        delete [] bgzfOut;
        bgzfOutSize=nBlocks*BGZF_MAX_BLOCK_SIZE;
        bgzfOut=new char [bgzfOutSize];
    };

    uint64 outBytes=0;
    for (uint64 inStart=0; inStart<binBytes1; inStart+=BGZF_BLOCK_SIZE) {//compress one BGZF block, as bgzf_compress() does
        uint32 sLen=(uint32) min((uint64) BGZF_BLOCK_SIZE, binBytes1-inStart);
        unsigned char *dst=(unsigned char*) bgzfOut+outBytes;

        deflateReset(&bgzfZs);
        bgzfZs.next_in=(Bytef*) bamArray+inStart;
        bgzfZs.avail_in=sLen;
        bgzfZs.next_out=dst+hLen;
        bgzfZs.avail_out=BGZF_MAX_BLOCK_SIZE-hLen-fLen;
        if (deflate(&bgzfZs, Z_FINISH) != Z_STREAM_END) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: BGZF block compression failed while writing the BAM file\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
        };
        uint32 dLen=(uint32) (bgzfZs.total_out+hLen+fLen);

        memcpy(dst, bgzfMagic, hLen);
        bgzfPack16(dst+16, dLen-1);
        bgzfPack32(dst+dLen-8, (uint32) crc32(crc32(0L, Z_NULL, 0), (Bytef*) bamArray+inStart, sLen));
        bgzfPack32(dst+dLen-4, sLen);
        outBytes += dLen;
    };

    if (g_threadChunks.threadBool) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
    ssize_t nW=bgzf_raw_write(bgzfBAM, bgzfOut, outBytes);
    if (g_threadChunks.threadBool) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
    if (nW != (ssize_t) outBytes) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not write to the BAM file (disk full?)\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_FILE_WRITE, P);
    };

    binBytes1=0;
};

void BAMoutput::coordOneAlign (char *bamIn, uint bamSize, uint iRead) {

    uint32 *bamIn32;
    uint alignG;
    uint32 iBin=0;

    if (bamSize==0) {
        return; //no output, could happen if one of the mates is not mapped
    } else {
        //determine which bin this alignment belongs to
        bamIn32=(uint32*) bamIn;
        alignG=( ((uint) bamIn32[1]) << 32 ) | ( (uint)bamIn32[2] );
        if (bamIn32[1] == ((uint32) -1) ) {//unmapped
            iBin=P.outBAMcoordNbins-1;
        } else if (nBins>1) {//bin starts have already been determined
            iBin=binarySearch1a <uint64> (alignG, P.outBAMsortingBinStart, (int32) (nBins-1));
        };
    };

    //write buffer is filled
    if (binBytes[iBin]+bamSize+sizeof(uint) > ( (iBin>0 || nBins>1) ? binSize : binSize1) ) {//write out this buffer
        if ( nBins>1 || iBin==(P.outBAMcoordNbins-1) ) {//normal writing, bins have already been determined
            binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
            binBytes[iBin]=0;//rewind the buffer
        } else {//the first chunk of reads was written in one bin, need to determine bin sizes, and re-distribute reads into bins
            coordBins();
            coordOneAlign (bamIn, bamSize, iRead);//record the current align into the new bins
            return;
        };
    };

    //record this alignment in its bin
    memcpy(binStart[iBin]+binBytes[iBin], bamIn, bamSize);
    binBytes[iBin] += bamSize;
    memcpy(binStart[iBin]+binBytes[iBin], &iRead, sizeof(uint));
    binBytes[iBin] += sizeof(uint);
    binTotalBytes[iBin] += bamSize+sizeof(uint);
    binTotalN[iBin] += 1;
    return;
};

void BAMoutput::coordBins() {//define genomic starts for bins
    nBins=P.outBAMcoordNbins;//this is the true number of bins

    //mutex here
    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexBAMsortBins);
    if (P.outBAMsortingBinStart[0]!=0) {//it's set to 0 only after the bin sizes are determined
        //extract coordinates and sort
        uint *startPos = new uint [binTotalN[0]+1];//array of aligns start positions
        for (uint ib=0,ia=0;ia<binTotalN[0];ia++) {
            uint32 *bamIn32=(uint32*) (binStart[0]+ib);
            startPos[ia]  =( ((uint) bamIn32[1]) << 32) | ( (uint)bamIn32[2] );
            ib+=bamIn32[0]+sizeof(uint32)+sizeof(uint);//note that size of the BAM record does not include the size record itself
        };
        qsort((void*) startPos, binTotalN[0], sizeof(uint), funCompareUint1);

        //determine genomic starts of the bins
        P.inOut->logMain << "BAM sorting: "<<binTotalN[0]<< " mapped reads\n";
        P.inOut->logMain << "BAM sorting bins genomic start loci:\n";

        P.outBAMsortingBinStart[0]=0;
        for (uint32 ib=1; ib<(nBins-1); ib++) {
            P.outBAMsortingBinStart[ib]=startPos[binTotalN[0]/(nBins-1)*ib];
            P.inOut->logMain << ib <<"\t"<< (P.outBAMsortingBinStart[ib]>>32) << "\t" << ((P.outBAMsortingBinStart[ib]<<32)>>32) <<endl;
            //how to deal with equal boundaries???
        };
        delete [] startPos;
    };
    //mutex here
    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexBAMsortBins);

    //re-allocate binStart
    uint binTotalNold=binTotalN[0];
    char *binStartOld=new char [binSize1];
    memcpy(binStartOld,binStart[0],binBytes[0]);

    binBytes[0]=0;
    binTotalN[0]=0;
    binTotalBytes[0]=0;

    //re-bin all aligns
    for (uint ib=0,ia=0;ia<binTotalNold;ia++) {
        uint32 *bamIn32=(uint32*) (binStartOld+ib);
        uint ib1=ib+bamIn32[0]+sizeof(uint32);//note that size of the BAM record does not include the size record itself
        coordOneAlign (binStartOld+ib, (uint) (bamIn32[0]+sizeof(uint32)), *((uint*) (binStartOld+ib1)) );
        ib=ib1+sizeof(uint);//iRead at the end of the BAM record
    };
    delete [] binStartOld;
    return;
};

void BAMoutput::coordFlush () {//flush all alignments
    if (nBins==1) {
        coordBins();
    };
    for (uint32 iBin=0; iBin<nBins; iBin++) {
        binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
        binStream[iBin]->flush();
        binBytes[iBin]=0;//rewind the buffer
    };
};

void BAMoutput::coordUnmappedPrepareBySJout () {//flush all alignments
    uint iBin=P.outBAMcoordNbins-1;
    binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
    binStream[iBin]->flush();
    binBytes[iBin]=0;//rewind the buffer
    binStream[iBin]->close();
    binStream[iBin]->open((bamDir +"/"+to_string(iBin)+".BySJout").c_str());
};
