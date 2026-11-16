#ifndef CODE_BAMoutput
#define CODE_BAMoutput

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
#include "Parameters.h"
#include <zlib.h>

class BAMoutput {//
public:
    //sorted output
    BAMoutput (int iChunk, string tmpDir, Parameters &Pin);
    void coordOneAlign (char *bamIn, uint bamSize, uint iRead);
    void coordBins ();
    void coordFlush ();
    //unsorted output
    BAMoutput (BGZF *bgzfBAMin, Parameters &Pin);
    void unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2);
    void unsortedFlush ();
    void coordUnmappedPrepareBySJout();

    uint32 nBins; //number of bins to split genome into
    uint* binTotalN; //total number of aligns in each bin
    uint* binTotalBytes;//total size of aligns in each bin
private:
    uint64 bamArraySize; //this size will be allocated
    char* bamArray; //large array to store the bam alignments, pre-sorted
    uint64 binSize, binSize1;//storage size of each bin
    uint64 binGlen;//bin genomic length
    char **binStart; //pointers to starts of the bins
    uint64 *binBytes, binBytes1;//number of bytes currently written to each bin
    ofstream **binStream;//output streams for each bin
    BGZF *bgzfBAM;
    Parameters &P;
    string bamDir;

    //unsorted output: each mapping thread compresses its own buffer into BGZF blocks outside
    //mutexOutSAM, and only appends the finished blocks under the lock (see unsortedWrite)
    bool bgzfParallel; //false: old path, bgzf_write under the lock
    int bgzfLevel; //zlib level, taken from bgzfBAM so it matches bgzf_write exactly
    z_stream bgzfZs; //reused for every block (deflateReset gives the same output as a fresh deflateInit2)
    char *bgzfOut; //compressed blocks of the current buffer
    uint64 bgzfOutSize; //allocated size of bgzfOut
    void unsortedWrite(); //write bamArray[0..binBytes1) to bgzfBAM and rewind
};

#endif
