//#include "blocksOverlap.h"
#include "ChimericDetection.h"
#include "ChimericSegment.h"
#include "ReadAlign.h"

int chimericAlignScore (ChimericSegment & seg1, ChimericSegment & seg2)
{
    int chimScore=0;
    uint chimOverlap = seg2.roS>seg1.roS ?  (seg2.roS>seg1.roE ? 0 : seg1.roE-seg2.roS+1) : (seg2.roE<seg1.roS ? 0 : seg2.roE-seg1.roS+1);
    bool diffMates=(seg1.roE < seg1.align.readLength[0] && seg2.roS >= seg1.align.readLength[0]) || (seg2.roE < seg1.align.readLength[0] && seg1.roS >= seg1.align.readLength[0]);

    //segment lengths && (different mates || small gap between segments)
    if (seg1.roE > seg1.P.pCh.segmentMin + seg1.roS + chimOverlap && seg2.roE > seg1.P.pCh.segmentMin + seg2.roS + chimOverlap  \
        && ( diffMates || ( (seg1.roE + seg1.P.pCh.segmentReadGapMax + 1) >= seg2.roS && (seg2.roE + seg1.P.pCh.segmentReadGapMax + 1) >= seg1.roS ) ) )
    {
       chimScore = seg1.align.maxScore + seg2.align.maxScore - (int)chimOverlap; //subtract overlap to avoid double counting
    };

    return chimScore;
};

/////////////////////////////////////////////////////////////
bool ChimericDetection::chimericDetectionMult(uint nW, uint *readLength, int maxNonChimAlignScore, ReadAlign *PEunmergedRA) {

//     for (uint ii=0;ii<chimAligns.size();ii++) {//deallocate aligns
//         if (chimAligns.at(ii).stitchingDone) {//al1,al2 were allocated
//             delete chimAligns.at(ii).al1;
//             delete chimAligns.at(ii).al2;
//         };
//     };

    for (auto cAit=chimAligns.begin(); cAit<chimAligns.end(); cAit++) {//deallocate aligns
        if (cAit->stitchingDone) {//al1,al2 were allocated
            delete cAit->al1;
            delete cAit->al2;
        };
    };

    chimAligns.clear();
    int chimScoreBest=0;
    std::size_t bestChimAlign=0; // points to element of chimAligns with highest chimScoreBest

    int maxPossibleAlignScore = (int)(readLength[0]+readLength[1]);
    int minScoreToConsider = P.pCh.scoreMin;
    if (maxNonChimAlignScore >= minScoreToConsider)
        minScoreToConsider = maxNonChimAlignScore + 1;
    if ((maxPossibleAlignScore - P.pCh.scoreDropMax) > minScoreToConsider)
        minScoreToConsider = maxPossibleAlignScore - P.pCh.scoreDropMax;

    //build each segment once and keep the good ones, in window/align order.
    //A segment depends only on its alignment, which nothing below modifies (stitching works on copies),
    //so pairs i1<i2 of good segments are the same pairs, in the same order, as the nested window/align loops gave.
    segGood.clear();
    for (uint iW=0; iW<nW; iW++) {//cycle windows
        for (uint iA=0; iA<nWinTr[iW]; iA++) {//cycle aligns in the window
            ChimericSegment seg(P,*trAll[iW][iA]);
            if (seg.segmentCheck())
                segGood.push_back(seg);
        };
    };

    for (std::size_t i1=0; i1<segGood.size(); i1++) {//cycle good segments
        ChimericSegment &seg1=segGood[i1];

        for (std::size_t i2=i1+1; i2<segGood.size(); i2++) {//cycle over the following good segments
            ChimericSegment &seg2=segGood[i2];

            if (seg1.str!=0 && seg2.str!=0 && seg2.str!=seg1.str)
                continue; //chimeric segments have to have consistent strands. TODO: make this optional

            int chimScore=chimericAlignScore(seg1,seg2);

            if (chimScore >= minScoreToConsider) {//candidate chimera
                ChimericAlign chAl(seg1, seg2, chimScore, outGen, RA);

                if (!chAl.chimericCheck())
                    continue; //check chimeric alignment

                //re-calculated chimScoreBest includes non-canonical penalty, so the re-calculated score is lower, in some cases it goes to 0 if some checks are not passed
                chAl.chimericStitching(outGen.G, Read1);
                // rescore after stitching.
                if (chAl.chimScore >= minScoreToConsider) { // survived stitching.
                    chimAligns.push_back(chAl);//add this chimeric alignment

                    if (chimAligns.back().chimScore > chimScoreBest) {
                        chimScoreBest=chimAligns.back().chimScore;
                        bestChimAlign = chimAligns.size()-1;
                        if ((chimScoreBest - (int)P.pCh.multimapScoreRange) > minScoreToConsider)
                            // best score increased, so subsequent alignment candidates must score higher
                            minScoreToConsider = chimScoreBest - (int)P.pCh.multimapScoreRange;
                    };
                } // endif stitched chimera survived.
                else {
                    // al1, al2 allocated during stitching
                    delete chAl.al1;
                    delete chAl.al2;
                };

            }; // endif meets chim score criteria
        };//cycle over seg2
    };//cycle over seg1

    if (chimScoreBest==0)
        return false;

    uint chimN=0;
    for (auto cAit=chimAligns.begin(); cAit<chimAligns.end(); cAit++) {
        //scan all chimeras, find the number within score range
        if (cAit->chimScore >= minScoreToConsider)
            ++chimN;
    };

    if (chimN > P.pCh.multimapNmax) //too many loci
        return false;

    uint iTr = 0; //keep track of "HI" SAM attribute
    for (std::size_t i = 0; i < chimAligns.size(); i++) {//output chimeras within score range
        if (chimAligns[i].chimScore >= minScoreToConsider) {

            if (P.pCh.out.junctions)
                chimAligns[i].chimericJunctionOutput(*ostreamChimJunction, chimN, maxNonChimAlignScore, PEunmergedRA != NULL, chimScoreBest, maxPossibleAlignScore);

            if (P.pCh.out.bam) {
                // convert merged SE chimera to PE chimera if this is a merged chimera
                if (PEunmergedRA != NULL) {
                    chimAligns[i].RA = PEunmergedRA;
                    chimAligns[i].RA->peOverlapChimericSEtoPE(chimAligns[i].al1, chimAligns[i].al2, chimAligns[i].al1, chimAligns[i].al2);
                };
                chimAligns[i].chimericBAMoutput(chimAligns[i].al1, chimAligns[i].al2, chimAligns[i].RA, iTr, chimN, i == bestChimAlign, P);
            };
            iTr++;

        };
    };

    return chimN > 0;
};//END
