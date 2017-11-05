//
// Created by 薄井光生 on 2017/10/15.
//

#ifndef SHUFFLE_GENOME_SHUFFLE_GENOME_H
#define SHUFFLE_GENOME_SHUFFLE_GENOME_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <random>

#include <seqan/arg_parse.h>
#include "myseqan.h"
#include "myutil.h"
#include "geneticcode.h"
#include "gtest/gtest.h"

using std::cout;
using std::cerr;
using std::endl;


template<typename TSeq>
void shuffle_base(TSeq &seq, std::mt19937 &mt) {
    unsigned bpLength = seqan::length(seq);
    std::vector<unsigned> v(bpLength);
    std::iota(v.begin(), v.end(), 0);
    std::shuffle(v.begin(), v.end(), mt);

    seqan::DnaString tmp;//TBI: adjust to template
    seqan::resize(tmp, bpLength);
    for (unsigned i = 0; i < bpLength; i++) {
        tmp[i] = seq[v[i]];
    }
    seqan::replace(seq, 0, bpLength, tmp);
    return;
}

template<typename TSeq>
void shuffle_codon(TSeq &seq, std::mt19937 &mt) {
    unsigned bpLength = seqan::length(seq);
    assert(bpLength % 3 == 0);
    unsigned aaLength = bpLength / 3;
    std::vector<unsigned> v(aaLength);
    std::iota(v.begin(), v.end(), 0);
    std::shuffle(v.begin(), v.end(), mt);

    seqan::DnaString tmp;//TBI: adjust to template
    seqan::resize(tmp, bpLength);
    for (unsigned i = 0; i < aaLength; i++) {
        for (unsigned j = 0; j < 3; j++) {
            tmp[3 * i + j] = seq[3 * v[i] + j];
        }
    }
    seqan::replace(seq, 0, bpLength, tmp);
    return;
}

void shuffle_synonymous(seqan::DnaString &seq, GeneticCode const &gc, std::mt19937 &mt) {
    unsigned bpLength = seqan::length(seq);
    assert(bpLength % 3 == 0);
    unsigned aaLength = bpLength / 3;

    for (unsigned i = 0; i < aaLength; i++) {
        seqan::Infix<seqan::DnaString>::Type codon = seqan::infix(seq, 3 * i, 3 * (i + 1));
        gc.synonymous_sub(codon, mt);
    }
    return;
}


//------------------------------------------------------------
// wrapper class for seqan GffRecord
//------------------------------------------------------------
class MyCDS {
    FRIEND_TEST(genome_shuffle, test_classify);

    seqan::GffRecord gff;
    int label;
    int phase = -1; //defined only if this CDS is typical

    int classify(seqan::DnaString const &seq,
                 GeneticCode const &geneticCode) {
        if (gff.endPos > seqan::length(seq)) {
            return 1;
        }

        int length = gff.endPos - gff.beginPos;
        if (length < 6 | (length % 3) != 0) {
            return 2;
        }

        seqan::DnaString subSeq = seqan::infix(seq, gff.beginPos, gff.endPos); //create new object in order to use Infix.

        if (gff.strand == '-') {
            seqan::reverseComplement(subSeq);
        }

        int aaLength = length / 3;
        for (int i = 0; i < aaLength - 1; i++) {
            seqan::Infix<seqan::DnaString>::Type codon = seqan::infix(subSeq, 3 * i, 3 * (i + 1));
            if (geneticCode.is_stop_codon(codon))
                return 4;
        }

        seqan::Infix<seqan::DnaString>::Type codon = seqan::infix(subSeq, 3 * (aaLength - 1), 3 * aaLength);
        if (!geneticCode.is_stop_codon(codon)) {
            return 5;
        }
        return 0;
    }

public:

    MyCDS() {}

    MyCDS(seqan::GffRecord const gff,
          seqan::DnaString const seq,
          GeneticCode const geneticCode) {
        this->gff = gff;
        this->label = classify(seq, geneticCode);
        if (is_typical()) {
            this->phase = get_start() % 3;
        }
    }

    int get_start() const {
        return gff.beginPos;
    }

    int get_end() const {
        return gff.endPos;
    }

    int get_phase() const {
        return phase;
    }

    bool is_typical() const {
        return (label == 0) ? true : false;
    }

    bool is_forward() const {
        return (gff.strand == '+') ? true : false;
    }
};

void get_myCDS_vecvec(std::vector< std::vector<MyCDS> > &myCDS_vecvec,
                      seqan::String<seqan::GffRecord> const &gffs,
                      seqan::StringSet<seqan::DnaString> const &seqs,
                      seqan::String<seqan::CharString> const &seqIds,
                      GeneticCode const geneticCode) {

    int numSeqs = seqan::length(seqs);
    for (int seqIdx = 0; seqIdx < numSeqs; seqIdx++) {
        std::vector<MyCDS> myCDS_vec;
        myCDS_vecvec.push_back(myCDS_vec);
    }

    int numGffs = seqan::length(gffs);
    for (int gffIdx = 0; gffIdx < numGffs; gffIdx++) {
        seqan::GffRecord const gff = gffs[gffIdx];
        if (gff.type == "CDS") {
            int seqIdx = -1;
            for (int i = 0; i < numSeqs; i++) {
                if (gff.ref == seqIds[i]) {
                    seqIdx = i;
                    break;
                }
            }

            if (seqIdx != -1) {
                MyCDS myCDS = MyCDS(gff, seqs[seqIdx], geneticCode);
                myCDS_vecvec[seqIdx].push_back(myCDS);
            }
        }
    }
}

//============================================================
//structure to convey shuffle command
//============================================================
struct ShuffleRegion {
    int seqIdx = -1;
    int start = -1;
    int end = -1;
    int shuffleMode = -1;
    bool isForward = true; //only required for shuffleMode = 2,3
    int phase = -1; //only required for shuffleMode = 2, 3
};

class ShuffleRegionFactory {
    bool inProcess = false;
    int cdsIdxInProcess = -1;
    std::vector<ShuffleRegion> shuffleRegions;

public:
    ShuffleRegionFactory() {}

    void initIntergenic(int seqIdx, int start, int shuffleMode) {
        assert(!inProcess);
        inProcess = true;
        ShuffleRegion shuffleRegion;
        shuffleRegion.seqIdx = seqIdx;
        shuffleRegion.start = start;
        shuffleRegion.shuffleMode = shuffleMode;
        shuffleRegions.push_back(shuffleRegion);
    }

    void completeIntergenic(int seqIdx, int end, int shuffleMode) {
        assert(inProcess);
        //some additional assertions

        if (end > shuffleRegions[shuffleRegions.size() - 1].start) {
            shuffleRegions[shuffleRegions.size() - 1].end = end;
        } else {
            shuffleRegions.pop_back();
        }
        inProcess = false;
    }

    void initGenic(int seqIdx, int tmpStart, MyCDS myCDS, int shuffleMode,
                   int cdsIdx) { //tmpStart because we need to adjust to phase
        assert(!inProcess);
        inProcess = true;

        ShuffleRegion shuffleRegion;
        shuffleRegion.seqIdx = seqIdx;
        shuffleRegion.isForward = myCDS.is_forward();
        shuffleRegion.shuffleMode = shuffleMode;

        int start;
        if (tmpStart % 3 == myCDS.get_phase()) {
            start = tmpStart;
        } else if (tmpStart % 3 < myCDS.get_phase()) {
            start = tmpStart + (myCDS.get_phase() - tmpStart % 3);
        } else if (tmpStart % 3 > myCDS.get_phase()) {
            start = tmpStart + (3 + myCDS.get_phase() - tmpStart % 3);
        }
        assert(start % 3 == myCDS.get_phase());
        shuffleRegion.start = start;

        shuffleRegions.push_back(shuffleRegion);
        cdsIdxInProcess = cdsIdx;
    }

    void completeGenic(int seqIdx, int tmpEnd, MyCDS myCDS, int shuffleMode, int cdsIdx) {
        assert(inProcess);
        assert(cdsIdx == cdsIdxInProcess);

        int end;
        if (tmpEnd % 3 == myCDS.get_phase()) {
            end = tmpEnd;
        } else if (tmpEnd % 3 < myCDS.get_phase()) {
            end = tmpEnd - (3 + tmpEnd % 3 - myCDS.get_phase());
        } else if (tmpEnd % 3 > myCDS.get_phase()) {
            end = tmpEnd - (tmpEnd % 3 - myCDS.get_phase());
        }
        assert(end % 3 == myCDS.get_phase());

        if (end > shuffleRegions[shuffleRegions.size() - 1].start) {
            shuffleRegions[shuffleRegions.size() - 1].end = end;
        } else {
            shuffleRegions.pop_back();
        }
        cdsIdxInProcess = -1;
        inProcess = false;
    }

    bool in_process() {
        return inProcess;
    }

    std::vector<ShuffleRegion> get_shuffleRegions() {
        return shuffleRegions;
    }
};

void update_genetic_code(GeneticCode &geneticCode,
                         std::vector<std::vector<MyCDS> > const myCDS_vecvec,
                         seqan::StringSet<seqan::DnaString> const &seqs) {

    int numSeqs = seqan::length(seqs);
    for (int seqIdx = 0; seqIdx < numSeqs; seqIdx++) {
        seqan::DnaString const seq = seqs[seqIdx];

        for (MyCDS const myCDS : myCDS_vecvec[seqIdx]) {
            if (myCDS.is_typical()) {
                seqan::DnaString sub = seqan::infix(seq, myCDS.get_start() + 3, myCDS.get_end() - 3);
                if (!myCDS.is_forward()) {
                    seqan::reverseComplement(sub);
                }
                geneticCode.update_count(sub);
            }
        }
    }
}

struct StateSwitch {
    int cdsIdx;
    int position;
    bool isStart;

    StateSwitch(int cdsIdx, int position, bool isStart) {
        this->cdsIdx = cdsIdx;
        this->position = position;
        this->isStart = isStart;
    }

    friend bool operator<(StateSwitch const &ss1, StateSwitch const &ss2) {
        if (ss1.position < ss2.position) {
            return true;
        } else if (ss1.position == ss2.position){
            if (!ss1.isStart && ss2.isStart){
                return true;
            } else if (ss1.isStart && !ss2.isStart){
                return false;
            } else { //need to be deterministic
                return ss1.cdsIdx < ss2.cdsIdx;
            }
        } else {
            return false;
        }
    }
};

void get_shuffle_region(std::vector<ShuffleRegion> &shuffleRegions,
                        int *shuffleMode,
                        seqan::StringSet<seqan::DnaString> const &seqs, //only needed for their length. Passing seqs is a little too much.
                        std::vector<std::vector<MyCDS> > const myCDS_vecvec) {


    shuffleRegions.clear();
    int numSeqs = seqan::length(seqs);
    for (int seqIdx = 0; seqIdx < numSeqs; seqIdx++) {
        int seqLength = seqan::length(seqs[seqIdx]);

        std::vector<StateSwitch> stateSwitchs;
        int numCDSs = myCDS_vecvec[seqIdx].size();
        for (int cdsIdx = 0; cdsIdx < numCDSs; cdsIdx++) {
            MyCDS const myCDS = myCDS_vecvec[seqIdx][cdsIdx];
            StateSwitch startSwitch = StateSwitch(cdsIdx, myCDS.get_start(), true);
            StateSwitch endSwitch = StateSwitch(cdsIdx, myCDS.get_end(), false);
            stateSwitchs.push_back(startSwitch);
            stateSwitchs.push_back(endSwitch);
        }

        std::sort(stateSwitchs.begin(), stateSwitchs.end());

        ShuffleRegionFactory shuffleRegionFactory;
        std::vector<int> inProcess; //list cdsIdx

        if (shuffleMode[0] == 1) {
            shuffleRegionFactory.initIntergenic(seqIdx, 0, 1);
        }
        for (StateSwitch stateSwitch : stateSwitchs) {
            if (stateSwitch.isStart) {
                if (inProcess.size() == 0) { //0 -> 1
                    if (shuffleMode[0] == 1) {
                        shuffleRegionFactory.completeIntergenic(seqIdx, stateSwitch.position, shuffleMode[0]);
                    }

                    MyCDS const myCDS = myCDS_vecvec[seqIdx][stateSwitch.cdsIdx];
                    if (myCDS.is_typical()) {
                        shuffleRegionFactory.initGenic(seqIdx, stateSwitch.position + 3, myCDS, shuffleMode[1],
                                                       stateSwitch.cdsIdx);
                    }
                } else if (inProcess.size() == 1) { //1 -> 2
                    MyCDS const myCDS = myCDS_vecvec[seqIdx][inProcess[0]];
                    if (myCDS.is_typical()) {
                        shuffleRegionFactory.completeGenic(seqIdx, stateSwitch.position, myCDS, shuffleMode[1],
                                                           inProcess[0]);
                    }
                }
                inProcess.push_back(stateSwitch.cdsIdx);
            } else {
                std::vector<int>::iterator position = std::find(inProcess.begin(), inProcess.end(), stateSwitch.cdsIdx);
                assert (position != inProcess.end()); // startSwitch should have come before
                inProcess.erase(position);

                if (inProcess.size() == 1) { // 2 -> 1
                    MyCDS const myCDS = myCDS_vecvec[seqIdx][inProcess[0]];
                    if (myCDS.is_typical()) {
                        shuffleRegionFactory.initGenic(seqIdx, stateSwitch.position, myCDS, shuffleMode[1],
                                                       inProcess[0]);
                    }
                } else if (inProcess.size() == 0) { // 1 -> 0
                    MyCDS const myCDS = myCDS_vecvec[seqIdx][stateSwitch.cdsIdx];
                    if (myCDS.is_typical()) {
                        shuffleRegionFactory.completeGenic(seqIdx, stateSwitch.position - 3, myCDS, shuffleMode[1],
                                                           stateSwitch.cdsIdx);
                    }
                    if (shuffleMode[0] == 1) {
                        shuffleRegionFactory.initIntergenic(seqIdx, stateSwitch.position, shuffleMode[0]);
                    }
                }
            }
        }

        if (shuffleRegionFactory.in_process()) {
            shuffleRegionFactory.completeIntergenic(seqIdx, seqLength, 1);
        }

        std::vector<ShuffleRegion> shuffleRegionsForSeq = shuffleRegionFactory.get_shuffleRegions();
        shuffleRegions.insert(shuffleRegions.end(), shuffleRegionsForSeq.begin(), shuffleRegionsForSeq.end());
    }

    return;
}

void output_shuffle_regions(std::vector<ShuffleRegion> const &shuffleRegions,
                            seqan::String<seqan::CharString> const &seqIds,
                            std::string const &outFilepath) {

    std::ofstream ofs(outFilepath);
    if (!ofs) {
        std::cerr << "ERROR: Could not open " << outFilepath << std::endl;
    }

    for(ShuffleRegion const &shuffleRegion : shuffleRegions){
        ofs<< seqIds[shuffleRegion.seqIdx]<<"\t"
           << shuffleRegion.start << "\t"
           << shuffleRegion.end << "\t"
           << shuffleRegion.shuffleMode << "\t"
           << 0 << "\t"
           << (shuffleRegion.isForward ? '+' : '-') << endl;
    }
    ofs.close();
    return;
}

void shuffle_region(seqan::DnaString & seq,
                    ShuffleRegion const shuffleRegion,
                    std::mt19937 & mt,
                    GeneticCode const & geneticCode) {
    seqan::DnaString subSeq = seqan::infix(seq, shuffleRegion.start, shuffleRegion.end);
    if(!shuffleRegion.isForward){
        seqan::reverseComplement(subSeq);
    }

    switch(shuffleRegion.shuffleMode) {
        case 1:
            shuffle_base(subSeq, mt);
            break;
        case 2:
            shuffle_codon(subSeq, mt);
            break;
        case 3:
            shuffle_synonymous(subSeq, geneticCode, mt);
            break;
    }

    if(!shuffleRegion.isForward){
        seqan::reverseComplement(subSeq);
    }
    seqan::replace(seq, shuffleRegion.start, shuffleRegion.end, subSeq);
    return;
}

void shuffle_genome(seqan::StringSet<seqan::DnaString> &seqs,
                    std::vector<ShuffleRegion> const &shuffleRegions,
                    GeneticCode const &geneticCode) {

    std::random_device rd;
    std::mt19937 mt(rd());

    for (ShuffleRegion shuffleRegion : shuffleRegions) {
        std::cout<<shuffleRegion.start<<","<<shuffleRegion.end<<endl;
        shuffle_region(seqs[shuffleRegion.seqIdx], shuffleRegion, mt, geneticCode);
    }
    return;
}

#endif //SHUFFLE_GENOME_SHUFFLE_GENOME_H
