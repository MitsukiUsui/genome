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

#define SEQAN_ENABLE_DEBUG 1

#include <seqan/arg_parse.h>
#include "myseqan.h"
#include "myutil.h"

#include "geneticcode.h"


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

template<typename TSeq>
void shuffle_synonymous(TSeq &seq, GeneticCode const &gc, std::mt19937 &mt) {
    unsigned bpLength = seqan::length(seq);
    assert(bpLength % 3 == 0);
    unsigned aaLength = bpLength / 3;

    for (unsigned i = 0; i < aaLength; i++) {
        seqan::Infix<seqan::DnaString>::Type codon = seqan::infix(seq, 3 * i, 3 * (i + 1));
        gc.synonymous_sub(codon, mt);
    }
    return;
}


//============================================================
//structure to convey shuffle command
//============================================================
struct ShuffleRegion {
    seqan::CharString seqId;
    int start;
    int end;
    bool isForward;
    int shuffleMode;

    ShuffleRegion(std::string seqId, int start, int end, bool isForward, int shuffleMode) {
        this->seqId = seqId;
        this->start = start;
        this->end = end;
        this->isForward = isForward;
        this->shuffleMode = shuffleMode;
    }
};

//template <typename TSeq>
//void overwrite_region(TSeq & seq, ShuffleRegion const & sr, GeneticCode const & gc, std::mt19937 & mt){
//	seqan::DnaString subseq=seqan::infix(seq, sr.start, sr.end);//create new DnaString object
//	if(!sr.isForward){
//		seqan::reverseComplement(subseq);
//	}
//
//	switch (sr.mode){
//		case 1:
//			shuffle_base(subseq, mt);
//			break;
//		case 2:
//			shuffle_codon(subseq, mt);
//			break;
//		case 3:
//			shuffle_synonymous(subseq, gc, mt);
//			break;
//		default:
//			cerr<<"ERROR: undefined shuffle mode "<<sr.mode<<endl;
//			std::exit(1);
//	}
//
//	if(!sr.isForward){
//		seqan::reverseComplement(subseq);
//	}
//
//	//overwrite
//	seqan::replace(seq, sr.start, sr.end, subseq);
//	return;
//}



//============================================================
//shuffle genome according to cdss end report
//according to cdss and shuffleMode, construct SuffleRegion instance and call overwrite_region() method
//============================================================
//template <typename TSeqs>
//void shuffle_genome(TSeqs & seqs, CDSs const  & cdss, GeneticCode const & gc, int * shuffleMode){
//	std::random_device rd;
//	std::mt19937 mt(rd());
//
//	std::vector<int> shufflableTypes{0};
//
//	cout<<"\tSHUFFLE REPORT : "<<endl;
//	for(unsigned refIdx=0; refIdx<seqan::length(seqs);refIdx++){
//		int countShuffle[4] = { 0, 0, 0, 0 };
//
//		int endMax = 0;
//		for (auto itr=cdss.cdss_vec[refIdx].begin(); itr!=cdss.cdss_vec[refIdx].end()-1; itr++){//skip last sentinel
//			int		 start = itr->startPos;
//			int		   end = itr->endPos;
//			int		  type = itr->type;
//			bool isForward = itr->isForward;
//			int  startNext = (itr+1)->startPos;
//
//			assert(start<=startNext);//regions should be sorted
//			//shuffle intergenic region
//			if(shuffleMode[0] == 1){
//				if (start > endMax){
//					int shuffleStart=endMax;
//					int shuffleEnd=start;
//					ShuffleRegion sr(shuffleStart, shuffleEnd, true, 1);
//					overwrite_region(seqs[refIdx], sr, gc, mt);
//					countShuffle[1]+=(sr.end-sr.start);
//				}
//			}
//			//shuffle genetic region
//			if (shuffleMode[1] > 0){
//				if (is_in(type, shufflableTypes)){
//					int cand1 = start + 3 * ceil((double)(endMax - start) / 3);
//					int cand2 = start + 3;
//					int shuffleStart = std::max(cand1, cand2);
//
//					cand1 = end - 3 * ceil((double)(end - startNext) / 3);
//					cand2 = end - 3;
//					int shuffleEnd = std::min(cand1, cand2);
//
//					if(shuffleEnd > shuffleStart){
//						ShuffleRegion sr(shuffleStart, shuffleEnd, isForward, shuffleMode[1]);
//						overwrite_region(seqs[refIdx], sr, gc, mt);
//						countShuffle[shuffleMode[1]]+=(sr.end-sr.start);
//					}
//				}
//			}
//			endMax=std::max(end, endMax);
//		}
//
//		//calculate unshuffled region
//		countShuffle[0] = seqan::length(seqs[refIdx]);
//		for (int i = 1; i < 4; i++){
//			countShuffle[0] -= countShuffle[i];
//		}
//		//shuffle report
//		cout<<"\t\tREFIDX: "<<refIdx<<endl;
//		cout<<"\t\t\tTotal: "<<seqan::length(seqs[refIdx])<<" bp."<<endl;
//		for (int i = 0; i < 4; i++){
//			cout << "\t\t\t\tmode" << i << " shuffled: " << countShuffle[i] << " bp." << endl;
//		}
//	}
//}


void update_genetic_code(GeneticCode &geneticCode,
                         seqan::String<seqan::CharString> const &seqIds,
                         seqan::StringSet<seqan::DnaString> const &seqs,
                         std::vector<int> const &gffLabels,
                         seqan::String<seqan::GffRecord> const &gffs) {

    assert (seqan::length(seqIds) == seqan::length(seqs));
    int numSeqs = seqan::length(seqIds);
    for (int i = 0; i < numSeqs; i++) {
        seqan::CharString const seqId = seqIds[i];
        seqan::DnaString const seq = seqs[i];

        assert (seqan::length(gffLabels) == seqan::length(gffs));
        int numGffs = seqan::length(gffs);

        for (int j = 0; j < numGffs; j++) {
            seqan::GffRecord const gff = gffs[j];
            int gffLabel = gffLabels[j];

            if (gff.ref == seqId && gffLabel == 0) {
                int start = gff.beginPos;
                int end = gff.endPos;
                bool isForward = gff.phase == '+' ? true : false;

                seqan::DnaString sub = seqan::infix(seq, start, end);
                if (isForward) {
                    geneticCode.update_count(sub);
                } else {
                    seqan::reverseComplement(sub);
                    geneticCode.update_count(sub);
                }
            }
        }
    }
    return;
}

void classify_gff(std::vector<int> gffLabels,
                  seqan::String<seqan::GffRecord> const &gffs,
                  seqan::String<seqan::CharString> const &seqIds,
                  seqan::StringSet<seqan::DnaString> const &seqs,
                  GeneticCode const &geneticCode) {

    gffLabels.resize(seqan::length(gffs));
    for (int label : gffLabels) {
        label = -1;
    }
    return;
}

void get_shuffle_region(std::vector<ShuffleRegion> &shuffleRegions,
                        int *shuffleMode,
                        seqan::String<seqan::CharString> const &seqIds,
                        seqan::StringSet<seqan::DnaString> const &seqs,
                        std::vector<int> const &gffLabels,
                        seqan::String<seqan::GffRecord> const &gffs) {

    assert (seqan::length(seqIds) == seqan::length(seqs));
    int numSeqs = seqan::length(seqIds);
    for (int i = 0; i < numSeqs; i++) {
        seqan::CharString const seqId = seqIds[i];
        seqan::DnaString const seq = seqs[i];

        assert (seqan::length(gffLabels) == seqan::length(gffs));
        int numGffs = seqan::length(gffs);

        //create list first
        for (int j = 0; j < numGffs; j++) {
            seqan::GffRecord const gff = gffs[j];
            int gffLabel = gffLabels[j];

            if (gff.ref == seqId && gffLabel == 0) {
                //WRITE ME
            }
        }

        //sort

        //push_back
    }

    shuffleRegions.clear();
    ShuffleRegion shuffleRegion("test", 0, 0, true, 0);
    shuffleRegions.push_back(shuffleRegion);
    return;
}

void shuffle_region(seqan::DnaString seq, ShuffleRegion shuffleRegion) {
    //TODO WIRTE ME
    return;
}

void shuffle_genome(std::vector<ShuffleRegion> const &shuffleRegions,
                    seqan::String<seqan::CharString> const &seqIds,
                    seqan::StringSet<seqan::DnaString> const &seqs) {

    assert (seqan::length(seqIds) == seqan::length(seqs));
    int numSeqs = seqan::length(seqIds);
    for (int i = 0; i < numSeqs; i++) {
        seqan::CharString seqId = seqIds[i];
        seqan::DnaString seq = seqs[i];

        for (ShuffleRegion shuffleRegion : shuffleRegions) {
            if (shuffleRegion.seqId == seqId) {
                shuffle_region(seq, shuffleRegion);
            }
        }
    }
    return;
}

#endif //SHUFFLE_GENOME_SHUFFLE_GENOME_H
