#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <random>

#define SEQAN_ENABLE_DEBUG 1

#include "mygenome.h"
#include "myutil.h"

//#include "/Users/mitsuki/sandbox/genome/genomeshuffle/src/mygenome.h"
//#include "/Users/mitsuki/sandbox/genome/genomeshuffle/src/myutil.h"


using std::cout;
using std::cerr;
using std::endl;

seqan::Dna5String shuffle_base(seqan::Dna5String & seq){
	seqan::Dna5String retVal;
	seqan::resize(retVal, seqan::length(seq));
	for(unsigned i=0;i<seqan::length(seq);i++){
		retVal[i]='A';
	}
	return retVal;
}

seqan::Dna5String shuffle_codon(seqan::Dna5String & seq){
	seqan::Dna5String retVal;
	seqan::resize(retVal, seqan::length(seq));
	for(unsigned i=0;i<seqan::length(seq);i++){
		retVal[i]='C';
	}
	return retVal;
}

seqan::Dna5String shuffle_synonymous(seqan::Dna5String & seq){
	seqan::Dna5String retVal;
	seqan::resize(retVal, seqan::length(seq));
	for(unsigned i=0;i<seqan::length(seq);i++){
		retVal[i]='G';
	}
	return retVal;
}


struct ShuffleRegion{
	int start;
	int end;
	bool isForward;
	int mode;

	ShuffleRegion(int s, int e, bool i, int m){
		start=s;
		end=e;
		isForward=i;
		mode=m;
	}
};

template <typename TSeq>
void overwrite_region(TSeq & seq, ShuffleRegion const & region, GeneticCode & gc, std::mt19937 & mt){
	//do overwrite
}




//shuffle genome according to cdss end report
template <typename TSeqs>
void shuffle_genome(TSeqs & seqs, CDSs const  & cdss, GeneticCode & gc, int * shuffleMode){
	std::random_device rd;
	std::mt19937 mt(rd());

	std::vector<int> shufflableTypes{0}; 

	for(unsigned refIdx=0; refIdx<seqan::length(seqs);refIdx++){
		int countShuffle[4] = { 0, 0, 0, 0 };
		
		int endMax = 0;
		for(int i=0;i<cdss.cdss_vec[refIdx].size()-1;i++){
			
			int      start = cdss.cdss_vec[refIdx][i].startPos;	
			int        end = cdss.cdss_vec[refIdx][i].endPos;
			int       type = cdss.cdss_vec[refIdx][i].type;
			bool isForward = cdss.cdss_vec[refIdx][i].isForward;
			int startNext = cdss.cdss_vec[refIdx][i+1].startPos;	
			
			
			if(refIdx==2){
				cout<<start<<","<<startNext<<endl;
			}
			
			//if(start>startNext){
			//	cout<<start<<","<<startNext<<endl;
			//	std::exit(1);	
			//}

			assert(start<=startNext);//regions should be sorted

			if(shuffleMode[0] == 1){//shuffle intergenic region
				if (start > endMax){
					ShuffleRegion sr(endMax, start, true, 1);
					overwrite_region(seqs[refIdx], sr, gc, mt);
					countShuffle[1]+=(sr.end-sr.start);
				}
			}
			if (shuffleMode[1] > 0){//shuffle genetic region
				if (is_in(type, shufflableTypes)){
					int cand1 = start + 3 * ceil((double)(endMax - start) / 3);
					int cand2 = start + 3;
					int shuffleStart = std::max(cand1, cand2);
					
					cand1 = end - 3 * ceil((double)(end - startNext) / 3);
					cand2 = end - 3;
					int shuffleEnd = std::min(cand1, cand2);



					if(shuffleEnd > shuffleStart){
						ShuffleRegion sr(shuffleStart, shuffleEnd, isForward, shuffleMode[1]);
						overwrite_region(seqs[refIdx], sr, gc, mt);
						countShuffle[shuffleMode[1]]+=(sr.end-sr.start);
					}
				}
			}
			endMax=std::max(end, endMax);
		}

		countShuffle[0] = seqan::length(seqs[refIdx]);
		for (int i = 1; i < 4; i++){
			countShuffle[0] -= countShuffle[i];
		}

		cout<<"REFIDX: "<<refIdx<<endl;
		cout<<"\tTotal: "<<seqan::length(seqs[refIdx])<<" bp."<<endl;
		for (int i = 0; i < 4; i++){
			cout << "\t\tmode" << i << " shuffled: " << countShuffle[i] << " bp." << endl;
		}
		cout<<endl;
	}
}


//update codon count according to CDSs
void update_gc(GeneticCode & gc, seqan::StringSet<seqan::Dna5String> & seqs, CDSs const & cdss){/*{{{*/
	for(int refIdx=0;refIdx<seqan::length(seqs);refIdx++){
		for(int i=0;i<cdss.cdss_vec[refIdx].size()-1;i++){
			int      start = cdss.cdss_vec[refIdx][i].startPos;	
			int        end = cdss.cdss_vec[refIdx][i].endPos;
			int       type = cdss.cdss_vec[refIdx][i].type;
			bool isForward = cdss.cdss_vec[refIdx][i].isForward;

			if(type==0){
				if(isForward){
					seqan::Infix<seqan::Dna5String>::Type sub=seqan::infix(seqs[refIdx],start, end);
					gc.update_count(sub);
				}
				else{
					seqan::Dna5String sub=seqan::infix(seqs[refIdx],start,end);
					seqan::reverseComplement(sub);
					gc.update_count(sub);
				}
			}
		}	
	}
	gc.calc_freq();
}/*}}}*/


int main(int argc, char *argv[]){
	//read fasta file
	seqan::CharString seqFilepath="GCF_000022205.1_ASM2220v1_genomic.fna";
	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::Dna5String> seqs;
	read_fasta(ids, seqs, seqFilepath);
	
	//trim ids
	for (unsigned i=0;i<length(ids);i++){
		ids[i]=seqan::CharString(split(seqan::toCString(ids[i]),' ')[0]);
	}

	//read gff file
	seqan::CharString gffFilepath="GCF_000022205.1_ASM2220v1_genomic.gff";
	seqan::String<seqan::GffRecord> records;
	read_gff(records,gffFilepath);

	//initialize genetic code
	std::string str="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    GeneticCode gc(str);
	
	//create CDSs
	CDSs cdss(records, ids, seqs, gc);

	//update genetic code
	update_gc(gc, seqs, cdss);

	int shuffleMode[2]={1,3};
    
	//output property
	cout << "seqFilepath : " << seqFilepath << endl;
    cout << "gffFilepath : " << gffFilepath << endl;
    cout << "shuffle mode : " << shuffleMode[0] << ", " << shuffleMode[1] << endl;
    cout << endl;
	for(unsigned i=0;i<seqan::length(ids);i++){
		cout<<ids[i]<<" : "<<seqan::length(seqs[i])<<endl;
	}
	cout<<endl;

	cdss.__show();
	cout<<endl;
    //gc.__show();
	//cout<<endl;
	gc.__show_freq();
	cout<<endl;
	gc.__show_freq(false);
	cout<<endl;
    
	shuffle_genome(seqs, cdss, gc, shuffleMode);

	return 0;
}
