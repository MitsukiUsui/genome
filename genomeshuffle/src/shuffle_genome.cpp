#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <random>

#define SEQAN_ENABLE_DEBUG 1

#include <seqan/arg_parse.h>
#include "mygenome.h"
#include "myutil.h"

//#include "/Users/mitsuki/sandbox/genome/genomeshuffle/src/mygenome.h"
//#include "/Users/mitsuki/sandbox/genome/genomeshuffle/src/myutil.h"


using std::cout;
using std::cerr;
using std::endl;


template <typename TSeq>
void shuffle_base(TSeq & seq, std::mt19937 & mt){
	unsigned bpLength=seqan::length(seq);
	std::vector<unsigned> v(bpLength);
	std::iota(v.begin(), v.end(), 0);
	std::shuffle(v.begin(), v.end(), mt);

	seqan::DnaString tmp;//TBI: adjust to template
	seqan::resize(tmp, bpLength);
	for(unsigned i=0; i < bpLength; i++){
		tmp[i]=seq[v[i]];
	}
	seqan::replace(seq, 0, bpLength, tmp);
	return;
}

template <typename TSeq>
void shuffle_codon(TSeq & seq, std::mt19937 & mt){
	unsigned bpLength=seqan::length(seq);
	assert(bpLength % 3 == 0);
	unsigned aaLength = bpLength / 3;
	std::vector<unsigned> v(aaLength);
	std::iota(v.begin(), v.end(), 0);
	std::shuffle(v.begin(), v.end(), mt);

	seqan::DnaString tmp;//TBI: adjust to template
	seqan::resize(tmp, bpLength);
	for(unsigned i=0;i<aaLength;i++){
		for(unsigned j=0;j<3;j++){
			tmp[3*i+j]=seq[3*v[i]+j];
		}
	}
	seqan::replace(seq, 0, bpLength, tmp);
	return;
}

template <typename TSeq>
void shuffle_synonymous(TSeq & seq, GeneticCode const & gc, std::mt19937 & mt){
	unsigned bpLength=seqan::length(seq);
	assert(bpLength % 3 == 0);
	unsigned aaLength = bpLength / 3;

	for(unsigned i=0; i<aaLength; i++){
		seqan::Infix<seqan::DnaString>::Type codon=seqan::infix(seq, 3*i, 3*(i+1));
		gc.synonymous_sub(codon, mt);
	}
	return;
}



//============================================================
//structure to convey shuffle command
//============================================================
struct ShuffleRegion{/*{{{*/
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
};/*}}}*/

template <typename TSeq>
void overwrite_region(TSeq & seq, ShuffleRegion const & sr, GeneticCode const & gc, std::mt19937 & mt){
	seqan::DnaString subseq=seqan::infix(seq, sr.start, sr.end);//create new DnaString object
	if(!sr.isForward){
		seqan::reverseComplement(subseq);
	}
	
	switch (sr.mode){
		case 1:
			shuffle_base(subseq, mt);
			break;
		case 2:
			shuffle_codon(subseq, mt);
			break;
		case 3:
			shuffle_synonymous(subseq, gc, mt);
			break;
		default:
			cerr<<"ERROR: undefined shuffle mode "<<sr.mode<<endl;
			std::exit(1);
	}

	if(!sr.isForward){
		seqan::reverseComplement(subseq);
	}

	//overwrite
	seqan::replace(seq, sr.start, sr.end, subseq);
	return;
}



//============================================================
//shuffle genome according to cdss end report
//according to cdss and shuffleMode, construct SuffleRegion instance and call overwrite_region() method
//============================================================
template <typename TSeqs>
void shuffle_genome(TSeqs & seqs, CDSs const  & cdss, GeneticCode const & gc, int * shuffleMode){
	std::random_device rd;
	std::mt19937 mt(rd());

	std::vector<int> shufflableTypes{0}; 

	cout<<"\tSHUFFLE REPORT : "<<endl;
	for(unsigned refIdx=0; refIdx<seqan::length(seqs);refIdx++){
		int countShuffle[4] = { 0, 0, 0, 0 };
		
		int endMax = 0;
		for (auto itr=cdss.cdss_vec[refIdx].begin(); itr!=cdss.cdss_vec[refIdx].end()-1; itr++){//skip last sentinel
			int		 start = itr->startPos;	
			int		   end = itr->endPos;
			int		  type = itr->type;
			bool isForward = itr->isForward;
			int  startNext = (itr+1)->startPos; 
			
			assert(start<=startNext);//regions should be sorted
			//shuffle intergenic region
			if(shuffleMode[0] == 1){
				if (start > endMax){
					int shuffleStart=endMax;
					int shuffleEnd=start;
					ShuffleRegion sr(shuffleStart, shuffleEnd, true, 1);
					overwrite_region(seqs[refIdx], sr, gc, mt);
					countShuffle[1]+=(sr.end-sr.start);
				}
			}
			//shuffle genetic region
			if (shuffleMode[1] > 0){
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
	
		//calculate unshuffled region
		countShuffle[0] = seqan::length(seqs[refIdx]);
		for (int i = 1; i < 4; i++){
			countShuffle[0] -= countShuffle[i];
		}
		//shuffle report
		cout<<"\t\tREFIDX: "<<refIdx<<endl;
		cout<<"\t\t\tTotal: "<<seqan::length(seqs[refIdx])<<" bp."<<endl;
		for (int i = 0; i < 4; i++){
			cout << "\t\t\t\tmode" << i << " shuffled: " << countShuffle[i] << " bp." << endl;
		}
	}
}

//------------------------------------------------------------
//update GeneticCode (codonCount) according to CDSs
//------------------------------------------------------------
void update_gc(GeneticCode & gc, seqan::StringSet<seqan::Dna5String> & seqs, CDSs const & cdss){//TBI: utilize template/*{{{*/
	for(int refIdx=0;refIdx<seqan::length(seqs);refIdx++){
		for(auto itr=cdss.cdss_vec[refIdx].begin(); itr!=cdss.cdss_vec[refIdx].end()-1;itr++){
			int		 start = itr->startPos;	
			int		   end = itr->endPos;
			int		  type = itr->type;
			bool isForward = itr->isForward;
			
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


//------------------------------------------------------------
//parse arguments 
//------------------------------------------------------------
void set_parser(seqan::ArgumentParser & parser, int argc, char ** argv){/*{{{*/
	seqan::addUsageLine(parser,
		"\"seqFilepath\" \"gffFilepath\"");
	seqan::addDescription(parser,
		"This program . "
		);
	addArgument(parser, seqan::ArgParseArgument(
		seqan::ArgParseArgument::INPUT_FILE, "seqFilepath"));
	addArgument(parser, seqan::ArgParseArgument(
		seqan::ArgParseArgument::INPUT_FILE, "gffFilepath"));
	addArgument(parser, seqan::ArgParseArgument(
		seqan::ArgParseArgument::OUTPUT_FILE, "outFilepath"));
	addArgument(parser, seqan::ArgParseArgument(
		seqan::ArgParseArgument::INTEGER, "shuffleMode1"));
	addArgument(parser, seqan::ArgParseArgument(
		seqan::ArgParseArgument::INTEGER, "shuffleMode2"));
	addOption(parser, seqan::ArgParseOption(
		"d", "dryrun", "Do dryrun"));
	
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
	if (res != seqan::ArgumentParser::PARSE_OK){
		std::exit(1);
	}
}/*}}}*/


int main(int argc, char ** argv){
	//------------------------------------------------------------
	//parse arguments
	//------------------------------------------------------------
	seqan::ArgumentParser parser("shuffle_genome");
	set_parser(parser, argc, argv);
	
	seqan::CharString seqFilepath;
	seqan::getArgumentValue(seqFilepath, parser, 0);
	seqan::CharString gffFilepath;
	seqan::getArgumentValue(gffFilepath, parser, 1);
	seqan::CharString outFilepath;
	seqan::getArgumentValue(outFilepath, parser, 2);
	int shuffleMode[2]={};
	seqan::getArgumentValue(shuffleMode[0], parser, 3);
	seqan::getArgumentValue(shuffleMode[1], parser, 4);
	bool dryrun=seqan::isSet(parser, "dryrun");
	
	//configuration
	cout<<"\tCONFIGURATION"<<endl;
	cout << "\t\tseqFilepath	: " << seqFilepath << endl;
	cout << "\t\tgffFilepath	: " << gffFilepath << endl;
	cout << "\t\toutFilepath	: " << outFilepath << endl;
	cout << "\t\tshuffle mode : " << shuffleMode[0] << ", " << shuffleMode[1] << endl;
	if(dryrun){
		cout<<"\t\texecution	  : DRYRUN"<<endl;
	}

	//------------------------------------------------------------
	//process fasta 
	//------------------------------------------------------------
	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::Dna5String> seqs;
	read_fasta(ids, seqs, seqFilepath);
	//trim ids
	for (unsigned i=0;i<length(ids);i++){
		ids[i]=seqan::CharString(split(seqan::toCString(ids[i]),' ')[0]);
	}
	//summarize 
	cout<<"\tDONE reading "<<seqFilepath<<endl;
	for(unsigned i=0;i<seqan::length(ids);i++){
		cout<<"\t\t"<<ids[i]<<" : "<<seqan::length(seqs[i])<<" bp"<<endl;
	}

	//------------------------------------------------------------
	//process gff and genetic code
	//------------------------------------------------------------
	std::string str="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	GeneticCode gc(str);//genetic code needs to be initialized first, in order to judge cds types
	seqan::String<seqan::GffRecord> records;
	read_gff(records,gffFilepath);
	CDSs cdss(records, ids, seqs, gc);
	//summarize 
	cout<<"\tDONE reading "<<gffFilepath<<endl;
	for(unsigned i=0;i<seqan::length(ids);i++){
		cout<<"\t\t"<<ids[i]<<" : "<<cdss.cdss_vec[i].size()<<" cds records"<<endl;
	}
	update_gc(gc, seqs, cdss);//update genetic code

	//cdss.__show();
	//cout<<endl;
	//gc.__show();
	//cout<<endl;
	//gc.__show_freq();
	//cout<<endl;
	//gc.__show_freq(false);
	//cout<<endl;
	
	shuffle_genome(seqs, cdss, gc, shuffleMode);

	//gc.clear_count();
	//gc.__show_freq(false);
	//cout<<endl;
	//update_gc(gc, seqs, cdss);
	//gc.__show_freq(false);
	//cout<<endl;
	

	write_fasta(ids, seqs, outFilepath);	
	
	return 0;
}
