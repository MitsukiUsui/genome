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

//------------------------------------------------------------
//update GeneticCode (codonCount) according to CDSs
//------------------------------------------------------------
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
	seqan::ArgumentParser parser("shuffle_report");
	set_parser(parser, argc, argv);
	
	seqan::CharString seqFilepath;
	seqan::getArgumentValue(seqFilepath, parser, 0);
	seqan::CharString gffFilepath;
	seqan::getArgumentValue(gffFilepath, parser, 1);
	int shuffleMode[2]={};
	seqan::getArgumentValue(shuffleMode[0], parser, 2);
	seqan::getArgumentValue(shuffleMode[1], parser, 3);
	bool dryrun=seqan::isSet(parser, "dryrun");
	//summarize 
	cout<<endl<<"PROCESSING..."<<endl;
	cout << "\tseqFilepath  : " << seqFilepath << endl;
    cout << "\tgffFilepath  : " << gffFilepath << endl;
    cout << "\tshuffle mode : " << shuffleMode[0] << ", " << shuffleMode[1] << endl;
	if(dryrun){
		cout<<"\texecution    : DRYRUN"<<endl;
	}
    cout << endl;
	

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
	cout<<"DONE reading "<<seqFilepath<<endl;
	for(unsigned i=0;i<seqan::length(ids);i++){
		cout<<"\t"<<ids[i]<<" : "<<seqan::length(seqs[i])<<endl;
	}
	cout<<endl;

	//------------------------------------------------------------
	//process gff and genetic code
	//------------------------------------------------------------
	std::string str="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    GeneticCode gc(str);//genetic code needs to be initialized first, in order to judge cds types
	seqan::String<seqan::GffRecord> records;
	read_gff(records,gffFilepath);
	CDSs cdss(records, ids, seqs, gc);
	//summarize 
	cout<<"DONE reading "<<gffFilepath<<endl;
	for(unsigned i=0;i<seqan::length(ids);i++){
		cout<<"\t"<<ids[i]<<" : "<<cdss.cdss_vec[i].size()<<endl;
	}
	cout<<endl;
	update_gc(gc, seqs, cdss);//update genetic code

	cdss.__show();
	cout<<endl;
    gc.__show();
	cout<<endl;
	gc.__show_freq();
	cout<<endl;
	gc.__show_freq(false);
	cout<<endl;
    
	shuffle_genome(seqs, cdss, gc, shuffleMode);

	return 0;
}
