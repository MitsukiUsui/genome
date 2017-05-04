#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <random>

#define SEQAN_ENABLE_DEBUG 1

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include <../../genomeshuffle/src/mygenome.h>
#include <../../genomeshuffle/src/myutil.h>

using std::cout;
using std::cerr;
using std::endl;


//------------------------------------------------------------
//Define labels (0: chromosome, 1:plasmid) according to each ids(heders)
//------------------------------------------------------------
void judge_type(std::vector<int> & labels, seqan::StringSet<seqan::CharString> const & ids){/*{{{*/
    labels.resize(seqan::length(ids));
	for(unsigned i=0;i<seqan::length(ids);i++){
        std::string id_str=seqan::toCString(ids[i]);
        std::string::size_type index=id_str.find("plasmid");
        if(index==std::string::npos){
            labels[i]=0;
        }
        else{
			labels[i]=1;
        }
    }
    return;
}/*}}}*/


//------------------------------------------------------------
//fill N randomly
//------------------------------------------------------------
void fill_N(seqan::Dna5String & seq, std::mt19937 & mt){/*{{{*/
	char aa[4]={'T','C','A','G'};
	double aaFreq[4]={0.25,0.25,0.25,0.25};

	for(unsigned i=0;i<seqan::length(seq);i++){
		if(seq[i]=='N'){
			std::uniform_real_distribution<double> dist(0.0,1.0);
			double d=dist(mt);	
			double cumSum=0;
			for(int idx = 0;idx < 4;idx++){
				cumSum+=aaFreq[idx];
				if(cumSum>=d){
					seq[i]=aa[idx];
					break;
				}
			}
		}
	}
	return;
}/*}}}*/


//------------------------------------------------------------
//output fasta separately according to labels
//if convertFlag > 0, conversion into Dna5String or DnaString is carried out simultaneously
//------------------------------------------------------------
void split_fasta(seqan::StringSet<seqan::CharString> const & outFilepaths, seqan::StringSet<seqan::CharString> const & ids, seqan::StringSet<seqan::IupacString> & seqs, std::vector<int> const & labels, int convertFlag){/*{{{*/
	
	assert(seqan::length(outFilepaths)==2);
	std::random_device rd;
	std::mt19937 mt(rd());


	for(unsigned label=0;label<2;label++){
		//open file
		seqan::SeqFileOut seqFileout;	
		if(!seqan::open(seqFileout, seqan::toCString(outFilepaths[label]))){
			cerr<<"ERROR: Could not open "<<outFilepaths[label]<<endl;;
			std::exit(1);
		}
		//output according to labels & convertFlag
		for(unsigned i=0;i<seqan::length(ids);i++){
			if(labels[i]==label){
				switch (convertFlag){
					case 0:{
						seqan::writeRecord(seqFileout,ids[i],seqs[i]);
						break;
					}
					case 1:{
						seqan::Dna5String dna5=seqs[i];
						seqan::writeRecord(seqFileout,ids[i],dna5);
						break;
					}
					case 2:{
						seqan::Dna5String dna5=seqs[i];
						fill_N(dna5, mt);
						seqan::DnaString dna=dna5;
						seqan::writeRecord(seqFileout,ids[i],dna);
						break;
					}
				}
			}
		}
	}
}/*}}}*/

//------------------------------------------------------------
//initialize parser
//------------------------------------------------------------
void set_parser(seqan::ArgumentParser & parser, int argc, char ** argv){/*{{{*/
    seqan::addUsageLine(parser,
        "\"basename\" \"seqFilepath\" \"chrFilepath\" \"plsFilepath\" \"logFilepath\" [OPTIONS]");
    seqan::addDescription(parser,
        "This program split input fasta file into 2 output fasta (chromosome, plasmid). "
        "The summary of spliting is written in log file. "
        "Beside, when given convert option, convert sequences into Dna5 or Dna sequence while splitting. "
        );
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "basename"));
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "seqFilepath"));
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::OUTPUT_FILE, "chrFilepath"));
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::OUTPUT_FILE, "plsFilepath"));
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::OUTPUT_FILE, "logFilepath"));
    addOption(parser, seqan::ArgParseOption(
        "c", "convert", "Whether to convert (0:no conversion, 1:to Dna5, 2:to Dna)",
		seqan::ArgParseArgument::INTEGER, "INT"));
	seqan::setDefaultValue(parser, "convert", "0");
	
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK){
    	std::exit(1);
	}
}/*}}}*/


int main(int argc, char ** argv){
	//argument setting
	seqan::ArgumentParser parser("fasta_split");
    set_parser(parser,argc,argv);
	seqan::CharString basename;
	seqan::getArgumentValue(basename, parser, 0);
	seqan::CharString seqFilepath;
	seqan::getArgumentValue(seqFilepath, parser, 1);
	seqan::StringSet<seqan::CharString> outFilepaths;
	seqan::resize(outFilepaths, 2);
	seqan::getArgumentValue(outFilepaths[0], parser, 2);
	seqan::getArgumentValue(outFilepaths[1], parser, 3);
	seqan::CharString logFilepath;
	seqan::getArgumentValue(logFilepath, parser, 4);
	int convertFlag;
	seqan::getOptionValue(convertFlag, parser,"convert");


    //read input file
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::IupacString> seqs;
    read_fasta(ids,seqs,seqFilepath);
    cout<<"\tDONE reading "<<seqan::length(ids)<<" seqs from "<<seqFilepath<<endl;

	//classify chromosome and plasmid
	std::vector<std::string> legends;
	legends={"chromosome", "plasmid"};
	std::vector<int> labels;
    judge_type(labels,ids);

	int labelCounts[2]={};
	for(unsigned i=0;i<seqan::length(labels);i++){
		labelCounts[labels[i]]++;
	}

	//write fasta
	split_fasta(outFilepaths, ids, seqs, labels, convertFlag);
    for(unsigned label=0;label<2;label++){
		cout<<"\tDONE writing "<<labelCounts[label]<<" seqs to   "<<outFilepaths[label]<<" with conversion "<<convertFlag<<endl;
	}
	
	//output log
	std::ofstream ofs(seqan::toCString(logFilepath));
    if (!ofs){
		cerr<<"ERROR: Could not open "<<logFilepath<<endl;;
		return 1;
    }
	for(unsigned label=0;label<2;label++){
		std::stringstream ss;//aggregate length of sequences
		for(unsigned i=0;i<seqan::length(ids);i++){
			if(labels[i]==label){
				ss<<seqan::length(seqs[i])<<":";
			}
		}
    	std::string length=ss.str();
		if(labelCounts[label]>0){//drop last ':'
			length.pop_back();
		}
    	ofs<<basename<<","<<legends[label]<<","<<outFilepaths[label]<<","<<labelCounts[label]<<","<<length<<endl;
	} 
	cout<<"\tDONE writing log to "<<logFilepath<<endl;
	
	return 0;
}
