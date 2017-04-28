#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>
#include <seqan/stream.h>

#include "myutil.h"

using std::cout;
using std::cerr;
using std::endl;

//0...chromosome, 1...plasmid
void judge_type(std::vector<int> & labels, seqan::StringSet<seqan::CharString> const & ids){
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
}



//do not use template because seqs are required to be able to be converted into Dna5String
void split_fasta(seqan::StringSet<seqan::CharString> const & outFilepaths, seqan::StringSet<seqan::CharString> const & ids, seqan::StringSet<seqan::IupacString> & seqs, std::vector<int> const & labels, bool doModification){
	assert(seqan::length(outFilepaths)==2);

	for(unsigned label=0;label<2;label++){
		//open file
		seqan::SeqFileOut seqFileout;	
		if(!seqan::open(seqFileout, seqan::toCString(outFilepaths[label]))){
			cerr<<"ERROR: Could not open "<<outFilepaths[label]<<endl;;
			std::exit(1);
		}
		//output according to labels & doModification flag
		for(unsigned i=0;i<seqan::length(ids);i++){
			if(labels[i]==label){
				if(doModification){
					seqan::Dna5String mod=seqs[i];
					seqan::writeRecord(seqFileout,ids[i],mod);
				}
				else{
					seqan::writeRecord(seqFileout,ids[i],seqs[i]);
				}
			}
		}
	}
}


int main(int argc, char** argv){
	std::string basename=argv[1];
    seqan::CharString seqFilepath=argv[2];
	seqan::StringSet<seqan::CharString> outFilepaths;
	seqan::appendValue(outFilepaths, argv[3]);
	seqan::appendValue(outFilepaths, argv[4]);
    std::string logFilepath=argv[5];
	bool doModification=false;
	
	if(argc==7){
		std::string tmp=argv[6];
		if(tmp=="1"){
			doModification=true;
			cout<<"\tWITH modification"<<endl;
		}
	}

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
	split_fasta(outFilepaths, ids, seqs, labels, doModification);
    for(unsigned label=0;label<2;label++){
		cout<<"\tDONE writing "<<labelCounts[label]<<" seqs to   "<<outFilepaths[label]<<endl;
	}
	
	//output log
	std::ofstream ofs(logFilepath);
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
