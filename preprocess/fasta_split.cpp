#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>

#include <seqan/seq_io.h>
#include <seqan/modifier.h>
#include <seqan/stream.h>

using std::cout;
using std::cerr;
using std::endl;

template <typename TSeqs>
void read_fasta(seqan::StringSet<seqan::CharString> & ids, TSeqs & seqs, seqan::CharString const & seqFilepath){
	seqan::SeqFileIn seqFileIn;
    if(!seqan::open(seqFileIn, seqan::toCString(seqFilepath))){
		cerr<<"ERROR: Could not open "<<seqFilepath<<endl;
		std::exit(1);
	}
	try{
		seqan::readRecords(ids,seqs,seqFileIn);
	}
	catch (seqan::Exception const & e){
		cout<<"ERROR: "<<e.what()<<endl;
		std::exit(2);
	}
    return;
}

//0...chromosome, 1...plasmid
void judge_type(std::vector<int> & labels, seqan::StringSet<seqan::CharString> const & ids){
    for(unsigned i=0;i<seqan::length(ids);i++){
        std::string id_str=seqan::toCString(ids[i]);
        std::string::size_type index=id_str.find("plasmid");
        if(index==std::string::npos){
            labels.push_back(0);
        }
        else{
            labels.push_back(1);
        }
    }
    return;
}


struct MyFunctor:
    public std::unary_function<char, char>
{
    inline char operator()(char x) const
    {
        char r;
        if(('a'<=x) && (x <= 'z')){
            r = x + ('A'-'a');
        }else{
            r = x;
        }

        switch(r){
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                return r;
            default:
                return 'N';
        }
    }
};


void split_fasta(seqan::StringSet<seqan::CharString> const & outFilepaths, seqan::StringSet<seqan::CharString> const & ids, seqan::StringSet<seqan::CharString> & seqs, std::vector<int> const & labels, bool doModification){
	assert(seqan::length(outFilepaths)==2);
	for(unsigned label=0;label<2;label++){
		seqan::SeqFileOut seqFileout;	
		if(!seqan::open(seqFileout, seqan::toCString(outFilepaths[label]))){
			cerr<<"ERROR: Could not open "<<outFilepaths[label]<<endl;;
			std::exit(1);
		}
		for(unsigned i=0;i<seqan::length(ids);i++){
			if(labels[i]==label){
				if(doModification){
					seqan::ModifiedString< seqan::CharString, seqan::ModView<MyFunctor> > seqMod(seqs[i]);
					seqan::writeRecord(seqFileout,ids[i],seqMod);
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
	std::string tmp=argv[6];
	bool doModification=false;
	if(tmp=="1"){
		doModification=true;
		cout<<"\tWITH modification"<<endl;
	}

    //read input file
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> seqs;
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
		std::stringstream ss;
		for(unsigned i=0;i<seqan::length(ids);i++){
			if(labels[i]==label){
				ss<<seqan::length(seqs[i])<<":";
			}
		}

    	std::string length=ss.str();
		if(labelCounts[label]>0){
			length.pop_back();
		}
    	ofs<<basename<<","<<legends[label]<<","<<outFilepaths[label]<<","<<labelCounts[label]<<","<<length<<endl;
	} 

	cout<<"\tDONE writing log to "<<logFilepath<<endl;
	return 0;
}
