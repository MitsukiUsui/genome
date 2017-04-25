#include <iostream>
#include <sstream>
#include <string>
#include <seqan/seq_io.h>

using std::cout;
using std::cerr;
using std::endl;


int main(int argc, char** argv){
//  seqan::String<char> seqFilepath="test.fasta";
//	seqan::String<char> chrFilepath="test_chr.fasta";
//  seqan::String<char> plsFilepath="test_pls.fasta";
    seqan::String<char> seqFilepath=argv[1];
	seqan::String<char> chrFilepath=argv[2];
    seqan::String<char> plsFilepath=argv[3];
    std::string logFilepath=argv[4];
    
    //read input file
    seqan::StringSet<seqan::CharString > ids;
    seqan::StringSet<seqan::Dna5String > seqs;
	seqan::SeqFileIn seqFileIn;
    if(!open(seqFileIn, seqan::toCString(seqFilepath))){
		cerr<<"ERROR: Could not open "<<seqFilepath<<"."<<endl;
		return 1;
	}
	try{
		seqan::readRecords(ids,seqs,seqFileIn);
	}
	catch (seqan::Exception const & e){
		cout<<"ERROR: "<<e.what()<<endl;
		return 1;
	}
    cout<<"DONE reading "<<seqan::length(ids)<<" seqs from"<<seqFilepath<<"."<<endl;


    //open output file
    seqan::SeqFileOut chrFileOut;
    if(!open(chrFileOut, seqan::toCString(chrFilepath))){
		cerr<<"ERROR: Could not open "<<chrFilepath<<"."<<endl;;
		return 1;
	}
    seqan::SeqFileOut plsFileOut;
    if(!open(plsFileOut, seqan::toCString(plsFilepath))){
		cerr<<"ERROR: Could not open "<<plsFilepath<<"."<<endl;
		return 1;
	}

    //output
    int chrCount=0;
    int plsCount=0;
    std::stringstream chrSs;
    std::stringstream plsSs;

    for(unsigned i=0;i<length(ids);i++){
        std::string id_str=seqan::toCString(ids[i]);
        std::string::size_type index=id_str.find("plasmid");

        if(index==std::string::npos){
		    seqan::writeRecord(chrFileOut,ids[i],seqs[i]);
            chrCount++;
            chrSs<<seqan::length(seqs[i])<<":";
        }
        else{
		    seqan::writeRecord(plsFileOut,ids[i],seqs[i]);
            plsCount++;
            plsSs<<seqan::length(seqs[i])<<":";
        }
    }
    cout<<"DONE outputing "<<chrCount<<" seqs to "<<chrFilepath<<"."<<endl;
    cout<<"DONE outputing "<<plsCount<<" seqs to "<<plsFilepath<<"."<<endl;

    //delete last ':'
    std::string chrLength=chrSs.str();
    std::string plsLength=plsSs.str();
    if(chrCount>0){
        chrLength.pop_back();
    }
    if(plsCount>0){
        plsLength.pop_back();
    }

    //write log function
    std::ofstream ofs(logFilepath);
    if (!ofs){
		cerr<<"ERROR: Could not open "<<logFilepath<<"."<<endl;;
		return 1;
    }
    
    ofs<<"test"<<","<<"chromosome"<<","<<chrFilepath<<","<<chrCount<<","<<chrLength<<endl;
    ofs<<"test"<<","<<"plasmid"<<","<<plsFilepath<<","<<plsCount<<","<<plsLength<<endl;

    return 0;
}
