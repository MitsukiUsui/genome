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

using std::cout;
using std::cerr;
using std::endl;


//------------------------------------------------------------
//Define labels (0: chromosome, 1:plasmid) according to each ids(heders)
//------------------------------------------------------------
void classify(std::vector<int> & labels,
                seqan::StringSet<seqan::CharString> const & ids){
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


//------------------------------------------------------------
//fill N randomly
//------------------------------------------------------------
int fill_N(seqan::Dna5String & seq,
            std::mt19937 & mt){
    char aa[4]={'T','C','A','G'};
    double aaFreq[4]={0.25,0.25,0.25,0.25};

    int fillCount=0;
    for(unsigned i=0;i<seqan::length(seq);i++){
        if(seq[i]=='N'){
            fillCount++;
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
    return fillCount;
}


//------------------------------------------------------------
//output fasta separately according to labels
//if convertFlag > 0, conversion into Dna5String or DnaString is carried out simultaneously
//------------------------------------------------------------
void split_fasta(seqan::StringSet<seqan::CharString> const & outFilepaths,
                 seqan::StringSet<seqan::CharString> const & ids,
                 seqan::StringSet<seqan::IupacString> & seqs, //!!! seqs are changed by fill_N? !!!!
                 std::vector<int> const & labels,
                 int convertFlag,
                 seqan::CharString const & logFilepath){

    assert(seqan::length(outFilepaths)==2);
    std::random_device rd;
    std::mt19937 mt(rd());

    //output log
    std::ofstream ofs(seqan::toCString(logFilepath));
    if (!ofs){
        cerr<<"ERROR: Could not open "<<logFilepath<<endl;;
        exit(1);
    }

    for(unsigned label=0;label<2;label++){
        //open file
        seqan::SeqFileOut seqFileout;
        if(!seqan::open(seqFileout, seqan::toCString(outFilepaths[label]))){
            cerr<<"ERROR: Could not open "<<outFilepaths[label]<<endl;
            std::exit(1);
        }
        //output according to labels & convertFlag
        for(unsigned i=0;i<seqan::length(ids);i++){
            if(labels[i]==label){
                int fillCount = 0;
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
                        fillCount = fill_N(dna5, mt);
                        seqan::DnaString dna=dna5;
                        seqan::writeRecord(seqFileout,ids[i],dna);
                        break;
                    }
                }

                std::string legend = (label==0 ? "chromosome" : "plasmid");
                ofs<<outFilepaths[label]<<",\""<<ids[i]<<"\","<<legend<<","<<seqan::length(seqs[i])<<endl;
                cout<<"DONE: output \""<<ids[i]<<"\" to "<<legend<<endl;
            }
        }
    }
}

//------------------------------------------------------------
//initialize parser
//------------------------------------------------------------
void set_parser(seqan::ArgumentParser & parser, int argc, char ** argv){/*{{{*/
    seqan::addUsageLine(parser,
        "\"seqFilepath\" \"chrFilepath\" \"plsFilepath\" \"logFilepath\" [OPTIONS]");
    seqan::addDescription(parser,
        "This program split input fasta file into 2 output fasta (chromosome, plasmid). "
        "The summary of spliting is written in log file. "
        "Beside, when given convert option, convert sequences into Dna5 or Dna sequence while splitting. "
        );
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
    seqan::CharString seqFilepath;
    seqan::getArgumentValue(seqFilepath, parser, 0);
    seqan::StringSet<seqan::CharString> outFilepaths;
    seqan::resize(outFilepaths, 2);
    seqan::getArgumentValue(outFilepaths[0], parser, 1);
    seqan::getArgumentValue(outFilepaths[1], parser, 2);
    seqan::CharString logFilepath;
    seqan::getArgumentValue(logFilepath, parser, 3);
    int convertFlag;
    seqan::getOptionValue(convertFlag, parser,"convert");


    //read input file
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::IupacString> seqs;
    read_fasta(ids,seqs,seqFilepath);
    cout<<"DONE: read "<<seqan::length(ids)<<" seqs from "<<seqFilepath<<endl;

    //classify chromosome and plasmid
    std::vector<int> labels;
    classify(labels,ids);

    //split fasta into 2 outFilepaths
    split_fasta(outFilepaths, ids, seqs, labels, convertFlag, logFilepath);
    return 0;
}
