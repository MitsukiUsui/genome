#include <iostream>

#include <seqan/arg_parse.h>

using std::cout;
using std::endl;

//initialize parser
void set_parser(seqan::ArgumentParser & parser, int argc, const char** argv){
    seqan::addUsageLine(parser,
        "\"basename\" \"seqFilepath\" \"chrFilepath\" \"plsFilepath\" \"logFilepath\" [OPTIONS]");
    seqan::addDescription(parser,
        "This program split fasta file into 2 (chromosome, plasmid). "
        "The summary of spliting is written in log file. "
        "Beside, when given modify option, convert sequences into ACGTN sequences simultaneously. "
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
        "m", "modify", "Whether to do modification (to upper & to N)"));

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK){
        std::cerr<<"ERORR: parse failed."<<endl;
        std::exit(1);
    }
}


int main(int argc, char const ** argv){
    seqan::ArgumentParser parser("test_argparse");
    set_parser(parser,argc,argv);

    // Extract option values and print them.
    for(unsigned i=0;i<5;i++){
        seqan::CharString text;
        getArgumentValue(text, parser, i);
        cout<<text<<endl;
    }
    bool toModify= seqan::isSet(parser, "modify");
    if(toModify){
        cout<<"set to modify string"<<endl;
    }else{
        cout<<"not set to modify"<<endl;
    }

    return 0;
}
