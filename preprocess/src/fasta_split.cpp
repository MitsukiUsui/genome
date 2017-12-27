//
// Created by 薄井光生 on 2017/10/10.
//

#define SEQAN_ENABLE_DEBUG 1

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include "myseqan.h"

#include "fasta_split.h"

using std::cout;
using std::cerr;
using std::endl;

//------------------------------------------------------------
//initialize parser
//------------------------------------------------------------
void set_parser(seqan::ArgumentParser &parser, int argc, char **argv) {/*{{{*/
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
    if (res != seqan::ArgumentParser::PARSE_OK) {
        std::exit(1);
    }
}/*}}}*/


int main(int argc, char **argv) {
    //argument setting
    seqan::ArgumentParser parser("fasta_split");
    set_parser(parser, argc, argv);
    seqan::CharString seqFilepath;
    seqan::getArgumentValue(seqFilepath, parser, 0);
    seqan::StringSet<seqan::CharString> outFilepaths;
    seqan::resize(outFilepaths, 2);
    seqan::getArgumentValue(outFilepaths[0], parser, 1);
    seqan::getArgumentValue(outFilepaths[1], parser, 2);
    seqan::CharString logFilepath;
    seqan::getArgumentValue(logFilepath, parser, 3);
    int convertFlag;
    seqan::getOptionValue(convertFlag, parser, "convert");

    //configuration
    cout << "DONE: configuration" << endl;
    cout << "\tseqFilepath: " << seqFilepath << endl;
    cout << "\tchrFilepath: " << outFilepaths[0] << endl;
    cout << "\tplsFilepath: " << outFilepaths[1] << endl;
    cout << "\tlogFilepath: " << logFilepath << endl;
    cout << "\tConvert mode: " << convertFlag << endl;

    //read input file
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::IupacString> seqs;
    read_fasta(ids, seqs, seqFilepath);
    cout << "DONE: read " << seqan::length(ids) << " seqs" << endl;

    //classify chromosome and plasmid
    std::vector<int> labels;
    classify(labels, ids);

    //split fasta into 2 outFilepaths
    split_fasta(outFilepaths, ids, seqs, labels, convertFlag, logFilepath);
    cout << "DONE: split" << endl;

    return 0;
}