#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <random>

#define SEQAN_ENABLE_DEBUG 1

#include <seqan/arg_parse.h>
#include "myseqan.h"
#include "myutil.h"

#include "geneticcode.h"
#include "shuffle_genome.h"

using std::cout;
using std::cerr;
using std::endl;


//------------------------------------------------------------
//parse arguments 
//------------------------------------------------------------
void set_parser(seqan::ArgumentParser &parser, int argc, char **argv) {/*{{{*/
    seqan::addUsageLine(parser,
                        "\"seqFilepath\" \"gffFilepath\"");
    seqan::addDescription(parser,
                          "This program shuffle genomes preseving CDS regions."
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
    if (res != seqan::ArgumentParser::PARSE_OK) {
        std::exit(1);
    }
}/*}}}*/


int main(int argc, char **argv) {
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
    int shuffleMode[2] = {};
    seqan::getArgumentValue(shuffleMode[0], parser, 3);
    seqan::getArgumentValue(shuffleMode[1], parser, 4);
    bool dryrun = seqan::isSet(parser, "dryrun");

    //configuration
    cout << "DONE: read configuration:" << endl;
    cout << "\tseqFilepath: " << seqFilepath << endl;
    cout << "\tgffFilepath: " << gffFilepath << endl;
    cout << "\toutFilepath: " << outFilepath << endl;
    cout << "\tshuffle mode: " << shuffleMode[0] << ", " << shuffleMode[1] << endl;
    if (dryrun) {
        cout << "\texecution: DRYRUN" << endl;
    }

    //------------------------------------------------------------
    //read fasta
    //------------------------------------------------------------
    seqan::String<seqan::CharString> seqIds;
    seqan::StringSet<seqan::DnaString> seqs;
    read_fasta(seqIds, seqs, seqFilepath);

    for (seqan::CharString  &seqId : seqIds) {
        seqId = seqan::CharString(split(seqan::toCString(seqId), ' ')[0]);//trim ids
    }
    cout << "DONE: read " << seqan::length(seqIds) << " seqs from " << seqFilepath << endl;

    //------------------------------------------------------------
    //read gff
    //------------------------------------------------------------
    seqan::String<seqan::GffRecord> gffs;
    read_gff(gffs, gffFilepath);
    cout << "DONE: read " << seqan::length(gffs) << " gffs from " << gffFilepath << endl;

    //------------------------------------------------------------
    //update GeneticCode with gffs
    //------------------------------------------------------------
    GeneticCode geneticCode(11);
    std::vector< std::vector<MyCDS> > myCDS_vecvec;

    get_myCDS_vecvec(myCDS_vecvec, gffs, seqs, seqIds, geneticCode);
    for (int seqIdx = 0; seqIdx < seqan::length(seqs); seqIdx++) {
        cout << "\tDONE: found " << myCDS_vecvec[seqIdx].size()<<" typical CDSs for "<<seqIds[seqIdx]<<endl;
    }

    update_genetic_code(geneticCode, myCDS_vecvec, seqs);
    geneticCode.__show_freq(true);

    std::vector<ShuffleRegion> shuffleRegions;
    get_shuffle_region(shuffleRegions, shuffleMode, seqs, myCDS_vecvec);
    cout << "DONE: find " << seqan::length(shuffleRegions) << " shuffle regions" << endl;

    shuffle_genome(seqs, shuffleRegions, geneticCode);

    geneticCode.clear_count();
    update_genetic_code(geneticCode, myCDS_vecvec, seqs);
    geneticCode.__show_freq(true);

    cout << "DONE: output to " << outFilepath << endl;
    return 0;
}
