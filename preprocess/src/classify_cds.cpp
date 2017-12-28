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

using std::cout;
using std::cerr;
using std::endl;

int classify_cds(seqan::GffRecord const &gff,
                 seqan::DnaString const &seq,
                 GeneticCode const &geneticCode) {
    if (gff.endPos > seqan::length(seq)) {
        return 1;
    }

    int length = gff.endPos - gff.beginPos;
    if (length < 6 | (length % 3) != 0) {
        return 2;
    }

    seqan::DnaString subSeq = seqan::infix(seq, gff.beginPos, gff.endPos); //create new object in order to use Infix.

    if (gff.strand == '-') {
        seqan::reverseComplement(subSeq);
    }

    int aaLength = length / 3;
    for (int i = 0; i < aaLength - 1; i++) {
        seqan::Infix<seqan::DnaString>::Type codon = seqan::infix(subSeq, 3 * i, 3 * (i + 1));
        if (geneticCode.is_stop_codon(codon))
            return 4;
    }

    seqan::Infix<seqan::DnaString>::Type codon = seqan::infix(subSeq, 3 * (aaLength - 1), 3 * aaLength);
    if (!geneticCode.is_stop_codon(codon)) {
        return 5;
    }
    return 0;
}

void update_gff(seqan::String<seqan::GffRecord> &gffs,
               seqan::StringSet<seqan::DnaString> const &seqs,
               seqan::String<seqan::CharString> const &seqIds,
               GeneticCode const geneticCode) {

    int numSeqs = seqan::length(seqs);
    int numGffs = seqan::length(gffs);

    int cdsCount=0;
    int typicalCount=0;

    for (int gffIdx = 0; gffIdx < numGffs; gffIdx++) {
        seqan::GffRecord &gff = gffs[gffIdx];
        if (gff.type == "CDS") {
            cdsCount += 1;

            int seqIdx = -1;
            for (int i = 0; i < numSeqs; i++) {
                if (gff.ref == seqIds[i]) {
                    seqIdx = i;
                    break;
                }
            }

            if (seqIdx != -1) {
                //update gff
                int cdsLabel = classify_cds(gff, seqs[seqIdx], geneticCode);
                seqan::appendValue(gff.tagNames, "cds_label");
                seqan::appendValue(gff.tagValues, std::to_string(cdsLabel));

                if (cdsLabel == 0) {
                    typicalCount += 1;
                }
            }
        }
    }

    cout << "DONE: classified " << cdsCount << " CDSs" << endl;
    double typicalPer=(double)typicalCount/cdsCount*100;
    cout << "\tDONE: found " <<typicalCount <<" ("<<typicalPer <<"%) typical CDSs"<<endl;
}

//------------------------------------------------------------
//parse arguments 
//------------------------------------------------------------
void set_parser(seqan::ArgumentParser &parser, int argc, char **argv) {/*{{{*/
    seqan::addUsageLine(parser,
                        "\"seqFilepath\" \"gffFilepath\"");
    seqan::addDescription(parser,
                          "This program classify cds and write the results to gff description field."
    );
    addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::INPUT_FILE, "gffFilepath"));
    addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::INPUT_FILE, "seqFilepath"));
    addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::OUTPUT_FILE, "outFilepath"));
    addOption(parser, seqan::ArgParseOption(
            "c", "code", "Genetic code (11: bacterial, 4:Mycoplasma)",
            seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "code", "11");

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK) {
        std::exit(1);
    }
}/*}}}*/


int main(int argc, char **argv) {
    //------------------------------------------------------------
    //parse arguments
    //------------------------------------------------------------
    seqan::ArgumentParser parser("classify_cds");
    set_parser(parser, argc, argv);

    seqan::CharString gffFilepath;
    seqan::getArgumentValue(gffFilepath, parser, 0);
    seqan::CharString seqFilepath;
    seqan::getArgumentValue(seqFilepath, parser, 1);
    seqan::CharString outFilepath;
    seqan::getArgumentValue(outFilepath, parser, 2);
    int geneticCodeId;
    seqan::getOptionValue(geneticCodeId, parser, "code");

    //configuration
    cout << "DONE: read configuration:" << endl;
    cout << "\tseqFilepath: " << seqFilepath << endl;
    cout << "\tgffFilepath: " << gffFilepath << endl;
    cout << "\toutFilepath: " << outFilepath << endl;
    cout << "\tgenetic code: " << geneticCodeId << endl;

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
    //edit gff
    //------------------------------------------------------------
    GeneticCode geneticCode(geneticCodeId);
    update_gff(gffs, seqs, seqIds, geneticCode);
    write_gff(gffs, outFilepath);
    cout << "DONE: output updated gff to " << outFilepath << endl;

    return 0;
}
