#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <random>


#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include "myutil.h"
#include "myseqan.h"

//------------------------------------------------------------
//Define labels (0: chromosome, 1:plasmid) according to each ids(heders)
//------------------------------------------------------------
void classify(std::vector<int> &labels,
              seqan::StringSet<seqan::CharString> const &ids) {
    labels.resize(seqan::length(ids));
    for (unsigned i = 0; i < seqan::length(ids); i++) {
        std::string id_str = seqan::toCString(ids[i]);
        std::string::size_type index = id_str.find("plasmid");
        if (index == std::string::npos) {
            labels[i] = 0;
        } else {
            labels[i] = 1;
        }
    }
    return;
}


//------------------------------------------------------------
//fill N randomly
//------------------------------------------------------------
int fill_N(seqan::Dna5String &seq,
           std::mt19937 &mt) {
    char base[4] = {'T', 'C', 'A', 'G'};
    double baseFreq[4] = {0.25, 0.25, 0.25, 0.25};

    int fillCount = 0;
    for (unsigned i = 0; i < seqan::length(seq); i++) {
        if (seq[i] == 'N') {
            fillCount++;
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            double d = dist(mt);
            double cumSum = 0;
            for (int idx = 0; idx < 4; idx++) {
                cumSum += baseFreq[idx];
                if (cumSum >= d) {
                    seq[i] = base[idx];
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
void split_fasta(seqan::StringSet<seqan::CharString> const &outFilepaths,
                 seqan::StringSet<seqan::CharString> const &ids,
                 seqan::StringSet<seqan::IupacString> &seqs, //!!! seqs are changed by fill_N? !!!!
                 std::vector<int> const &labels,
                 int convertFlag,
                 seqan::CharString const &logFilepath) {

    assert(seqan::length(outFilepaths) == 2);
    std::random_device rd;
    std::mt19937 mt(rd());

    Logger logger(seqan::toCString(logFilepath));

    for (unsigned label = 0; label < 2; label++) {
        //open file
        seqan::SeqFileOut seqFileout;
        if (!seqan::open(seqFileout, seqan::toCString(outFilepaths[label]))) {
            std::cerr << "ERROR: Could not open " << outFilepaths[label] << std::endl;
            std::exit(1);
        }
        //output according to labels & convertFlag
        for (unsigned i = 0; i < seqan::length(ids); i++) {
            if (labels[i] == label) {
                int fillCount = 0;
                switch (convertFlag) {
                    case 0: {
                        seqan::writeRecord(seqFileout, ids[i], seqs[i]);
                        break;
                    }
                    case 1: {
                        seqan::Dna5String dna5 = seqs[i];
                        seqan::writeRecord(seqFileout, ids[i], dna5);
                        break;
                    }
                    case 2: {
                        seqan::Dna5String dna5 = seqs[i];
                        fillCount = fill_N(dna5, mt);
                        seqan::DnaString dna = dna5;
                        seqan::writeRecord(seqFileout, ids[i], dna);
                        break;
                    }
                }

                std::string legend = (label == 0 ? "chromosome" : "plasmid");

                logger << std::string(seqan::toCString(outFilepaths[label]))
                       << std::string(seqan::toCString(ids[i]))
                       << legend
                       << seqan::length(seqs[i])
                       << fillCount;
                logger.flush();
            }
        }
    }
}
