#ifndef GENETICCODE_H_
#define GENETICCODE_H_

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <random>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

class GeneticCode {
    std::string id2aa = "*ACDEFGHIKLMNPQRSTVWY";//1 stop, 20 aa
    seqan::StringSet<seqan::DnaString> id2codon;//64 codon

    int codon2aa[64]; // codonId -> aaId
    std::vector< std::vector <int> > aa2codon = {{}}; //aaId -> codonId
    int codonCount[64] = {};
    double codonFrequency[64] = {};

    void calc_freq();
    void __GeneticCode(std::string str);

public:
    std::string geneticCode_str;

    //------------------------------------------------------------
    // visible for testing
    //------------------------------------------------------------
    template<typename TSeq> int codon_encode(TSeq const &codon) const;//return 0~63, return -1 if N included
    seqan::DnaString codon_decode(int codonId) const;//*** need to be fixed to just returning array value, without calculation
    int aa_encode(char aa) const; //return 0~21
    char aa_decode(int aaId) const;

    GeneticCode(std::string geneticCode_str);
    GeneticCode(int transl_table);
    ~GeneticCode(){};

    int get_count(int codonId) {
        assert (codonId>=0 && codonId<64);
        return codonCount[codonId];
    }

    double get_frequency(int codonId) {
        assert (codonId>=0 && codonId<64);
        return codonFrequency[codonId];
    }

    void update_count(seqan::Dna5String &seq);
    void clear_count();

    template<typename TSeq> void synonymous_sub(TSeq &codon, std::mt19937 &mt) const;
    template<typename TSeq> seqan::String<seqan::AminoAcid> translate(TSeq seq) const;//ISSUE does not support ambiguous base for now
    template<typename TSeq> bool is_stop_codon(TSeq codon) const;

    void __show() const;
    void __show_freq(bool isFreq = true) const;//true...freq, false...count;
};


//============================================================
//Private Methods
//============================================================

void GeneticCode::__GeneticCode(std::string geneticCode_str) {
    if (geneticCode_str.length() != 64) {
        std::cerr << "The length of genetic code is not 64, but " << geneticCode_str.length() << std::endl;
        exit(1);
    }
    this->geneticCode_str = geneticCode_str;

    //initialize id2codon
    seqan::resize(id2codon, 64);
    std::string bases = "TCAG";
    for (int i = 0; i<4; i++) {
        for (int j =0 ; j<4; j++) {
            for (int k =0 ; k<4; k++) {
                seqan::DnaString codon = "NNN";
                codon[0] = bases[i];
                codon[1] = bases[j];
                codon[2] = bases[k];
                id2codon[codon_encode(codon)]=codon;
            }
        }
    }

    //initialize codon2aa
    for (int codonId = 0 ; codonId < 64; codonId++) {
        codon2aa[codonId] = aa_encode(geneticCode_str[codonId]);
    }

    //initialize aa2codon
    for (int aaId =0 ; aaId < 21; aaId++ ) {
        std::vector<int> v;
        aa2codon = std::vector<int>;
    }

    //codonToAa, count in aaBucket
    int aaBucket[21] = {};
    for (int codonId = 0; codonId < 64; codonId++) {
        int aaId = aa_encode(geneticCode_str[codonId]);
        codonToAa[codonId] = aaId;
        aaBucket[aaId]++;
    }

    //convert bucket to index
    aaIndex[0] = 0;
    for (int i = 1; i < 22; i++) {
        aaIndex[i] = aaIndex[i - 1] + aaBucket[i - 1];
    }
    assert(aaIndex[21] == 64);
    //copy index to bucket
    for (int i = 0; i < 21; i++) {
        aaBucket[i] = aaIndex[i];
    }
    //aaToCodon with incrementing bucket
    for (int codonId = 0; codonId < 64; codonId++) {
        int aaId = aa_encode(geneticCode_str[codonId]);
        aaToCodon[aaBucket[aaId]++] = codonId;
    }
}

template<typename TSeq>
inline int GeneticCode::codon_encode(TSeq const &codon) const {
    assert(seqan::length(codon) == 3);

    int codonId = 0;
    int base = 1;
    for (int i = 0; i < 3; i++) {
        char c = codon[2 - i];
        switch (c) {
            case 'T':
                codonId += 0 * base;
                break;
            case 'C':
                codonId += 1 * base;
                break;
            case 'A':
                codonId += 2 * base;
                break;
            case 'G':
                codonId += 3 * base;
                break;
            default:
                std::cerr << "ERROR: Unusual nucleotide " << codon[2 - i] << " is included." << std::endl;
                return -1;
        }
        base *= 4;
    }

    return codonId;
}

inline seqan::DnaString GeneticCode::codon_decode(int codonId) const {
    assert(codonId >= 0 & codonId < 64);
    return id2codon[codonId];
}

inline int GeneticCode::aa_encode(char aa) const {
    int aaId = -1;
    for (int i = 0; i < geneticCode_str.length(); i++) {
        if (aa == geneticCode_str[i]) {
            aaId = i;
            break;
        }
    }
    if (aaId == -1) {
        std::cerr << "ERROR: Undefined amino acid " << aa << std::endl;
    }
    return aaId;
}

inline char GeneticCode::aa_decode(int aaId) const {/*{{{*/
    assert(aaId >= 0 & aaId < 21);
    return id2aa[aaId];
}/*}}}*/



//============================================================
//(de)Constructor
//============================================================

GeneticCode::GeneticCode(std::string str) {
    __GeneticCode(str);
}

GeneticCode::GeneticCode(int transl_table) {
    switch (transl_table) {
        case 4:
            __GeneticCode("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
            break;
        case 11:
            __GeneticCode("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
            break;
        default:
            std::cerr << "ERROR: tranl_table = " << transl_table << " is not supported" << std::endl;
            exit(1);
    }
}

//============================================================
//Public Methods
//============================================================

void GeneticCode::clear_count() {/*{{{*/
    for (int i = 0; i < 64; i++) {
        codonCount[i] = 0;
        codonFreq[i] = 0;
    }
}/*}}}*/

void GeneticCode::update_count(seqan::Dna5String &seq) {/*{{{*/
    assert(seqan::length(seq) % 3 == 0);
    int aaLength = seqan::length(seq) / 3;
    for (int i = 0; i < aaLength; i++) {
        seqan::Infix<seqan::Dna5String>::Type codon = seqan::infix(seq, 3 * i, 3 * (i + 1));
        int codonId = codon_encode(codon);
        if (codonId != -1)
            codonCount[codonId]++;
    }
}/*}}}*/

void GeneticCode::calc_freq() {/*{{{*/
    for (int aaId = 0; aaId < 21; aaId++) {
        //calculate count sum
        int countSum = 0;
        for (int i = aaIndex[aaId]; i < aaIndex[aaId + 1]; i++) {
            countSum += codonCount[aaToCodon[i]];
        }
        //calculate frequency
        for (int i = aaIndex[aaId]; i < aaIndex[aaId + 1]; i++) {
            int codonId = aaToCodon[i];
            if (countSum == 0) {//to avoid division by 0
                codonFreq[codonId] = 0;
            } else {
                codonFreq[codonId] = (double) codonCount[codonId] / countSum;
            }
        }
    }
}/*}}}*/

template<typename TSeq>
void GeneticCode::synonymous_sub(TSeq &codon, std::mt19937 &mt) const {/*{{{*/
    int codonIdOld = codon_encode(codon);
    int aaId = codonToAa[codonIdOld];

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double d = dist(mt);
    double cumSum = 0;
    //calc codonIdNew
    int codonIdNew = -1;
    for (int i = aaIndex[aaId]; i < aaIndex[aaId + 1]; i++) {
        int codonId = aaToCodon[i];
        cumSum += codonFreq[codonId];
        if (cumSum >= d) {
            codonIdNew = codonId;
            break;
        }
    }

    if (codonIdNew == -1){
        std::cerr<<"ERROR: no frequency is defined."<<std::endl;
        codonIdNew = aaToCodon[aaIndex[aaId]];
    }

    seqan::DnaString codonNew = codon_decode(codonIdNew);
    for (int i = 0; i < 3; i++) {
        codon[i] = codonNew[i];
    }
}/*}}}*/

template<typename TSeq>
seqan::String<seqan::AminoAcid> GeneticCode::translate(TSeq seq) const {/*{{{*/
    assert(seqan::length(seq) % 3 == 0);
    int aaLength = seqan::length(seq) / 3;
    seqan::String<seqan::AminoAcid> seq_aa;
    seqan::resize(seq_aa, aaLength);
    for (int i = 0; i < aaLength; i++) {
        seqan::Infix<seqan::Dna5String>::Type codon = seqan::infix(seq, 3 * i, 3 * (i + 1));
        seq_aa[i] = aa_decode(codonToAa[codon_encode(codon)]);
    }
    return seq_aa;
}/*}}}*/

template<typename TSeq>
bool GeneticCode::is_stop_codon(TSeq codon) const {/*{{{*/
    int aaId = codonToAa[codon_encode(codon)];
    if (aaId == 0) {
        return true;
    } else {
        return false;
    }
}/*}}}*/

void GeneticCode::__show() const {/*{{{*/
    std::cout << "genetic code: " << geneticCode_str << std::endl;
    std::cout << "aa: " << aa_str << std::endl;

    //output codon:aa pair
    for (int codonId = 0; codonId < 64; codonId++) {
        std::cout << codon_decode(codonId) << "-" << aa_decode(codonToAa[codonId]) << ", ";
    }
    std::cout << std::endl;

    //output aa:codonLst pair
    for (int aaId = 0; aaId < 21; aaId++) {
        std::cout << aa_decode(aaId) << "-[";
        for (int i = aaIndex[aaId]; i < aaIndex[aaId + 1]; i++) {
            std::cout << codon_decode(aaToCodon[i]) << ",";
        }
        std::cout << "], ";
    }
    std::cout << std::endl;
}/*}}}*/

void GeneticCode::__show_freq(bool isFreq) const {/*{{{*/
    for (int aaId = 0; aaId < 21; aaId++) {
        std::cout << aa_decode(aaId) << "\t";

        for (int i = aaIndex[aaId]; i < aaIndex[aaId + 1]; i++) {
            int codonId = aaToCodon[i];
            if (isFreq) {
                std::cout << codon_decode(codonId) << ": " << codonFreq[codonId] << ", ";
            } else {
                std::cout << codon_decode(codonId) << ": " << codonCount[codonId] << ", ";
            }
        }
        std::cout << std::endl;
    }
}/*}}}*/


#endif //GENETICCODE_H_
