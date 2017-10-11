//
// Created by 薄井光生 on 2017/10/10.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <map>
#include "geneticcode.h"
#include "myutil.h"

#define SEQAN_ENABLE_DEBUG 1

using std::cout;
using std::cerr;
using std::endl;

namespace {
    class GeneticCodeTest: public testing::Test {
    public:
        GeneticCode gc = GeneticCode(11);
    };
}

TEST_F(GeneticCodeTest, test_codon_encode) {
    seqan::DnaString seq_ds= "GGG";
    seqan::Infix<seqan::DnaString>::Type seq_if=seqan::infix(seq_ds, 0, 3);
    std::string seq_str = "GGG";

    ASSERT_EQ(gc.codon_encode(seq_ds), 63);
    ASSERT_EQ(gc.codon_encode(seq_if), 63);
    ASSERT_EQ(gc.codon_encode(seq_str), 63);
}

TEST_F(GeneticCodeTest, test_codon_encode_fail) {
    std::string codon = "ggg";
    ASSERT_EQ(gc.codon_encode(codon), -1);

    codon = "NNN";
    ASSERT_EQ(gc.codon_encode(codon), -1);
}

TEST_F(GeneticCodeTest, test_codon_count) {
    seqan::Dna5String seq1 = "TTTGGG";
    seqan::Dna5String seq2 = "TTTNNNGGG";

    int codonId1 = gc.codon_encode("TTT");
    int codonId2 = gc.codon_encode("GGG");
    int codonId3 = gc.codon_encode("CCC");

    ASSERT_EQ(gc.get_count(codonId1), 0);
    ASSERT_EQ(gc.get_count(codonId2), 0);
    ASSERT_EQ(gc.get_count(codonId3), 0);

    gc.update_count(seq1);

    ASSERT_EQ(gc.get_count(codonId1), 1);
    ASSERT_EQ(gc.get_count(codonId2), 1);
    ASSERT_EQ(gc.get_count(codonId3), 0);

    gc.update_count(seq2);

    ASSERT_EQ(gc.get_count(codonId1), 2);
    ASSERT_EQ(gc.get_count(codonId2), 2);
    ASSERT_EQ(gc.get_count(codonId3), 0);

    gc.clear_count();
    ASSERT_EQ(gc.get_count(codonId1), 0);
    ASSERT_EQ(gc.get_count(codonId2), 0);
    ASSERT_EQ(gc.get_count(codonId3), 0);
}

TEST_F(GeneticCodeTest, test_codon_freq) {
    int codonId1 = gc.codon_encode("TCT");
    int codonId2 = gc.codon_encode("TCC");
    int codonId3 = gc.codon_encode("TCA");
    int codonId4 = gc.codon_encode("TCG");

    seqan::Dna5String seq1 = "TCTTCC";
    seqan::Dna5String seq2 = "TCTTCA";

    ASSERT_EQ(gc.get_freq(codonId1), 0);
    ASSERT_EQ(gc.get_freq(codonId2), 0);
    ASSERT_EQ(gc.get_freq(codonId3), 0);
    ASSERT_EQ(gc.get_freq(codonId4), 0);

    gc.update_count(seq1);

    ASSERT_EQ(gc.get_freq(codonId1), 0);
    ASSERT_EQ(gc.get_freq(codonId2), 0);
    ASSERT_EQ(gc.get_freq(codonId3), 0);
    ASSERT_EQ(gc.get_freq(codonId4), 0);

    gc.calc_freq();

    ASSERT_EQ(gc.get_freq(codonId1), 0.5);
    ASSERT_EQ(gc.get_freq(codonId2), 0.5);
    ASSERT_EQ(gc.get_freq(codonId3), 0);
    ASSERT_EQ(gc.get_freq(codonId4), 0);

    gc.update_count(seq2);
    gc.calc_freq();

    ASSERT_EQ(gc.get_freq(codonId1), 0.5);
    ASSERT_EQ(gc.get_freq(codonId2), 0.25);
    ASSERT_EQ(gc.get_freq(codonId3), 0.25);
    ASSERT_EQ(gc.get_freq(codonId4), 0);
}

TEST_F(GeneticCodeTest, test_synonymous_sub){
    seqan::Dna5String seq = "TCTTCTTCCTCA";
    gc.update_count(seq);
    gc.calc_freq();

    std::vector<int> ids;
    ids.push_back(gc.codon_encode("TCT"));
    ids.push_back(gc.codon_encode("TCC"));
    ids.push_back(gc.codon_encode("TCA"));
    ids.push_back(gc.codon_encode("TCG"));

    std::map<int, int> codonCount;
    for (int codonId : ids) {
        codonCount[codonId]=0;
    }

    std::random_device rd;
    std::mt19937 mt(rd());
    std::string codon = "TCT";
    for (int i =0;i < 1000; i++){
        gc.synonymous_sub(codon, mt);
        int codonId = gc.codon_encode(codon);
        ASSERT_TRUE(is_in(codonId, ids));
        codonCount[codonId]+=1;
    }

    ASSERT_TRUE(codonCount[ids[0]]>=450 && codonCount[ids[0]]<=550);//Expected 500
    ASSERT_TRUE(codonCount[ids[1]]>=200 && codonCount[ids[1]]<=300);//Expected 250
    ASSERT_TRUE(codonCount[ids[2]]>=200 && codonCount[ids[2]]<=300);//Expected 250
    ASSERT_EQ(codonCount[ids[3]], 0);
}
