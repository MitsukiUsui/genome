//
// Created by 薄井光生 on 2017/10/10.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "myseqan.h"
#include "fasta_split.h"

using std::cout;
using std::cerr;
using std::endl;

TEST(myseqan, test_read_fasta) {
    seqan::CharString seqFilepath = "../unittest/test.fasta";
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::IupacString> seqs;
    read_fasta(ids, seqs, seqFilepath);

    EXPECT_EQ(seqan::length(ids), 2);
    EXPECT_EQ(seqan::length(seqs), 2);

    EXPECT_TRUE(ids[0] == "seq1");
    EXPECT_TRUE(seqs[0] == "AaTt");
    EXPECT_TRUE(ids[1] == "seq2 plasmid");
    EXPECT_TRUE(seqs[1] == "GNCn");
}

TEST(fasta_split, test_classify) {
    seqan::StringSet<seqan::CharString> ids;
    seqan::CharString id1 = "id1 complete";
    seqan::CharString id2 = "id2 plasmid complete";
    seqan::CharString id3 = "id3 PLASMID complete";

    seqan::appendValue(ids, id1);
    seqan::appendValue(ids, id2);
    seqan::appendValue(ids, id3);

    std::vector<int> labels;
    classify(labels, ids);

    EXPECT_EQ(seqan::length(labels), 3);
    EXPECT_EQ(labels[0], 0);
    EXPECT_EQ(labels[1], 1);
    EXPECT_EQ(labels[2], 0); //currently classify is case sensitive
}

TEST(fasta_split, test_fill_N) {
    seqan::Dna5String seq1 = "ATGC";
    seqan::Dna5String seq2 = "NNNN";

    std::random_device rd;
    std::mt19937 mt(rd());

    fill_N(seq1, mt);
    EXPECT_TRUE(seq1 == "ATGC");

    fill_N(seq2, mt);
    //WRITE ME
}

TEST(fasta_split, test_auto_fill) {
    seqan::CharString seq1 = "NYgG";
    seqan::Dna5String seq2 = seq1;
    cout<<seq2<<endl;
    EXPECT_TRUE(seq2 == "NNgG");
}


//TEST(fasta_split, test_fasta_split) {
//    seqan::CharString seqFilepath = "/Users/mitsuki/sandbox/genome/preprocess/unittest/test.fasta";
//    seqan::StringSet<seqan::CharString> ids;
//    seqan::StringSet<seqan::IupacString> seqs;
//
//    read_fasta(ids, seqs, seqFilepath);
//
//    std::vector<int> labels;
//    classify(labels, ids);
//
//    seqan::StringSet<seqan::CharString> outFilepaths;
//    seqan::CharString chrFilepath = "/Users/mitsuki/sandbox/genome/preprocess/unittest/chr.fasta";
//    seqan::CharString plsFilepath = "/Users/mitsuki/sandbox/genome/preprocess/unittest/pls.fasta";
//    seqan::appendValue(outFilepaths, chrFilepath);
//    seqan::appendValue(outFilepaths, plsFilepath);
//
//    int convertFlag=2;
//    seqan::CharString logFilepath = "/Users/mitsuki/sandbox/genome/preprocess/unittest/test.log";
//    split_fasta(outFilepaths, ids, seqs, labels, convertFlag, logFilepath);
//
//    EXPECT_EQ(1, 0); // check the result manually
//}

