//
// Created by 薄井光生 on 2017/10/11.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include "myseqan.h"

#define SEQAN_ENABLE_DEBUG 1

using std::cout;
using std::cerr;
using std::endl;

TEST(myseqan, test_read_gff) {
    seqan::String<seqan::GffRecord> records;
    seqan::CharString gffFilepath = "/Users/mitsuki/sandbox/genome/genomeshuffle/unittest/test.gff";
    read_gff(records, gffFilepath);

    ASSERT_EQ(seqan::length(records), 3);

    ASSERT_TRUE(records[0].ref == "seq1");
    ASSERT_TRUE(records[0].type == "gene");
    ASSERT_EQ(records[0].beginPos, 99); //seqan converts "1-based closed intervals" to "0-based half-open intervals"
    ASSERT_EQ(records[0].endPos, 200);
    ASSERT_EQ(records[0].strand, '+');
}

TEST(myseqan, test_read_gff_error) {
    seqan::String<seqan::GffRecord> records;
    seqan::CharString gffFilepath = "/Users/mitsuki/rsyncdir/files/964_genomic.gff";
    read_gff(records, gffFilepath);
    cout<<seqan::length(records)<<endl;
}

TEST(seqan, test_infix) {
    seqan::Dna5String seq_ds= "NNN";
    seqan::DnaString seq_if = seqan::infix(seq_ds, 0, 3);
    ASSERT_TRUE(seq_if=="AAA"); //N will be convert to NNN
    ASSERT_TRUE(seq_ds=="NNN"); //won't change the base string

    std::string seq_str = "ATGC";
    seq_if = seqan::infix(seq_str, 0, 4);
    ASSERT_TRUE(seq_if=="ATGC"); //infix can be used for std::string
}

TEST(seqan, test_length) {
    seqan::DnaString seq_ds = "NNN";
    ASSERT_EQ(seqan::length(seq_ds), 3);

    std::string seq_str = "NNN";
    ASSERT_EQ(seqan::length(seq_str), 3);
}
