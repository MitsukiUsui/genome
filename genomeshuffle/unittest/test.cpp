//
// Created by 薄井光生 on 2017/10/11.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include "shuffle_genome.h"
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

seqan::GffRecord create_gffRecord(seqan::CharString ref,
                                  int beginPos,
                                  int endPos,
                                  char strand) {
    seqan::GffRecord gff;
    gff.type = "CDS";
    gff.ref=ref;
    gff.beginPos=beginPos;
    gff.endPos=endPos;
    gff.strand=strand;
    return gff;
}

TEST(genome_shuffle, test_classify_gff) {
    seqan::String<seqan::CharString> seqIds;
    seqan::appendValue(seqIds, "seq1");
    seqan::appendValue(seqIds, "seq2");

    seqan::StringSet<seqan::DnaString> seqs;
    seqan::appendValue(seqs, "ATGC");
    seqan::appendValue(seqs, "AAAA");

    seqan::String<seqan::GffRecord> gffs;
    seqan::appendValue(gffs, create_gffRecord("seq1", 0, 99, '+'));

    //WRITE ME
    ASSERT_EQ(0, 1);
}

TEST(genome_shuffle, test_get_shuffle_region) {
    seqan::String<seqan::GffRecord> gffs;
    seqan::appendValue(gffs, create_gffRecord("seq1", 3, 12, '+'));
//    seqan::appendValue(gffs, create_gffRecord("seq2", 0, 12, '+'));
//    seqan::appendValue(gffs, create_gffRecord("seq2", 7, 19, '+'));

    std::vector<int> gffLabels(seqan::length(gffs), 0);

    seqan::String<seqan::CharString> seqIds;
    seqan::StringSet<seqan::Dna5String> seqs;
    read_fasta(seqIds, seqs, "/Users/mitsuki/sandbox/genome/genomeshuffle/unittest/test.fasta");

    int shuffleMode[2] = {1, 3};
    std::vector<ShuffleRegion> shuffleRegions;
    get_shuffle_region(shuffleRegions, shuffleMode, seqIds, seqs, gffLabels, gffs);

    ASSERT_EQ(shuffleRegions.size(), 3);

    ASSERT_TRUE(shuffleRegions[0].seqId == "seq1");
    ASSERT_EQ(shuffleRegions[0].start, 0);
    ASSERT_EQ(shuffleRegions[0].end, 3);
    ASSERT_EQ(shuffleRegions[0].shuffleMode, 1);

    ASSERT_TRUE(shuffleRegions[1].seqId == "seq1");
    ASSERT_EQ(shuffleRegions[1].start, 6);
    ASSERT_EQ(shuffleRegions[1].end, 9);
    ASSERT_EQ(shuffleRegions[1].shuffleMode, 3);

    ASSERT_TRUE(shuffleRegions[2].seqId == "seq1");
    ASSERT_EQ(shuffleRegions[2].start, 9);
    ASSERT_EQ(shuffleRegions[2].end, 12);
    ASSERT_EQ(shuffleRegions[2].shuffleMode, 1);
}

