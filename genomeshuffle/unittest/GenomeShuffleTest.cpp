//
// Created by 薄井光生 on 2017/10/11.
//

#define SEQAN_ENABLE_DEBUG 1

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include "geneticcode.h"
#include "shuffle_genome.h"
#include "myseqan.h"


using std::cout;
using std::cerr;
using std::endl;

TEST(myseqan, test_read_gff) {
    seqan::String<seqan::GffRecord> records;
    seqan::CharString gffFilepath = "../unittest/test.gff";
    read_gff(records, gffFilepath);

    ASSERT_EQ(seqan::length(records), 2);

    ASSERT_TRUE(records[0].ref == "seq1");
    ASSERT_TRUE(records[0].type == "gene");
    ASSERT_EQ(records[0].beginPos, 99); //seqan converts "1-based closed intervals" to "0-based half-open intervals"
    ASSERT_EQ(records[0].endPos, 200);
    ASSERT_EQ(records[0].strand, '+');
}

TEST(seqan, test_infix) {
    seqan::Dna5String seq1= "NNN";
    seqan::DnaString ifx1 = seqan::infix(seq1, 0, 3); //infix copy
    ASSERT_TRUE(seq1=="NNN"); //won't change the base string
    ASSERT_TRUE(ifx1=="AAA"); //N will be convert to NNN

    seqan::Dna5String seq2 = "NNN";
    seqan::Infix<seqan::Dna5String>::Type ifx2 = seqan::infix(seq2, 0, 3); //infix
    for (int i=0;i<3;i++){
        ifx2[i] = 'C';
    }
    ASSERT_TRUE(seq2=="CCC"); //will change the base string
    ASSERT_TRUE(ifx2=="CCC"); //N will be convert to NNN

    std::string seq3 = "ATG";
    seqan::DnaString ifx3 = seqan::infix(seq3, 0, 3);
    ASSERT_TRUE(ifx3=="ATG"); //infix can be used for std::string
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

TEST(genome_shuffle, test_classify) {
    seqan::GffRecord gff;
    seqan::DnaString seq;
    GeneticCode geneticCode = GeneticCode(11);
    MyCDS myCDS;

    gff = create_gffRecord("seq1", 0, 6, '+');
    seq = "ATGTAA";
    myCDS = MyCDS(gff, seq, geneticCode);
    ASSERT_EQ(myCDS.label, 0);

    gff = create_gffRecord("seq1", 0, 6, '-');
    seqan::reverseComplement(seq);
    myCDS = MyCDS(gff, seq, geneticCode);
    ASSERT_EQ(myCDS.label, 0);


    gff = create_gffRecord("seq1", 0, 6, '-');
    seq = "AAAAAA";
    myCDS = MyCDS(gff, seq, geneticCode);
    ASSERT_EQ(myCDS.label, 5);
}

TEST(genome_shuffle, test_get_shuffle_region) {
    seqan::String<seqan::GffRecord> gffs;
    seqan::appendValue(gffs, create_gffRecord("seq1", 3, 12, '+'));
    seqan::appendValue(gffs, create_gffRecord("seq2", 0, 12, '+'));
    seqan::appendValue(gffs, create_gffRecord("seq2", 7, 19, '+'));
    seqan::appendValue(gffs, create_gffRecord("seq3", 0, 24, '+'));
    seqan::appendValue(gffs, create_gffRecord("seq3", 8, 17, '+'));

    seqan::String<seqan::CharString> seqIds;
    seqan::appendValue(seqIds, "seq1");
    seqan::appendValue(seqIds, "seq2");
    seqan::appendValue(seqIds, "seq3");
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::appendValue(seqs, "NNNATGAAATAANNN"); //15
    seqan::appendValue(seqs, "ATGAAAAAATAAAAAATAA"); //19
    seqan::appendValue(seqs, "ATGAAAAAAAAAAATAAAAAATAA"); //24

    GeneticCode geneticCode = GeneticCode(11);
    std::vector< std::vector<MyCDS> > myCDS_vecvec;
    get_myCDS_vecvec(myCDS_vecvec, gffs, seqs, seqIds, geneticCode);

    int shuffleMode[2] = {1, 3};
    std::vector<ShuffleRegion> shuffleRegions;
    get_shuffle_region(shuffleRegions, shuffleMode, seqs, myCDS_vecvec);

    ASSERT_EQ(shuffleRegions.size(), 7);

    ASSERT_TRUE(shuffleRegions[0].seqIdx == 0);
    ASSERT_EQ(shuffleRegions[0].start, 0);
    ASSERT_EQ(shuffleRegions[0].end, 3);
    ASSERT_EQ(shuffleRegions[0].shuffleMode, 1);

    ASSERT_TRUE(shuffleRegions[1].seqIdx == 0);
    ASSERT_EQ(shuffleRegions[1].start, 6);
    ASSERT_EQ(shuffleRegions[1].end, 9);
    ASSERT_EQ(shuffleRegions[1].shuffleMode, 3);

    ASSERT_TRUE(shuffleRegions[2].seqIdx == 0);
    ASSERT_EQ(shuffleRegions[2].start, 12);
    ASSERT_EQ(shuffleRegions[2].end, 15);
    ASSERT_EQ(shuffleRegions[2].shuffleMode, 1);
}

TEST(genome_shuffle, test_shuffle_base){
    seqan::DnaString seq = "AAATTTGGGCCC";
    std::random_device rd;
    std::mt19937 mt(rd());

    shuffle_base(seq, mt);
    cout<<seq<<endl;
    ASSERT_TRUE(false);
}

TEST(genome_shuffle, test_shuffle_codon){
    seqan::DnaString seq = "AAATTTGGGCCC";
    std::random_device rd;
    std::mt19937 mt(rd());

    shuffle_codon(seq, mt);
    cout<<seq<<endl;
    ASSERT_TRUE(false);
}

TEST(genome_shuffle, test_shuffle_synonymous){
    seqan::DnaString seq = "TCTTCTTCCTCA";
    GeneticCode geneticCode = GeneticCode(11);
    geneticCode.update_count(seq);
    std::random_device rd;
    std::mt19937 mt(rd());

    shuffle_synonymous(seq, geneticCode, mt);
    cout<<seq<<endl;
    ASSERT_TRUE(geneticCode.translate(seq) == "SSSS");
    ASSERT_TRUE(false);
}
