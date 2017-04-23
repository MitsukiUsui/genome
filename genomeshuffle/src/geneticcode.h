#ifndef __H_GENETICCODE
#define __H_GENETICCODE

#include <iostream>
#include <cassert>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using std::cout;
using std::endl;
using std::string;

class GeneticCode{
	int codonToAa[64]={};
	int codonCount[64]={};
	int codonFreq[64]={};
	int aaToCodon[64]={};
	int aaIndex[22]={};

int codon_encode(seqan::Dna5String codon){
	int codonId=0;
	return codonId;
}

seqan::Dna5String codon_decode(int codonId){
	seqan::Dna5String codon="ATG";
	return codon;
}

int aa_encode(char aa){
	int aaId=-1;
	for (int i=0;i<aa_str.length();i++){
		if (aa==aa_str[i]){
			aaId=i;
			break;
		}
	}
	if (aaId=-1){
		cerr<<"Undefined amino acid "<<aa<<endl;
		exit(1);
	}
	return aaId;
}

char aa_decode(){
	assert
	return 
}

public:
	string geneticCode_str;
	string aa_str="*ACDEFGHIKLMNPQRSTVWY";

	GeneticCode(string str){
		if(str.length()!=64){
			cerr<<"The length of genetic code is not 64."<<endl;
			exit(1);
		}
		geneticCode_str=str;

		//codonToAa
		for(int codonId;codonId<64;codonId++){
			codonToAa[codonId]=aa_encode(geneticCode_str[codonId]);
		}
		//codonCount,codonFreq
		clear_count();
		for (char aa:aa_str){
			aaToCodon
			for
		}

	}

	void clear_count(){
		for(int i; i<64;i++){
			codonCount[i]=0;
			codonFreq[i]=0;
		}	
	}


}

#endif /* __H_GENETICCODE/
