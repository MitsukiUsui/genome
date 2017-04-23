#ifndef __H_GENETICCODE
#define __H_GENETICCODE

#include <iostream>
#include <cassert>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

class GeneticCode{
	int codonToAa[64]={};
	int codonCount[64]={};
	double codonFreq[64]={};
	int aaToCodon[64]={};
	int aaIndex[22]={};//20 aa, 1 stop, 1 sentinal

int codon_encode(seqan::Dna5String codon){/*{{{*/
	if (seqan::length(codon) != 3){
		cerr << "The length of codon is not 3." << endl;
		exit(1);
	}

	int codonId=0;
	int base=1;
	for (int i = 0; i < 3; i++){
		switch (codon[2-i]){
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
			cerr << "Unusual nucleotide " << codon[2-i] << " is included." << endl;
			exit(1);
		}
		base*=4;
	}
	return codonId;
}/*}}}*/

seqan::Dna5String codon_decode(int codonId){/*{{{*/
	assert(codonId>=0 & codonId<64);

	int codon_int[3];
	seqan::Dna5String codon;
	seqan::resize(codon,3);

	codon_int[0] = codonId / 16;
	codon_int[1] = (codonId % 16) / 4;
	codon_int[2] = (codonId % 16) % 4;

	for (int i = 0; i < 3; i++){
		switch (codon_int[i]){
		case 0:
			codon[i] = 'T';
			break;
		case 1:
			codon[i] = 'C';
			break;
		case 2:
			codon[i] = 'A';
			break;
		case 3:
			codon[i] = 'G';
			break;
		}
	}
	return codon;
}/*}}}*/

int aa_encode(char aa){/*{{{*/
	int aaId=-1;
	for (int i=0;i<aa_str.length();i++){
		if (aa==aa_str[i]){
			aaId=i;
			break;
		}
	}
	if (aaId==-1){
		cerr<<"Undefined amino acid "<<aa<<endl;
		exit(1);
	}
	return aaId;
}/*}}}*/

char aa_decode(int aaId){/*{{{*/
	assert(aaId>=0 & aaId<21);
	return aa_str[aaId];
}/*}}}*/

public:
	string geneticCode_str;
	string aa_str="*ACDEFGHIKLMNPQRSTVWY";

	GeneticCode(string str){/*{{{*/
		if(str.length()!=64){
			cerr<<"The length of genetic code is not 64, but "<<str.length()<<endl;
			exit(1);
		}
		geneticCode_str=str;

		//codonToAa, count in aaBucket
		int aaBucket[21]={};
		for(int codonId;codonId<64;codonId++){
			int aaId=aa_encode(geneticCode_str[codonId]);
			codonToAa[codonId]=aaId;
			aaBucket[aaId]++;
		}
		//convert bucket to index
		aaIndex[0]=0;
		for(int i=1;i<22;i++){
			aaIndex[i]=aaIndex[i-1]+aaBucket[i-1];
		}
		assert(aaIndex[21]==64);
		//copy index to bucket
		for(int i=0;i<21;i++){
			aaBucket[i]=aaIndex[i];
		}
		//aaToCodon
		for(int codonId;codonId<64;codonId++){
			int aaId=aa_encode(geneticCode_str[codonId]);
			aaToCodon[aaBucket[aaId]++]=codonId;
		}
	}/*}}}*/

	~GeneticCode(){}

	void clear_count(){
		for(int i; i<64;i++){
			codonCount[i]=0;
			codonFreq[i]=0;
		}	
	}

	void update_count(seqan::Dna5String seq){
		assert(seqan::length(seq)%3==0);
		int loopCount = seqan::length(seq) / 3;
		for (int i = 0; i < loopCount; i++){
			seqan::Dna5String codon=seqan::infix(seq, 3*i, 3*(i+1));
			int codonId = codon_encode(codon);
			codonCount[codonId]++;
		}
	}

	void calc_freq(){
	}

	seqan::Dna5String synonymous_sub(seqan::Dna5String codon, mt19937 &mt){
		seqan::Dna5String retVal="ATG";
		return retVal;
	}
	seqan::String<seqan::AminoAcid> translate(seqan::Dna5String seq){
		seqan::String<seqan::AminoAcid> retVal="*";
		return retVal;
	}


	void __debug(){
		cout<<"genetic code: "<<geneticCode_str<<endl;
		cout<<"aa: ";
		for(int i=0;i<21;i++){
			cout<<i<<"-"<<aa_str[i]<<", ";
		}
		cout<<endl;
		for(int i=0;i<64;i++){
			cout<<codonToAa[i]<<", "; 
		}
		cout<<endl;
		for(int i=0;i<64;i++){
			cout<<aaToCodon[i]<<", ";
		}
		cout<<endl;
		for(int i=0;i<21;i++){
			cout<<aaIndex[i]<<", ";
		}
		cout<<endl;
	}
};

int main(){
	string str="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	GeneticCode gc(str);
	gc.__debug();
	return 0;
}

#endif //* __H_GENETICCODE/
