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
	if (aaId==-1){
		cerr<<"Undefined amino acid "<<aa<<endl;
		exit(1);
	}
	return aaId;
}

char aa_decode(){
	return '*';
}

public:
	string geneticCode_str;
	string aa_str="*ACDEFGHIKLMNPQRSTVWY";

	GeneticCode(string str){
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

	~GeneticCode(){}

	void clear_count(){
		for(int i; i<64;i++){
			codonCount[i]=0;
			codonFreq[i]=0;
		}	
	}
};

int main(){
	string str="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	GeneticCode gc(str);
	return 0;
}

#endif //* __H_GENETICCODE/
