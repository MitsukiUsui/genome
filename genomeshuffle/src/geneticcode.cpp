#ifndef __H_GENETICCODE
#define __H_GENETICCODE

#include <iostream>
#include <cassert>
#include <random>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

class GeneticCode{
	string aa_str="*ACDEFGHIKLMNPQRSTVWY";
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
		char c=codon[2-i];
		switch (c){
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

	void clear_count(){/*{{{*/
		for(int i; i<64;i++){
			codonCount[i]=0;
			codonFreq[i]=0;
		}	
	}/*}}}*/

	//get sequence(multiple of 3) and update codonCount
	void update_count(seqan::Dna5String seq){
		assert(seqan::length(seq)%3==0);
		int aaLength = seqan::length(seq) / 3;
		for (int i = 0; i < aaLength; i++){
			seqan::Dna5String codon=seqan::infix(seq, 3*i, 3*(i+1));
			int codonId = codon_encode(codon);
			codonCount[codonId]++;
		}
	}

	//calculate codonFreq from codonCount
	void calc_freq(){
		for(int aaId=0;aaId<21;aaId++){
			//calculate count sum
			int countSum=0;
			for(int i=aaIndex[aaId];i<aaIndex[aaId+1];i++){
				countSum+=codonCount[aaToCodon[i]];
			}
			//calculate frequency
			for(int i=aaIndex[aaId];i<aaIndex[aaId+1];i++){
				int codonId=aaToCodon[i];
				if(countSum==0){//to avoid division by 0
					codonFreq[codonId]=0;
				}else{
					codonFreq[codonId]=codonCount[codonId]/countSum;
				}
			}
		}
	}

	//return randomly synonymous substituted codon, according to codon usage.
	seqan::Dna5String synonymous_sub(seqan::Dna5String codon, std::mt19937 &mt){
		int codonIdOld=codon_encode(codon);
		int aaId=codonToAa[codonIdOld];

		std::uniform_real_distribution<double> dist(0.0,1.0);
		double d=dist(mt);
		double cumSum=0;
		//calc codonIdNew
		int codonIdNew=-1;
		for(int i=aaIndex[aaId];i<aaIndex[aaId+1];i++){
			int codonId=aaToCodon[i];
			cumSum+=codonFreq[codonId];
			if(cumSum>d){
				codonIdNew=codonId;
			}
		}
		assert(codonIdNew!=-1);
		return codon_decode(codonIdNew);
	}

	seqan::String<seqan::AminoAcid> translate(seqan::Dna5String seq){
		assert(seqan::length(seq)%3==0);
		int aaLength=seqan::length(seq)/3;
		seqan::String<seqan::AminoAcid> seq_aa;
		seqan::resize(seq_aa, aaLength);
		for(int i=0;i<aaLength;i++){
			seqan::Dna5String codon=seqan::infix(seq, 3*i, 3*(i+1));
			seq_aa[i]=aa_decode(codonToAa[codon_encode(codon)]);	
		}
		return seq_aa;
	}

	void __debug(){
		cout<<"genetic code: "<<geneticCode_str<<endl;
		cout<<"aa: "<<aa_str<<endl;
		
		//output codon:aa pair
		for(int i=0;i<64;i++){
			cout<<codon_decode(i)<<"-"<<aa_decode(codonToAa[i])<<", ";
		}
		cout<<endl;

		//output aa:codonLst pair
		for(int i=0;i<21;i++){
			cout<<aa_decode(i)<<"-[";
			for(int j=aaIndex[i];j<aaIndex[i+1];j++){
				cout<<codon_decode(aaToCodon[j])<<",";
			}
			cout<<"], ";
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
