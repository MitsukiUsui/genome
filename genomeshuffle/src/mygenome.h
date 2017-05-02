#ifndef MYGENOME_H_
#define MYGENOME_H_

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <random>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>
#include <seqan/stream.h>
#include <seqan/translation.h>


template <typename Tids, typename TSeqs, typename TFilepath>
void read_fasta(Tids & ids, TSeqs & seqs, TFilepath const & filepath){/*{{{*/
	seqan::SeqFileIn seqFileIn;
	if(!seqan::open(seqFileIn, seqan::toCString(filepath))){
		std::cerr<<"ERROR: Could not open "<<filepath<<std::endl;
		std::exit(1);
	}
	try{
		seqan::readRecords(ids,seqs,seqFileIn);
	}
	catch (seqan::Exception const & e){
		std::cout<<"ERROR: "<<e.what()<<std::endl;
		std::exit(2);
	}
	return;
}/*}}}*/



//modify: records...string of gff
//input: filepath...path to gff file
template <typename TFilepath>
void read_gff(seqan::String<seqan::GffRecord> & records,  TFilepath const & filepath){/*{{{*/
	seqan::GffFileIn gffIn;
	if(!seqan::open(gffIn,toCString(filepath))){
		std::cerr<<"ERROR: Could not open "<<filepath<<std::endl;
		std::exit(1);
	}

	seqan::GffRecord record;
	while(!seqan::atEnd(gffIn)){
		try{
			seqan::readRecord(record,gffIn);
		}
		catch(seqan::Exception const & e){
			std::cerr<<"ERROR: "<<e.what()<<std::endl;	
			continue;
		}
		seqan::appendValue(records,record);
	}
}/*}}}*/



//modify: records...string of gff (need to not to be const in order to use wrtieRecord function, but confusing)
//input: filepath...path to gff file
template <typename TFilepath>
void write_gff(seqan::String<seqan::GffRecord> & records,  TFilepath const & filepath){/*{{{*/
	seqan::GffFileOut gffOut;
	if(!seqan::open(gffOut,toCString(filepath))){
		std::cerr<<"ERROR: Could not open "<<filepath<<std::endl;
		std::exit(1);
	}

	for(unsigned i=0;i<seqan::length(records);i++){
		seqan::writeRecord(gffOut,records[i]);
	}

}/*}}}*/



class GeneticCode{
	std::string aa_str="*ACDEFGHIKLMNPQRSTVWY";//1 stop, 20 aa
	int codonToAa[64]={};
	int codonCount[64]={};
	double codonFreq[64]={};
	int aaToCodon[64]={};
	int aaIndex[22]={};//1 stop, 20 aa, 1 sentinal

	//if N included, return -1
	template <typename TSeq>
	int codon_encode(TSeq const & codon) const{/*{{{*/
		assert(seqan::length(codon)==3);
		
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
			case 'N':
				return -1;	
			default:
				std::cerr << "Unusual nucleotide " << codon[2-i] << " is included." << std::endl;
				exit(1);
			}
			base*=4;
		}

		return codonId;
	}/*}}}*/

	seqan::DnaString codon_decode(int codonId) const{/*{{{*/
		assert(codonId>=0 & codonId<64);

		int codon_int[3];
		seqan::DnaString codon;
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

	int aa_encode(char aa) const{/*{{{*/
		int aaId=-1;
		for (int i=0;i<aa_str.length();i++){
			if (aa==aa_str[i]){
				aaId=i;
				break;
			}
		}
		if (aaId==-1){
			std::cerr<<"Undefined amino acid "<<aa<<std::endl;
			exit(1);
		}
		return aaId;
	}/*}}}*/

	char aa_decode(int aaId) const{/*{{{*/
		assert(aaId>=0 & aaId<21);
		return aa_str[aaId];
	}/*}}}*/

public:
	std::string geneticCode_str;

	GeneticCode(std::string str){/*{{{*/
		if(str.length()!=64){
			std::cerr<<"The length of genetic code is not 64, but "<<str.length()<<std::endl;
			exit(1);
		}
		geneticCode_str=str;

		//codonToAa, count in aaBucket
		int aaBucket[21]={};
		for(int codonId=0; codonId<64;codonId++){
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
		//aaToCodon with incrementing bucket
		for(int codonId=0;codonId<64;codonId++){
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
	template <typename TSeq>//seqan::segment or seqan::String 
	void update_count(TSeq & seq){/*{{{*/
		assert(seqan::length(seq)%3==0);
		int aaLength = seqan::length(seq) / 3;
		for (int i = 0; i < aaLength; i++){
			seqan::Infix<seqan::Dna5String>::Type codon=seqan::infix(seq, 3*i, 3*(i+1));
			int codonId = codon_encode(codon);
			if(codonId!=-1)
				codonCount[codonId]++;
		}
	}/*}}}*/



	//calculate codonFreq from codonCount
	void calc_freq(){/*{{{*/
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
					codonFreq[codonId]=(double)codonCount[codonId]/countSum;
				}
			}
		}
	}/*}}}*/

	//randome inplace synonymous substitution, according to codon usage.
	template <typename TSeq>
	void synonymous_sub(TSeq codon, std::mt19937 &mt) const{/*{{{*/
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
		seqan::DnaString codonNew=codon_decode(codonIdNew);
		for(int i=0;i<3;i++){
			codon[i]=codonNew[i];	
		}
	}/*}}}*/

	
	//does not support ambiguous base for now
	seqan::String<seqan::AminoAcid> translate(seqan::Dna5String seq) const{/*{{{*/
		assert(seqan::length(seq)%3==0);
		int aaLength=seqan::length(seq)/3;
		seqan::String<seqan::AminoAcid> seq_aa;
		seqan::resize(seq_aa, aaLength);
		for(int i=0;i<aaLength;i++){
			seqan::Infix<seqan::Dna5String>::Type codon=seqan::infix(seq, 3*i, 3*(i+1));
			seq_aa[i]=aa_decode(codonToAa[codon_encode(codon)]);	
		}
		return seq_aa;
	}/*}}}*/

	template <typename TSeq>
	bool is_stop_codon(TSeq codon) const{/*{{{*/
		int aaId = codonToAa[codon_encode(codon)];
		if (aaId == 0){
			return true;
		}
		else{
			return false;
		}
	}/*}}}*/

	void __show() const{/*{{{*/
		std::cout<<"genetic code: "<<geneticCode_str<<std::endl;
		std::cout<<"aa: "<<aa_str<<std::endl;
		
		//output codon:aa pair
		for(int codonId=0; codonId<64; codonId++){
			std::cout<<codon_decode(codonId)<<"-"<<aa_decode(codonToAa[codonId])<<", ";
		}
		std::cout<<std::endl;

		//output aa:codonLst pair
		for(int aaId=0; aaId<21; aaId++){
			std::cout<<aa_decode(aaId)<<"-[";
			for(int i=aaIndex[aaId]; i<aaIndex[aaId+1];i++){
				std::cout<<codon_decode(aaToCodon[i])<<",";
			}
			std::cout<<"], ";
		}
		std::cout<<std::endl;
	}/*}}}*/

	//true...freq, false...count;
	void __show_freq(bool isFreq = true) const{/*{{{*/
		for (int aaId = 0; aaId < 21; aaId++){
			std::cout << aa_decode(aaId) << "\t";
		
			for(int i=aaIndex[aaId];i<aaIndex[aaId+1];i++){
				int codonId=aaToCodon[i];
				if(isFreq){
					std::cout<<codon_decode(codonId)<<": "<<codonFreq[codonId]<<", ";
				}
				else{
					std::cout<<codon_decode(codonId)<<": "<<codonCount[codonId]<<", ";
				}
			}
			std::cout << std::endl;
		}
	}/*}}}*/
};



class CDS{
public:
	int startPos;
	int endPos;
	bool isForward;
	int type=-1;

	CDS(seqan::GffRecord const & record){/*{{{*/
		startPos=record.beginPos;
		endPos=record.endPos;
		
		if(record.strand=='+'){
			isForward=true;
		}else if(record.strand=='-'){
			isForward=false;
		}else{
			std::cerr<<"ERROR: Undifined strand "<<record.strand<<std::endl;	
			exit(1);
		}
	}/*}}}*/

	//constractor for sentinal
	CDS(int pos){/*{{{*/
		startPos=pos;
		endPos=pos;
		isForward=true;
	}/*}}}*/

	~CDS(){/*{{{*/
	}/*}}}*/
};



class CDSs{
	//0...typical
	//1...shorter than 6
	//2...not multiple of 3
	//3...N included
	//4...stop codon inserted
	//5...stop codon not exists
	template <typename TSeq>//seqan::String or segment
	int judge_type(TSeq & seq, GeneticCode const & gc){/*{{{*/
		int length=seqan::length(seq);
		if(length < 6){
			return 1;
		}
		if((length%3) != 0){
			return 2;
		}
		for(int i=0;i<length;i++){
			if(seq[i]=='N'){
				return 3;
			}
		}
		int aaLength=length/3;
		for(int i=0;i<aaLength-1;i++){
			seqan::Infix<seqan::Dna5String>::Type codon=seqan::infix(seq,3*i,3*(i+1));//infix of infix is infix
			if(gc.is_stop_codon(codon))
				return 4;
		}
		
		seqan::Infix<seqan::Dna5String>::Type codon=seqan::infix(seq,3*(aaLength-1),3*aaLength);//infix of infix is infix
		if(!gc.is_stop_codon(codon))
			return 5;

		return 0;
	}/*}}}*/


	void set_types(seqan::StringSet<seqan::Dna5String> & seqs, GeneticCode const & gc){/*{{{*/
		for(int refIdx=0;refIdx<cdss_vec.size();refIdx++){
			for(auto itr=cdss_vec[refIdx].begin();itr!=cdss_vec[refIdx].end();itr++){
				if(itr->isForward){
					seqan::Infix<seqan::Dna5String>::Type sub=seqan::infix(seqs[refIdx], itr->startPos, itr->endPos);
					itr->type=judge_type(sub, gc);
				}
				else{
					seqan::Dna5String sub=seqan::infix(seqs[refIdx], itr->startPos, itr->endPos);
					seqan::reverseComplement(sub);
					itr->type=judge_type(sub, gc);
				}
			}
		}
		
	}/*}}}*/


public:
	std::vector< std::vector<CDS> > cdss_vec;

	//seqs need not to be const, because infix are used
	CDSs(seqan::String<seqan::GffRecord> const & records, seqan::StringSet<seqan::CharString> const & ids, seqan::StringSet<seqan::Dna5String> & seqs, GeneticCode const & gc){/*{{{*/
		cdss_vec.resize(seqan::length(ids));
		
		for (unsigned i=0;i<seqan::length(records);i++){
			if (records[i].type=="CDS"){
				//calc refIdx
				int refIdx=-1;
				seqan::CharString ref=records[i].ref;
				for(unsigned i=0;i<seqan::length(ids);i++){
					if(ids[i]==ref){
						refIdx=i;
						break;
					}
				}
				//add cds	
				if(refIdx!=-1){
					CDS c(records[i]);
					cdss_vec[refIdx].push_back(c);	
				}
			}
		}

		set_types(seqs, gc);
		
		//add sentinal
		for(int refIdx=0;refIdx<cdss_vec.size();refIdx++){
			CDS sentinal(seqan::length(seqs[refIdx]));
			cdss_vec[refIdx].push_back(sentinal);
		}

	}/*}}}*/


	~CDSs(){/*{{{*/
	}/*}}}*/
	

	void __show(){/*{{{*/
		//two counter	

		//count
		for(int refIdx=0;refIdx<cdss_vec.size();refIdx++){
			std::vector<int> strandCount(3,0);//backward, forward, sentinal
			std::vector<int> typeCount(6,0);//6 type 

			for(auto itr=cdss_vec[refIdx].begin();itr!=cdss_vec[refIdx].end();itr++){
				if(itr->type!=-1){
					strandCount[itr->isForward]++;
					typeCount[itr->type]++;
				}else{
					strandCount[2]++;
				}
			}	
			
			//output
			std::cout<<"REFIDX: "<<refIdx<<std::endl;
			std::cout<<"\tforward: "<<strandCount[1]<<std::endl;
			std::cout<<"\tbackward: "<<strandCount[0]<<std::endl;
			std::cout<<"\tsentinal: "<<strandCount[2]<<std::endl;
			std::cout<<std::endl;

			for(unsigned i=0;i<6;i++){
				std::cout<<"\ttype"<<i<<": "<<typeCount[i]<<std::endl;
			}
		}

	}/*}}}*/
};


#endif //MYGENOME_H_
