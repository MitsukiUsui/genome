#ifndef CDSS_H_
#define CDSS_H_

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <random>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>
#include <seqan/stream.h>

#include "geneticcode.h"

struct CDS{
	int startPos;
	int endPos;
	bool isForward;
	//------------------------------------------------------------
	//0...typical
	//1...shorter than 6
	//2...not multiple of 3
	//3...N included
	//4...stop codon inserted
	//5...stop codon not exists
	//------------------------------------------------------------
	int type=-1;

	//------------------------------------------------------------
	//basic constructor
	//------------------------------------------------------------
	CDS(seqan::GffRecord const & record);
	//------------------------------------------------------------
	//sentinal constructor
	//------------------------------------------------------------
	CDS(int pos);
	~CDS();
};


//============================================================
//Class to organize set of CDSs as vector of vector of CDS(cdss_vec)
//Each vector has set of CDSs from the same reference id.
//Last Element of each vector are sentinal, in order to possess reference length.
//============================================================
class CDSs{
	template<typename TSeq>
	int judge_type(TSeq & seq, 
				   GeneticCode const & gc);
	template<typename TSeqs>
	void set_types(TSeqs & seqs, 
	               GeneticCode const & gc);

public:
	std::vector< std::vector<CDS> > cdss_vec;

	CDSs(seqan::String<seqan::GffRecord> const & records, 
	     seqan::StringSet<seqan::CharString> const & ids, 
		 seqan::StringSet<seqan::Dna5String> & seqs, //seqs need not to be const, because of infix
		 GeneticCode const & gc);
	~CDSs();
	void __show();
};


//============================================================
//implementation of CDS
//============================================================

CDS::CDS(seqan::GffRecord const & record){/*{{{*/
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

CDS::CDS(int pos){/*{{{*/
	startPos=pos;
	endPos=pos;
	isForward=true;
}/*}}}*/

CDS::~CDS(){}


//============================================================
//implementation of CDSs
//============================================================

template <typename TSeq>//seqan::String or segment
int CDSs::judge_type(TSeq & seq, GeneticCode const & gc){/*{{{*/
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

template <typename TSeqs>
void CDSs::set_types(TSeqs & seqs, GeneticCode const & gc){/*{{{*/
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

CDSs::CDSs(seqan::String<seqan::GffRecord> const & records, seqan::StringSet<seqan::CharString> const & ids, seqan::StringSet<seqan::Dna5String> & seqs, GeneticCode const & gc){/*{{{*/
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

CDSs::~CDSs(){}
	
void CDSs::__show(){/*{{{*/
	for(int refIdx=0;refIdx<cdss_vec.size();refIdx++){
		//two counter	
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

#endif //CDSS_H_
