#ifndef CDSS_H_
#define CDSS_H_

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <random>
#include <algorithm>

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
	//-1...sentinel
	//0...typical
	//1...shorter than 6
	//2...not multiple of 3
	//3...spanning boundary
	//4...N included
	//5...stop codon inserted
	//6...stop codon not exists
	//------------------------------------------------------------
	int type=-1;

	//------------------------------------------------------------
	//comparison operator
	//------------------------------------------------------------
	friend bool operator < (CDS const & c1, CDS const & c2){
		return c1.startPos < c2.startPos;
	}

	//------------------------------------------------------------
	//constructors
	//------------------------------------------------------------
	CDS(seqan::GffRecord const & record);//basic
	CDS(int pos);//sentinel
	~CDS();
};


//============================================================
//Class to organize set of CDSs as vector of vector of CDS(cdss_vec)
//Each vector has set of CDSs from the same reference id.
//Last Element of each vector are sentinal, in order to possess reference length.
//============================================================
class CDSs{
	template<typename TItr, typename TSeq>
	int judge_type(TItr itr, TSeq & seq, GeneticCode const & gc);
	
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

template <typename TItr, typename TSeq>//seqan::String or segment
int CDSs::judge_type(TItr itr, TSeq & seq, GeneticCode const & gc){/*{{{*/
	assert(itr->startPos < itr->endPos);

	//3 position level classification
	if( itr->startPos < 0 | itr->endPos > seqan::length(seq)){
		return 1;
	}
	int length=(itr->endPos)-(itr->startPos);
	if(length < 6){
		return 2;
	}
	if((length%3) != 0){
		return 3;
	}

	//get subseq
	seqan::Dna5String subseq=seqan::infix(seq, itr->startPos, itr->endPos);
	if(!itr->isForward){
		seqan::reverseComplement(subseq);
	}
	
	for(int i=0;i<length;i++){
		if(subseq[i]=='N'){
			return 4;
		}
	}
	int aaLength=length/3;
	for(int i=0;i<aaLength-1;i++){
		seqan::Infix<seqan::Dna5String>::Type codon=seqan::infix(subseq,3*i,3*(i+1));//infix of infix is infix
		if(gc.is_stop_codon(codon))
			return 5;
	}
	
	seqan::Infix<seqan::Dna5String>::Type codon=seqan::infix(subseq,3*(aaLength-1),3*aaLength);//infix of infix is infix
	if(!gc.is_stop_codon(codon))
		return 6;

	return 0;
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

	//set types
	for(int refIdx=0; refIdx<cdss_vec.size();refIdx++){
		for(auto itr=cdss_vec[refIdx].begin(); itr!=cdss_vec[refIdx].end(); itr++){
			itr->type=judge_type(itr, seqs[refIdx], gc);
		}
	}

	//sort vector
	for(int refIdx=0;refIdx<cdss_vec.size();refIdx++){
		std::sort(cdss_vec[refIdx].begin(), cdss_vec[refIdx].end());
	}

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
		std::vector<int> typeCount(7,0);//7 type 

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

		for(int i=0;i<7;i++){
			std::cout<<"\ttype"<<i<<": "<<typeCount[i]<<std::endl;
		}
	}

}/*}}}*/

#endif //CDSS_H_
