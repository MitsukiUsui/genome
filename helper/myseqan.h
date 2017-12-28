#ifndef MYSEQAN_H_
#define MYSEQAN_H_

#include <iostream>
#include <cassert>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>
#include <seqan/stream.h>


//------------------------------------------------------------
//input:  filepath (path to fasta file)
//modify: ids (sequence ids), seqs(sequences)
//------------------------------------------------------------
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

//------------------------------------------------------------
//input:  filepath (path to output file)
//modify: ids (sequence ids), seqs(sequences)
//------------------------------------------------------------
template <typename Tids, typename TSeqs, typename TFilepath>
void write_fasta(Tids const & ids, TSeqs const & seqs, TFilepath const & filepath){/*{{{*/
	seqan::SeqFileOut seqFileOut;
	if(!seqan::open(seqFileOut, seqan::toCString(filepath))){
		std::cerr<<"ERROR: Could not open "<<filepath<<std::endl;
		std::exit(1);
	}
	try{
		seqan::writeRecords(seqFileOut, ids, seqs);
	}
	catch (seqan::Exception const & e){
		std::cout<<"ERROR: "<<e.what()<<std::endl;
		std::exit(2);
	}
	return;
}/*}}}*/

//------------------------------------------------------------
//input:  filepath (path to gff file)
//modify: records (set of GffRecord)
//------------------------------------------------------------
template <typename TFilepath>//*** adopt to template (String or StringSet)
void read_gff(seqan::String<seqan::GffRecord> & records,  TFilepath const & filepath){/*{{{*/
	seqan::GffFileIn gffIn;
	if(!seqan::open(gffIn,seqan::toCString(filepath))){
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


//------------------------------------------------------------
//input:  filepath (path to gff file)
//modify: records (set of GffRecord)
//	* records need not to be const in order to use writeRecord funcion, confusing....
//------------------------------------------------------------
template <typename TFilepath>//*** adopt to template 
void write_gff(seqan::String<seqan::GffRecord> & records,  TFilepath const & filepath){/*{{{*/
	seqan::GffFileOut gffOut;
	if(!seqan::open(gffOut,seqan::toCString(filepath))){
		std::cerr<<"ERROR: Could not open "<<filepath<<std::endl;
		std::exit(1);
	}

	for(unsigned i=0;i<seqan::length(records);i++){
		seqan::writeRecord(gffOut,records[i]);
	}

}/*}}}*/


#endif //MYSEQAN_H_
