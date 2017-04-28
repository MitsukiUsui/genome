#ifndef MYUTIL_H_
#define MYUTIL_H_

#include <iostream>
#include <seqan/seq_io.h>

template <typename Tids, typename TSeqs, typename TFilepath>
void read_fasta(Tids & ids, TSeqs & seqs, TFilepath const & filepath){
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
}

#endif //MYUTIL_H_
