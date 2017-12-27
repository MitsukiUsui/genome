project(src)

add_executable (fasta_split fasta_split.cpp)
target_link_libraries (fasta_split ${SEQAN_LIBRARIES})
