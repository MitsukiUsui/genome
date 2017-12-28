project(src)

add_executable (fasta_split fasta_split.cpp)
target_link_libraries (fasta_split ${SEQAN_LIBRARIES})

add_executable (classify_cds classify_cds.cpp)
target_link_libraries (classify_cds ${SEQAN_LIBRARIES})
