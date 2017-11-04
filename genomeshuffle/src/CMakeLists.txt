project(src)

add_executable (shuffle_genome shuffle_genome.cpp)
target_link_libraries (shuffle_genome ${SEQAN_LIBRARIES})
