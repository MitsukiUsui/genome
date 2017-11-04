project(src CXX)

# Add executable and link against SeqAn dependencies.
add_executable (shuffle_genome shuffle_genome.cpp)
target_link_libraries (shuffle_genome ${SEQAN_LIBRARIES})
