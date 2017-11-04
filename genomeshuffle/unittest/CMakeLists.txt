project(unittest)

add_subdirectory(lib/googletest-release-1.8.0)
add_executable(runGenomeShuffleTest GenomeShuffleTest.cpp)
target_link_libraries(runGenomeShuffleTest gtest gtest_main ${SEQAN_LIBRARIES})
#add_executable(runGeneticCodeTest GeneticCodeTest.cpp)
#target_link_libraries(runGenomeShuffleTest gtest gtest_main ${SEQAN_LIBRARIES})

