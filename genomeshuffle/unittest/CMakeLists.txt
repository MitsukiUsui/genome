project(unittest)

include_directories(../src)
add_subdirectory(lib/googletest-release-1.8.0)

add_executable(runGenomeShuffleTest GenomeShuffleTest.cpp)
target_link_libraries(runGenomeShuffleTest gtest gtest_main)

add_executable(runGeneticCodeTest GeneticCodeTest.cpp)
target_link_libraries(runGeneticCodeTest gtest gtest_main)
