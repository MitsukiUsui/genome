project(unittest)

add_subdirectory(lib/googletest-release-1.8.0)
include_directories(lib/googletest-release-1.8.0/googletest/include)
include_directories(lib/googletest-release-1.8.0/googlemock/include)
include_directories(/Users/mitsuki/sandbox/genome/genomeshuffle/src)
include_directories (/Users/mitsuki/sandbox/genome/helper)

add_executable(runTest test.cpp)
target_link_libraries(runTest gtest gtest_main)
