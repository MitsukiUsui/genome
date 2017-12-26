# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

SET(CMAKE_CXX_COMPILER /usr/local/bin/g++)
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build)

# Name of project and that it is C++ only.
project (shuffle_genome CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

#Add additional include directories.
include_directories(../helper)
include_directories(unittest/lib/googletest-release-1.8.0/googletest/include)
include_directories(unittest/lib/googletest-release-1.8.0/googlemock/include)

add_subdirectory(unittest)

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})
# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

add_subdirectory(src)
