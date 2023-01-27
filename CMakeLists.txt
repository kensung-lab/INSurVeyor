cmake_minimum_required(VERSION 2.8)
project(SurVeyor-INS)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(libs)

find_library(HTS_LIB hts PATHS htslib-1.13/lib NO_DEFAULT_PATH)
include_directories(htslib-1.13/include)
if (NOT HTS_LIB)
	message(STATUS "Included htslib was not built. Assuming htslib 1.13 or higher is installed system-wise.")
	find_library(HTS_LIB hts)
endif()

add_executable(reads_categorizer reads_categorizer.cpp)

add_executable(clip_consensus_builder clip_consensus_builder.cpp)

add_executable(call_insertions libs/ssw.c libs/ssw_cpp.cpp call_insertions.cpp)

add_executable(dc_remapper libs/ssw.c libs/ssw_cpp.cpp dc_remapper.cpp)

add_executable(add_filtering_info add_filtering_info.cpp)

add_executable(filter filter.cpp)

add_executable(normalise normalise.cpp)

target_link_libraries(reads_categorizer "${HTS_LIB}")
target_link_libraries(clip_consensus_builder "${HTS_LIB}")
target_link_libraries(call_insertions "${HTS_LIB}")
target_link_libraries(dc_remapper "${HTS_LIB}")
target_link_libraries(add_filtering_info "${HTS_LIB}")
target_link_libraries(filter "${HTS_LIB}")
target_link_libraries(normalise "${HTS_LIB}")
