cmake_minimum_required(VERSION 3.0)
project(WVTICs)

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -O2 -Wall -fopenmp")

set(SOURCE_FILES
    src/problems/constant_density.c
    src/problems/magneticum.c
    src/aux.c
    src/globals.h
    src/ids.c
    src/io.c
    src/io.h
    src/macro.h
    src/main.c
    src/peano.c
    src/peano.h
    src/positions.c
    src/proto.h
    src/setup.c
    src/sort.c
    src/sort.h
    src/sph.c
    src/sph.h
    src/tree.c
    src/tree.h
    src/unit.c
    src/wvt_relax.c)

add_executable(WVTICs ${SOURCE_FILES})

target_link_libraries(WVTICs -lm -lgsl -lgslcblas )