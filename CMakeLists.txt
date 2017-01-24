cmake_minimum_required(VERSION 2.8)
project(WVTICs)

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -O2 -Wall -fopenmp")

option(SAVE_WVT_STEPS "Save every relaxation step to file" ON)
option(SPH_CUBIC_SPLINE "Use cubic spline instead of wendland c6" OFF)
option(REJECTION_SAMPLING "Use von Neumann rejection sampling to improve initial random positions" ON)

if(SAVE_WVT_STEPS)
    add_definitions(-DSAVE_WVT_STEPS)
endif(SAVE_WVT_STEPS)

if(SPH_CUBIC_SPLINE)
    add_definitions(-DSPH_CUBIC_SPLINE)
endif(SPH_CUBIC_SPLINE)

if(REJECTION_SAMPLING)
    add_definitions(-DREJECTION_SAMPLING)
endif(REJECTION_SAMPLING)


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
    src/sph.c
    src/sph.h
    src/tree.c
    src/tree.h
    src/wvt_relax.c
    src/problems/sawtooth.c
    src/problems/tophat.c
    src/problems/sinewave.c
        src/problems/gradient.c)

add_executable(WVTICs ${SOURCE_FILES})

target_link_libraries(WVTICs -lm -lgsl -lgslcblas )