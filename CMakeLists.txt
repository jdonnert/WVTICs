cmake_minimum_required(VERSION 2.8)
project(WVTICs)

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -O2 -Wall -fopenmp")

option(SAVE_WVT_STEPS "Save every relaxation step to file" ON)
option(SPH_CUBIC_SPLINE "Use cubic spline instead of wendland c6" OFF)
option(SPH_WC2 "Use wendland c2 instead of wendland c6" OFF)
option(REJECTION_SAMPLING "Use von Neumann rejection sampling to improve initial random positions" ON)
option(PEANO_SAMPLING "Use peano curve based sampling to improve initial random positions" OFF)
option(EAT_PNG "Possibility to use .png files as input for density" ON)
option(TWO_DIM "Set the code to do 2D initial conditions instead of 3D -> z component is set to 0" OFF)
option(BRUTE_FORCE_NGB "Use a brute force neighbour finder instead of the tree based one" OFF)
option(OUTPUT_DIAGNOSTICS "Output extra diagnostics to file" ON)

if(SAVE_WVT_STEPS)
    add_definitions(-DSAVE_WVT_STEPS)
endif(SAVE_WVT_STEPS)

if(SPH_CUBIC_SPLINE)
    add_definitions(-DSPH_CUBIC_SPLINE)
endif(SPH_CUBIC_SPLINE)

if(SPH_WC2)
    add_definitions(-DSPH_WC2)
endif(SPH_WC2)

if(PEANO_SAMPLING)
    add_definitions(-DPEANO_SAMPLING)
endif(PEANO_SAMPLING)

if(REJECTION_SAMPLING)
    add_definitions(-DREJECTION_SAMPLING)
endif(REJECTION_SAMPLING)

if(EAT_PNG)
    add_definitions(-DEAT_PNG)

    find_package(PNG REQUIRED)
    add_definitions(${PNG_DEFINITIONS})
    include_directories(${PNG_INCLUDE_DIRS})
endif(EAT_PNG)

if(TWO_DIM)
    add_definitions(-DTWO_DIM)
endif(TWO_DIM)

if(BRUTE_FORCE_NGB)
    add_definitions(-DBRUTE_FORCE_NGB)
endif(BRUTE_FORCE_NGB)

if(OUTPUT_DIAGNOSTICS)
    add_definitions(-DOUTPUT_DIAGNOSTICS)
endif(OUTPUT_DIAGNOSTICS)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

get_directory_property( DirDefs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
message(STATUS "Compiling with ${CMAKE_C_FLAGS}\n   and definitions ${DirDefs}")


set(SOURCE_FILES
        src/aux.c
        src/bias_correction.c
        src/ids.c
        src/io.c
        src/main.c
        src/peano.c
        src/positions.c
        src/setup.c
        src/sph.c
        src/tree.c
        src/wvt_relax.c
        src/external/readpng.c
        src/external/morton_utils.c
        src/peanowalk.c
    	src/kernel.c
    	src/diagnostics.c
    	src/redistribution.c
        src/problems/constant_density.c
        src/problems/magneticum.c
        src/problems/png.c
        src/problems/sawtooth.c
        src/problems/tophat.c
        src/problems/sinewave.c
        src/problems/gradient.c
        src/problems/doubleshock.c
        src/problems/user.c
        src/problems/sodshock.c
        src/problems/sedov.c
        src/problems/kelvin_helmholtz.c
        src/problems/box.c
        src/problems/gresho.c
        src/problems/orszag_tang.c
        src/problems/rotor.c
        src/problems/alfven.c
        src/problems/boss.c
        src/problems/galaxy_cluster.c
        src/problems/rayleigh_taylor.c
        src/problems/ryu_jones.c
        src/problems/strong_blast.c
        src/problems/blob.c
        src/problems/disk.c
        src/problems/evrard.c
        src/problems/zeldovich_pancake.c
        src/problems/keplerian_ring.c)

add_executable(WVTICs ${SOURCE_FILES})

target_link_libraries(WVTICs -lm ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )

if(EAT_PNG)
    target_link_libraries(WVTICs ${PNG_LIBRARIES})
endif(EAT_PNG)
