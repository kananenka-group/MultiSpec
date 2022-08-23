cmake_minimum_required(VERSION 3.15)
project(MultiSpec 
        VERSION 1.3
        DESCRIPTION "MultiSpec is a set of tools for calculating vibrational spectra of condensed-phase systems."
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#########################################
# Make default build type as Release
########################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set (CMAKE_CXX_FLAGS "-Wall -Wextra")
set (CMAKE_CXX_FLAGS_DEBUG "-g")
set (CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")

#######################################
# Boost
#######################################
find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIR})
else()
   message(FATAL_ERROR "Boost not found")
endif()

#######################################
# BLAS
#######################################
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  include_directories(${BLAS_INCLUDE_DIR})
else()
  message(FATAL_ERROR "BLAS not found")
endif()

#####################################
#  OpenMP
#####################################
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###################################
# SOURCE
###################################
file(GLOB SOURCES_EXC "exciton/*.cpp" "time/*.cpp" "util/*.cpp" "fft/*.cpp" "math/*.cpp")
file(GLOB SOURCES_WAT "water/*.cpp" "util/*.cpp" "const/const.h" "vec/vecManip.h" "traj/*.cpp")
file(GLOB SOURCES_AMD "amideI/*.cpp" "util/*.cpp" "const/const.h" "vec/vecManip.h" "traj/*.cpp" "system/*.cpp")

include_directories(${PROJECT_SOURCE_DIR}/software/fftw3/build/include)
include_directories(${PROJECT_SOURCE_DIR}/software/xdrfile/build/include/xdrfile)

link_directories(${CMAKE_SOURCE_DIR}/software/xdrfile/build/lib)

########################################
#  ADD EXECUTABLES HERE
########################################
add_executable(excp ${SOURCES_EXC})
add_executable(water_gen ${SOURCES_WAT})
add_executable(amideI_gen ${SOURCES_AMD})

target_link_libraries(excp ${Boost_LIBRARIES}) #Boost::program_options)
target_link_libraries(excp ${BLAS_LIBRARIES})

target_link_libraries(water_gen ${Boost_LIBRARIES}) #Boost::program_options)
target_link_libraries(water_gen xdrfile)

target_link_libraries(amideI_gen ${Boost_LIBRARIES}) #Boost::program_options)
target_link_libraries(amideI_gen xdrfile)

