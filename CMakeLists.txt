cmake_minimum_required(VERSION 3.5)
project(MultiSpec VERSION 1.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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
file(GLOB SOURCES_EXC "src/exciton/*.cpp" "src/time/*.cpp" "src/util/*.cpp" "src/fft/*.cpp" "src/math/*.cpp")
file(GLOB SOURCES_WAT "src/water/*.cpp" "src/util/*.cpp" "src/const/const.h" "src/vec/vecManip.h" "src/traj/*.cpp")

include_directories(${PROJECT_SOURCE_DIR}/src/software/fftw3/build/include)
include_directories(${PROJECT_SOURCE_DIR}/src/software/xdrfile/build/include/xdrfile)

link_directories(${CMAKE_SOURCE_DIR}/src/software/xdrfile/build/lib)

########################################
#  ADD EXECUTABLES HERE
########################################
add_executable(excp ${SOURCES_EXC})
add_executable(water_gen ${SOURCES_WAT})

target_link_libraries(excp ${Boost_LIBRARIES}) #Boost::program_options)
target_link_libraries(excp ${BLAS_LIBRARIES})
target_link_libraries(water_gen ${Boost_LIBRARIES}) #Boost::program_options)
target_link_libraries(water_gen xdrfile)

