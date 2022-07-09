cmake_minimum_required(VERSION 3.13)
project(MultiSpec)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
find_package(BLAS REQUIRED)

file(GLOB SOURCES_EXC "src/exciton/*.cpp" "src/time/*.cpp" "src/util/*.cpp" "src/fft/*.cpp" "src/math/*.cpp")
file(GLOB SOURCES_WAT "src/water/*.cpp" "src/util/*.cpp" "src/const/const.h" "src/vec/vecManip.h" "src/traj/*.cpp")

include_directories(${BLAS_INCLUDE_DIR})
include_directories(/work/akanane/sw/MultiSpec/src/software/fftw3/build/include)
include_directories(/work/akanane/sw/MultiSpec/src/software/xdrfile/build/include/xdrfile)

add_executable(excp ${SOURCES_EXC})
add_executable(water_gen ${SOURCES_WAT})

target_link_libraries(excp Boost::program_options)
target_link_libraries(excp ${BLAS_LIBRARIES})

target_link_libraries(water_gen Boost::program_options)
target_link_libraries(water_gen /work/akanane/sw/MultiSpec/src/software/xdrfile/build/lib/libxdrfile.a)

