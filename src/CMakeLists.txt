#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -fabi-version=6 -fopenmp -fopenmp-simd")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -fopenmp-simd")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp-simd")

add_executable(readmapper readmapper.cpp 
                          #alignment.cpp
                          sequence.cpp 
                          fastareader.cpp
                          sw.cpp)

install (TARGETS readmapper DESTINATION bin)