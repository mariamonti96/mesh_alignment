cmake_minimum_required(VERSION 3.9)
project(first_project)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

link_libraries(../openGL/OPENGL32)
link_libraries(../openGL/glut32)
link_libraries(../openGL/GLU32)
link_libraries(../openGL/GLAUX)
link_libraries(../glew_2_1_0/lib/Release/Win32/glew32s)
link_libraries(../glfw-3.2.1.bin.WIN32/lib-vc2015/glfw3)

#Header Locations

set(HEADER_DIRS
        tinyply/
        libigl/include/
        eigen-eigen-5a0156e40feb/
        glew_2_1_0/include/
        glfw-3.2.1/include/
        nanoflann-master/include/
        )


# Library
add_library(tinyply tinyply/tinyply.cpp)


#Executable
add_executable(first_project main.cpp calc_barycentre.cpp icp.cpp ICP_p2plane.cpp point_to_plane.cpp)
target_include_directories(first_project PRIVATE ${HEADER_DIRS})
target_link_libraries(first_project PRIVATE tinyply)