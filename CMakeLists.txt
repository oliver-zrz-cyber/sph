cmake_minimum_required(VERSION 3.17)
project(RzSPH LANGUAGES CXX)
find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
set(POLYSCOPE_BACKEND_OPENGL3_GLFW "ON" CACHE BOOL "Enable openGL3_glfw backend")
set(POLYSCOPE_BACKEND_OPENGL_MOCK "ON" CACHE BOOL "Enable openGL_mock backend")

### Do anything needed for dependencies and bring their stuff in to scope
add_subdirectory(deps/polyscope)

set(US_SOURCES 
        	"src/main.cpp"
        	"src/func.cpp"
        	"src/func.h"
	"src/dfn.h"
        )
        
add_executable(RzSPH ${US_SOURCES})
set_target_properties(RzSPH PROPERTIES
                      CXX_STANDARD 17
                      CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(RzSPH PRIVATE /W4 /wd4201 /wd4840 /wd4100 /fp:fast /MP )
    set_target_properties(RzSPH PROPERTIES LINK_FLAGS /STACK:"100000000")

else()
    target_compile_options(RzSPH PRIVATE -Wall -Wextra -Wno-reorder -Wno-unused-parameter )
    set_target_properties(RzSPH PROPERTIES LINK_FLAGS -Wl,--stack, 100000000)
endif()

target_link_libraries(RzSPH  PRIVATE polyscope)
target_include_directories(RzSPH  
    PRIVATE "src/"
)