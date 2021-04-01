cmake_minimum_required(VERSION 3.1)
project(DevelopableStitching)

set(SOURCE_DIR "./src/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(INCLUDES "${SOURCE_DIR}/../ext/gco-v3.0/")


# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"    ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
FILE(GLOB SRCFILES ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.hpp ${SOURCE_DIR}/*.c ${SOURCE_DIR}/*.h)

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDES})

FILE(GLOB GCOSRCFILES ${SOURCE_DIR}/../ext/gco-v3.0/*.cpp)
add_library(gco STATIC ${GCOSRCFILES})
target_compile_definitions(gco PRIVATE -DGCO_ENERGYTYPE="long long")

target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl_glfw_imgui gco)