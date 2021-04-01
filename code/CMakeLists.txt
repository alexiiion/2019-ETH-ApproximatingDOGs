cmake_minimum_required(VERSION 3.1)
project(ApproximatingWithDOGs)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# set(DEPENDENCIES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies")
set(APPROX_DEPENDENCIES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies" CACHE PATH "dependencies directory")
# message("DEPENDENCIES_DIR: " ${DEPENDENCIES_DIR})


# libigl settings
option(LIBIGL_USE_STATIC_LIBRARY 	 "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            ON)



# Include libigl
set(APPROX_LIBIGL_DIR "${APPROX_DEPENDENCIES_DIR}/libigl" CACHE PATH "libigl directory")
list(APPEND CMAKE_MODULE_PATH "${APPROX_LIBIGL_DIR}/cmake")
include(libigl)



# Include fitpack spline library
set(APPROX_FITPACK_DIR "${APPROX_DEPENDENCIES_DIR}/fitpack" CACHE PATH "fitpack (splines) directory")

include_directories(${APPROX_FITPACK_DIR}/fitpackpp)
link_directories(${APPROX_FITPACK_DIR}/binaries_win_x64)
link_libraries(fitpackpp.lib fitpack.lib)

set(APPROX_FORTRAN_COMPILER_DIR "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows" CACHE PATH "directory where the fortan compiler is found")
link_directories(
	"${APPROX_FORTRAN_COMPILER_DIR}/compiler/lib/intel64"
	"${APPROX_FORTRAN_COMPILER_DIR}/redist/intel64_win/compiler" 
)



# Include multi-label graph cut library
set(APPROX_GRAPHCUT_DIR "${APPROX_DEPENDENCIES_DIR}/gco-v3.0" CACHE PATH "graph cut directory")
include_directories(${APPROX_GRAPHCUT_DIR})
list(APPEND APPROX_GRAPHTCUT_SRC 
	${APPROX_GRAPHCUT_DIR}/GCoptimization.cpp
	${APPROX_GRAPHCUT_DIR}/graph.cpp
	${APPROX_GRAPHCUT_DIR}/LinkedBlockList.cpp
	${APPROX_GRAPHCUT_DIR}/maxflow.cpp
)
message("including graph cut dir: ${APPROX_GRAPHCUT_DIR}")
message("adding graph cut src: ${APPROX_GRAPHTCUT_SRC}")

foreach(GRAPHTCUTSRC_FILE IN ITEMS ${APPROX_GRAPHTCUT_SRC})
	source_group("__External Sources\\GraphCutOptimization" FILES ${GRAPHTCUTSRC_FILE})
endforeach()



# Add Pardiso solver
set(APPROX_PARDISO_DIR "${APPROX_DEPENDENCIES_DIR}/pardiso-win/" CACHE PATH "directory where pardiso lib is found")
set(APPROX_PARDISO_LIB "libpardiso600-WIN-X86-64" CACHE STRING "name of pardiso lib (without extension)")
find_library(APPROX_PARDISO ${APPROX_PARDISO_LIB} PATHS ${APPROX_PARDISO_DIR})

if (${APPROX_PARDISO_LIB} MATCHES "NOTFOUND")
	message("ERROR: Could not find Pardiso solver.\n")
else()
	message("Found PARDISO library at ${APPROX_PARDISO}")
	link_directories(${APPROX_PARDISO_DIR})
	link_libraries(${APPROX_PARDISO_LIB})
endif()




# Add source code such that it mirrors the structure on the file system in the project

set(SRC_ROOT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/")
# message("SRC_ROOT_PATH: " ${SRC_ROOT_PATH})


# Add header files recursively
file(
    GLOB_RECURSE HEADER_FILES 
    LIST_DIRECTORIES true
    "${SRC_ROOT_PATH}/*.h*"
)

foreach(H_FOLDER IN ITEMS ${HEADER_FILES})
    get_filename_component(H_FOLDER_PATH "${H_FOLDER}" PATH)
    #file(RELATIVE_PATH H_FOLDER_PATH_RELATIVE "${SRC_ROOT_PATH}" "${H_FOLDER_PATH}")
    #string(REPLACE "/" "\\" ESCAPED_FOLDER_PATH "${H_FOLDER_PATH_RELATIVE}")
    string(REPLACE "/" "\\" ESCAPED_FOLDER_PATH "${H_FOLDER_PATH}")
	include_directories(${ESCAPED_FOLDER_PATH})
endforeach()


# Add source files recursively
file(
    GLOB_RECURSE SRCFILES 
    LIST_DIRECTORIES false
    "${SRC_ROOT_PATH}/*.c*"
    "${SRC_ROOT_PATH}/*.h*"
)

# Group source files for IDE (Visual Studio) -- remove this if it causes problems for other IDEs
foreach(SRC_FILE IN ITEMS ${SRCFILES})
    get_filename_component(SRC_FILE_PATH "${SRC_FILE}" PATH)
    file(RELATIVE_PATH SRC_PATH_RELATIVE "${SRC_ROOT_PATH}" "${SRC_FILE_PATH}")
    string(REPLACE "/" "\\" GROUP_PATH "${SRC_PATH_RELATIVE}")
    source_group("${GROUP_PATH}" FILES "${SRC_FILE}")
endforeach()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -MD")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -MD")

add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${APPROX_GRAPHTCUT_SRC})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui igl::xml)#







############################################
# TRIED adding FITPACK with find_library() #
############################################



# set(FITPACK_DIR ${DEPENDENCIES_DIR}/fitpack)
# set(FITPACK_BIN ${FITPACK_DIR}/binaries_win_x64/${Configuration})
# message("FITPACK_BIN: ${FITPACK_BIN}")
# set(FITPACK_LIB "fitpack.lib")
# set(FITPACKPP_LIB "fitpackpp.lib")
# find_library(FITPACK ${FITPACK_LIB} PATHS ${FITPACK_BIN})
# find_library(FITPACKPP ${FITPACKPP_LIB} PATHS ${FITPACK_BIN})

# if (${FITPACK_LIB} MATCHES "NOTFOUND" OR ${FITPACKPP_LIB} MATCHES "NOTFOUND")
	# message("ERROR: Could not find FITPACK lib.")
	# message("       FITPACK_LIB: ${FITPACK_LIB}")
	# message("       FITPACKPP_LIB: ${FITPACKPP_LIB}")
# else()
	# message("Found FITPACK library at ${FITPACK} and ${FITPACKPP}")
	
	# include_directories(${FITPACK_DIR}/fitpackpp)
	# link_directories(${FITPACK_BIN})
	# link_libraries(${FITPACK_LIB} ${FITPACKPP_LIB})

	# ## splines are using a fortran compiler. link your fortran compiler here (I use Intel)
	# set(FORTRAN_COMPILER_DIR "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows")
	# link_directories(
		# "${FORTRAN_COMPILER_DIR}/compiler/lib/intel64"
		# "${FORTRAN_COMPILER_DIR}/redist/intel64_win/compiler" 
	# )
	# #C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\compiler\lib\intel64
	# #C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64_win\compiler
# endif()
