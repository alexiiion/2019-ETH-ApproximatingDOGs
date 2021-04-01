## How to build

This code was only tested on **Windows** 10, 64-bit with MSVC build system, using Visual Studio 2017/2019.

This project uses CMake as a build system and is (almost) completley automated. 
Run cmake in this folder, i.e. where CMakeLists.txt is located.
I use CMake GUI and the following settings:
- Generator: Visual Studio 16 2019  
- Platform for the generator: x64
- 'Use default native compilers'
The pardiso solver lib (`libpardiso600-WIN-X86-64.dll`) still has to be copied manually to the build folder where the .exe will be generated.

In case of problems, I stored all the relevant Visual Studio project settings in `/__project-config-WIN/## current project properties.txt`

## Dependencies
Used libraries are located in the folder `/dependencies/`
Some binaries are compiled for Windows. To build this project on other platforms, the dependencies need to be rebuild and the CMake file needs to be updated to point to them. All platform-dependent libraries are marked accordingly and the CMakeLists.txt file is commented.

Additionally, for math libraries, this project has been used with the **Intel Performance Libraries** (also contained in Intel Parallel Studio). It also uses the **Pardiso** solver, for which the user needs to request a (free academic) licence. 

## Contact
For questions, contact Alexandra Ion at alexandraion@cmu.edu