![Continuous Integration](https://github.com/njoy/resonanceReconstruction/workflows/Continuous%20Integration/badge.svg)

# R2: resonance reconstruction

Toolkit for reconstructing cross sections from resonance parameters in the resolved and unresolved resonance region using various R-matrix formalisms. This toolkit provides a full C++ library along with python bindings.

## R2 in python

The python bindings for R2 are still work in progress and should be used accordingly. Please report any issues encountered while using the python bindings using the issue tracker on this repository.

### Installing R2 for python

First of all, a user should clone the resonanceReconstruction repository and build the python bindings:
```
git clone https://github.com/njoy/resonanceReconstruction
cd ENDFtk
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make resonanceReconstruction.python -j8
cp _deps/endftk-build/*.so .
cp _deps/elementary-build/*.so .
```

The R2 python bindings require the ENDFtk and elementary python bindings as well. These are also compiled as part of the `resonanceReconstruction.python` but they are not copied into the build directory (they can be found in the `_deps` subdirectory of the build directory). The last two lines in the commands given above copy the python bindings of ENDFtk and elementary to the top level of the build directory.

The compilation will produce a number of dynamic libraries linked to the python libraries on the user's computer (these will be named something like `< component >.cpython-37m-darwin.so` with `< component >` the name of the component). The names of these dynamic libraries will also indicate which version of the python libraries they are linked against. This is important since you will need to use the associated python version along with them.

R2 in python requires python 3.x so you will need to have at least one python 3.x installed. When multiple python versions are installed, it may be beneficial to include `-DPYTHON_EXECUTABLE=$(which python3)` in the cmake configuration step so that the default python3 version will be picked.

In order to use the R2 python package, the user should make sure that the libraries are within the python path. This can be done in multiple ways. You can set that up by adding the R2 build path to the python path `$PYTHONPATH` environmental variable on your machine, or by using the following in your python code:
```
import sys
sys.path.append( < dynamic-library-location-path > )
```
where `< dynamic-library-location-path >` is the path to the python dynamic libraries.

When running python in the build directory directly, none of these steps are required.

#### Troubleshooting ####

##### CMake doesn’t detect the right Python version #####

Taken from the pybind11 FAQ.

The CMake-based build system will try to automatically detect the installed version of Python and link against that. When this fails, or when there are multiple versions of Python and it finds the wrong one, delete CMakeCache.txt and then add -DPYTHON_EXECUTABLE=$(which python) to your CMake configure line. (Replace $(which python) with a path to python if your prefer.)

A version of python 3.x is preferred.

##### importError cannot import name <sysconfig> #####

This error sometimes comes up when running the cmake command. This appears to be related to an incomplete/corrupted python installation. It can be rectified by installing the distutils package for the python version that is being used. On a linux system, the following command should install the distutils package:
```
sudo apt install python3-distutils
```

##### cannot find python.h #####

When compiling the python bindings, this error indicates that the python header files and static library we need to link to are not installed on your system. This appears to be related to an incomplete python installation. It can be rectified by installing the python3-dev package (when using python 3). On a linux system, the following command should install the header files:
```
sudo apt install python3-dev
```

### A minimal user guide:

#### Building a compound system using Reich-Moore and general R-matrix

R-matrix theory gives us observables for a compound nucleus, for instance cross sections, secondary angular distributions and particle spectra. These are obtained for the compound nucleus, regardless of how the compound system was created. For example: O17 can be created by n+O16, p+N16 or a+C13. As a result, the R-matrix analysis for the O17 compound system can thus give us cross sections for O16(n,a)C13, O16(n,n)O16, O16(n,p)N16, C13(a,n)O16, etc.

The basic building blocks for R-matrix are particle channels defined by:
- A pair of particles making up the channel, e.g. n,O16
- A set of quantum numbers composed of the orbital momentum l, channel spin s and the angular momentum and parity Jpi

## LICENSE
This software is copyrighted by Los Alamos National Laboratory and distributed
according to the conditions in the accompanying [LICENSE](LICENSE) file.
