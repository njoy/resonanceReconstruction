// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// namespace aliases
namespace python = pybind11;

// declarations
// void wrapFromENDF( python::module& );

int add( int i, int j ) {

    return i + j;
}

/**
 *  @brief ENDFtk python bindings
 *
 *  The name given here (elementary) must be the same as the name
 *  set on the PROPERTIES OUTPUT_NAME in the CMakeLists.txt file.
 */
PYBIND11_MODULE( resonanceReconstruction, module ) {

  m.doc() = "pybind11 example plugin"; // optional module docstring

  m.def("add", &add, "A function which adds two numbers");

  // wrapFromENDF( module );
}
