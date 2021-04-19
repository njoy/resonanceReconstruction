// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// namespace aliases
namespace python = pybind11;

// declarations
namespace rmatrix {

  void wrapChannelQuantumNumbers( python::module& );
  void wrapParticle( python::module& );
  void wrapParticlePair( python::module& );
  void wrapChannelRadii( python::module& );
}

/**
 *  @brief ENDFtk python bindings
 *
 *  The name given here (elementary) must be the same as the name
 *  set on the PROPERTIES OUTPUT_NAME in the CMakeLists.txt file.
 */
PYBIND11_MODULE( resonanceReconstruction, module ) {

  rmatrix::wrapChannelQuantumNumbers( module );
  rmatrix::wrapParticle( module );
  rmatrix::wrapParticlePair( module );
  rmatrix::wrapChannelRadii( module );
}
