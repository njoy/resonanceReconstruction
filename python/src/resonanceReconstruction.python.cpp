// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// internal includes
#include "views.hpp"

// namespace aliases
namespace python = pybind11;

// declarations
namespace rmatrix {

  void wrapChannelQuantumNumbers( python::module&, python::module& );
  void wrapParticle( python::module&, python::module& );
  void wrapParticlePair( python::module&, python::module& );
  void wrapChannelRadii( python::module&, python::module& );
  void wrapParticleChannels( python::module&, python::module& );
  void wrapParticleChannelData( python::module&, python::module& );
  void wrapResonance( python::module&, python::module& );
  void wrapResonanceTable( python::module&, python::module& );
}

/**
 *  @brief ENDFtk python bindings
 *
 *  The name given here (elementary) must be the same as the name
 *  set on the PROPERTIES OUTPUT_NAME in the CMakeLists.txt file.
 */
PYBIND11_MODULE( resonanceReconstruction, module ) {

  // create the views submodule
  python::module viewmodule = module.def_submodule(

    "sequence",
    "sequence - resonance reconstruction sequences (internal use only)"
  );

  // wrap some basic recurring views
  // none of these are supposed to be created directly by the user
  wrapBasicRandomAccessAnyViewOf< double >(
      viewmodule,
      "any_view< double, random_access >" );

  rmatrix::wrapChannelQuantumNumbers( module, viewmodule );
  rmatrix::wrapParticle( module, viewmodule );
  rmatrix::wrapParticlePair( module, viewmodule );
  rmatrix::wrapChannelRadii( module, viewmodule );
  rmatrix::wrapParticleChannels( module, viewmodule );
  rmatrix::wrapParticleChannelData( module, viewmodule );
  rmatrix::wrapResonance( module, viewmodule );
  rmatrix::wrapResonanceTable( module, viewmodule );
}
