// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/options.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapBoundary( python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::Boundary;

  // wrap views created by this component

  // create the component
  python::enum_< Component > component(

    module,
    "Boundary",
    "The R-matrix boundary condition option enumerator",
    python::arithmetic()
  );

  // wrap the component
  component
  .value( "ShiftFactor", Component::ShiftFactor )
  .value( "Constant", Component::Constant );
}

} // namespace rmatrix
