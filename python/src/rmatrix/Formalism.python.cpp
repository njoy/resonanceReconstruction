// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/options.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapFormalism( python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::Formalism;

  // wrap views created by this component

  // create the component
  python::enum_< Component > component(

    module,
    "Formalism",
    "The R-matrix formalism option enumerator",
    python::arithmetic()
  );

  // wrap the component
  component
  .value( "SingleLevelBreitWigner", Component::SingleLevelBreitWigner )
  .value( "MultiLevelBreitWigner", Component::MultiLevelBreitWigner )
  .value( "ReichMoore", Component::ReichMoore )
  .value( "GeneralRMatrix", Component::GeneralRMatrix );
}

} // namespace rmatrix
