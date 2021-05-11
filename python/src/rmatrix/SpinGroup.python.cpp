// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/options.hpp"
#include "resonanceReconstruction/rmatrix/SpinGroup.hpp"
#include "conversion.hpp"
#include "views.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

template < typename Formalism, typename BoundaryOption >
void wrapSpinGroup( const std::string& name, python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::SpinGroup< Formalism, BoundaryOption >;
  using ParticleChannel = njoy::resonanceReconstruction::rmatrix::ParticleChannel;
  using ParticleChannelData = njoy::resonanceReconstruction::rmatrix::ParticleChannelData;
  using ResonanceTable = njoy::resonanceReconstruction::rmatrix::ResonanceTable;

  // wrap views created by this component

  // create the component
  python::class_< Component > component(

    module,
    name.c_str(),
    "A spin group corresponding to a Jpi quantum number set"
  );

  // wrap the component
  component
  .def(

    python::init( [] ( std::vector< ParticleChannel >&& channels,
                       ResonanceTable table )
                     { return Component( std::move( channels ),
                                         std::move( table ) ); } ),
    python::arg( "channels" ), python::arg( "table" ),
    "Initialise the spin group\n\n"
    "Arguments:\n"
    "    self       the spin group\n"
    "    channels   the channels involved in the spin group\n"
    "    table      the table of resonance parameters"
  )
  .def(

    python::init< std::vector< ParticleChannelData >&& >(),
    python::arg( "channels" ),
    "Initialise the spin group\n\n"
    "Arguments:\n"
    "    self       the spin group\n"
    "    channels   the channel data involved in the spin group"
  );
}

void wrapSpinGroups( python::module& module, python::module& viewmodule ) {

  // formalism options
  using ReichMoore = njoy::resonanceReconstruction::rmatrix::ReichMoore;

  // boundary condition options
  using ShiftFactor = njoy::resonanceReconstruction::rmatrix::ShiftFactor;
  using Constant = njoy::resonanceReconstruction::rmatrix::Constant;

  // wrap spin groups
  wrapSpinGroup< ReichMoore, Constant >( "SpinGroup<ReichMoore,Constant>", module, viewmodule );
  wrapSpinGroup< ReichMoore, ShiftFactor >( "SpinGroup<ReichMoore,ShiftFactor>", module, viewmodule );
}

} // namespace rmatrix
