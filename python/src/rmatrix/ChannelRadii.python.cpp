// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/ChannelRadii.hpp"
#include "conversion.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapChannelRadii( python::module& module ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::ChannelRadii;
  using ChannelRadiusTable = njoy::resonanceReconstruction::rmatrix::ChannelRadiusTable;
  using Variant = Component::ChannelRadiusVariant;

  //! @todo the use of a channel radius table is currently not allowed

  // create the component
  python::class_< Component > component(

    module,
    "ChannelRadii",
    "Channel radii used in wave function calculations\n\n"
    "The penetrability P, shift factor S and phase shift phi require knowledge\n"
    "of the channel radius in their calculation. The ChannelRadii class provides\n"
    "these radii for each on of these."
  );

  // wrap the component
  component
  .def(

    python::init( [] ( double radius )
                     { return Component(
                                Variant( toChannelRadius( radius ) ) ); } ),
    python::arg( "radius" ),
    "Initialise the channel radii using a single radius\n\n"
    "Arguments:\n"
    "    self      the channel radii\n"
    "    radius    the channel radius to be used for P, S and phi"
  )
  .def(

    python::init( [] ( double trueRadius, double effectiveRadius )
                     { return Component(
                                Variant( toChannelRadius( trueRadius ) ),
                                Variant( toChannelRadius( effectiveRadius ) ) ); } ),
    python::arg( "true_radius" ), python::arg( "effective_radius" ),
    "Initialise the channel radii using a true and effective radius\n\n"
    "Arguments:\n"
    "    self                the channel radii\n"
    "    true_radius         the channel radius to be used for P and S\n"
    "    effective_radius    the channel radius to be used for phi"
  )
  .def(

    python::init( [] ( double penetrability, double shiftFactor,
                       double phaseShift )
                     { return Component(
                                Variant( toChannelRadius( penetrability ) ),
                                Variant( toChannelRadius( shiftFactor ) ),
                                Variant( toChannelRadius( phaseShift ) ) ); } ),
    python::arg( "penetrability_radius" ), python::arg( "shift_factor_radius" ),
    python::arg( "phase_shift_radius" ),
    "Initialise the channel radii using three different radii\n\n"
    "Arguments:\n"
    "    self                    the channel radii\n"
    "    penetrability_radius    the channel radius to be used for P\n"
    "    shift_factor_radius     the channel radius to be used for S\n"
    "    phase_shift_radius      the channel radius to be used for phi"
  )
  .def(

    "penetrability_radius",
    [] ( const Component& self, double energy )
       { return removeUnit( self.penetrabilityRadius( toEnergy( energy ) ) ); },
    python::arg( "energy" ),
    "Return the channel radius for the penetrability P\n\n"
    "Arguments:\n"
    "    self     the channel radii\n"
    "    energy   the energy for which the radius must be given"
  )
  .def(

    "shift_factor_radius",
    [] ( const Component& self, double energy )
       { return removeUnit( self.shiftFactorRadius( toEnergy( energy ) ) ); },
    python::arg( "energy" ),
    "Return the channel radius for the shift factor S\n\n"
    "Arguments:\n"
    "    self     the channel radii\n"
    "    energy   the energy for which the radius must be given"
  )
  .def(

    "phase_shift_radius",
    [] ( const Component& self, double energy )
       { return removeUnit( self.phaseShiftRadius( toEnergy( energy ) ) ); },
    python::arg( "energy" ),
    "Return the channel radius for the phase shift phi\n\n"
    "Arguments:\n"
    "    self     the channel radii\n"
    "    energy   the energy for which the radius must be given"
  );
}

} // namespace rmatrix
