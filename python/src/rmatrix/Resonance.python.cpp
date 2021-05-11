// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/Resonance.hpp"
#include "conversion.hpp"
#include "views.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapResonance( python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::Resonance;

  // wrap views created by this component

  // create the component
  python::class_< Component > component(

    module,
    "Resonance",
    "Resolved resonance parameters for a single resonance\n\n"
    "The Resonance class contains the reduced widths (given in sqrt(eV)) for a\n"
    "number of channels. This class is both compatible with the general R matrix\n"
    "theory  and with the Reich-Moore approximation with a single eliminated\n"
    "capture channel. When using the general R matrix theory, the eliminated\n"
    "width is zero."
  );

  // wrap the component
  component
  .def(

    python::init( [] ( double energy, std::vector< double >&& widths,
                       double eliminated )
                     { return Component( toEnergy( energy ),
                                         toReducedWidthArray( widths ),
                                         toReducedWidth( eliminated ) ); } ),
    python::arg( "energy" ), python::arg( "widths" ),
    python::arg( "eliminated" ) = 0.,
    "Initialise the Resonance\n\n"
    "Both energy values and reduced widths may be negative. When the energy is\n"
    "negative, the penetrability, shift factor, phase shift, etc. are calculated\n"
    "at the absolute value of the energy.\n\n"
    "When using this Resonance class with general R matrix theory, the reduced\n"
    "eliminated capture width does not need to be specified (it has a default\n"
    "value of zero).\n\n"
    "Arguments:\n"
    "energy       the resonance energy (in eV, may be negative)\n"
    "widths       the reduced widths (in sqrt(eV))\n"
    "eliminated   the reduced eliminated capture width (in sqrt(eV))"
  )
  .def_property_readonly(

    "energy",
    [] ( const Component& self ) -> decltype(auto)
       { return removeUnit( self.energy() ); },
    "The resonance energy (in eV)"
  )
  .def_property_readonly(

    "eliminated_width",
    [] ( const Component& self ) -> decltype(auto)
       { return removeUnit( self.eliminatedWidth() ); },
    "The eliminated capture width (in sqrt(eV))"
  )
  .def_property_readonly(

    "widths",
    [] ( const Component& self ) -> DoubleRange
       { return removeArrayUnit( self.widths() ); },
    "The reduced widths for this resonance (in sqrt(eV))"
  );
}

} // namespace rmatrix
