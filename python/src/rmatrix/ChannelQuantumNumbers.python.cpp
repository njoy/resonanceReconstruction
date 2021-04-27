// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/ChannelQuantumNumbers.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapChannelQuantumNumbers( python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::ChannelQuantumNumbers;
  using OrbitalAngularMomentum = njoy::resonanceReconstruction::rmatrix::OrbitalAngularMomentum;
  using Spin = njoy::resonanceReconstruction::rmatrix::Spin;
  using TotalAngularMomentum = njoy::resonanceReconstruction::rmatrix::TotalAngularMomentum;
  using Parity = njoy::resonanceReconstruction::rmatrix::Parity;

  // create the component
  python::class_< Component > component(

    module,
    "ChannelQuantumNumbers",
    "The (l,s,J,parity) quantum numbers of a reaction channel"
  );

  // wrap the component
  component
  .def(

    python::init< const OrbitalAngularMomentum&, const Spin&,
                  const TotalAngularMomentum&, const Parity& >(),
    python::arg( "l" ), python::arg( "s" ),
    python::arg( "J" ), python::arg( "parity" ),
    "Initialise the quantum numbers\n\n"
    "Arguments:\n"
    "    self      the quantum numbers\n"
    "    l         the orbital angular momentum\n"
    "    s         the channel spin\n"
    "    J         the total angular momentum\n"
    "    parity    the parity"
  )
  .def_property_readonly(

    "orbital_angular_momentum",
    &Component::orbitalAngularMomentum,
    "The orbital angular momentum l of the channel"
  )
  .def_property_readonly(

    "spin",
    &Component::spin,
    "The spin s of the channel"
  )
  .def_property_readonly(

    "total_angular_momentum",
    &Component::totalAngularMomentum,
    "The total angular momentum J of the channel"
  )
  .def_property_readonly(

    "parity",
    &Component::parity,
    "The parity of the channel"
  )
  .def(

    "__repr__",
    [] ( const Component& self ) { return self.toString(); },
    "Convenience function for printing the identifier"
  );
}

} // namespace rmatrix
