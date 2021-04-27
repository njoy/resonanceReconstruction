// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/Particle.hpp"
#include "conversion.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapParticle( python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::Particle;
  using ParticleID = njoy::elementary::ParticleID;
  using Spin = njoy::resonanceReconstruction::rmatrix::Spin;
  using Parity = njoy::resonanceReconstruction::rmatrix::Parity;

  // create the component
  python::class_< Component > component(

    module,
    "Particle",
    "Particle information\n\n"
    "The Particle class contains specific information for a particle as used\n"
    "during resonance reconstruction. The Particle has an atomic mass, an\n"
    "electrical charge, a spin and a parity (either + or -).\n\n"
    "These variables are used to calculate quantities like the wave number k for\n"
    "an incident reaction channel."
  );

  // wrap the component
  component
  .def(

    python::init( [] ( const ParticleID& id, double mass, double charge,
                       const Spin& spin, const Parity& parity )
                     { return Component( id, toAtomicMass( mass ),
                                         toElectricalCharge( charge ),
                                         spin, parity ); } ),
    python::arg( "id" ), python::arg( "mass" ), python::arg( "charge" ),
    python::arg( "spin" ), python::arg( "parity" ),
    "Initialise the particle\n\n"
    "Arguments:\n"
    "    self      the particle\n"
    "    id        the particle id or name (e.g. n, U235, U235_e1)\n"
    "    mass      the particle or nuclide mass (in amu or dalton)\n"
    "    charge    the charge of the particle or nuclide (in Coulomb)\n"
    "    spin      the particle spin\n"
    "    parity    the particle parity"
  )
  .def(

    python::init( [] ( const ParticleID& id, double mass,
                       const Spin& spin, const Parity& parity )
                     { return Component( id, toAtomicMass( mass ),
                                         spin, parity ); } ),
    python::arg( "id" ), python::arg( "mass" ),
    python::arg( "spin" ), python::arg( "parity" ),
    "Initialise the particle without specifying the electrical charge\n\n"
    "Arguments:\n"
    "    self      the particle\n"
    "    id        the particle id or name (e.g. n, U235, U235_e1)\n"
    "    mass      the particle or nuclide mass (in amu or dalton)\n"
    "    spin      the particle spin\n"
    "    parity    the particle parity"
  )
  .def_property_readonly(

    "particle_id",
    &Component::particleID,
    "The identifier of the particle"
  )
  .def_property_readonly(

    "mass",
    [] ( const Component& self ) -> decltype(auto)
       { return removeUnit( self.mass() ); },
    "The atomic mass of the particle"
  )
  .def_property_readonly(

    "charge",
    [] ( const Component& self ) -> decltype(auto)
       { return removeUnit( self.charge() ); },
    "The electrical charge of the particle"
  )
  .def_property_readonly(

    "spin",
    &Component::spin,
    "The spin of the particle"
  )
  .def_property_readonly(

    "parity",
    &Component::parity,
    "The particle parity"
  );
}

} // namespace rmatrix
