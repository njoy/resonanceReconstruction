// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/ParticlePair.hpp"
#include "conversion.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapParticlePair( python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::ParticlePair;
  using Particle = njoy::resonanceReconstruction::rmatrix::Particle;
  using ParticlePairID = njoy::elementary::ParticlePairID;

  // create the component
  python::class_< Component > component(

    module,
    "ParticlePair",
    "Particle pair information\n\n"
    "A ParticlePair represents the two particles involved in a entrance or exit\n"
    "reaction channel (we assume that the reaction is a two-body reaction). The\n"
    "pair consists of a \"small\" incident or outgoing particle (e.g. a neutron,\n"
    "photon, alpha, etc.) and a \"larger\" target or residual nucleus (e.g. H1,\n"
    "He4, U235, etc.).\n\n"
    "The ParticlePair class gives us access to information related to the\n"
    "pair of particles such as the mass ratio and the reduced mass."
  );

  // wrap the component
  component
  .def(

    python::init< const Particle&, const Particle&, const ParticlePairID& >(),
    python::arg( "particle" ), python::arg( "residual" ), python::arg( "id" ),
    "Initialise the particle pair with a custom identitifer\n\n"
    "In some cases the identifier of the particle pair that is generated\n"
    "automatically using the particle identifiers makes no sense, as would be\n"
    "the case for fission. This constructor can be used to override the\n"
    "automatically generated identifier.\n\n"
    "Arguments:\n"
    "    self        the particle pair\n"
    "    particle    the incident or outgoing particle\n"
    "    residual    the target or residual nuclide\n"
    "    id          a unique identifier for the particle pair"
  )
  .def(

    python::init< const Particle&, const Particle& >(),
    python::arg( "particle" ), python::arg( "residual" ),
    "Initialise the particle pair\n\n"
    "The identifier of the particle pair is generated automatically using the\n"
    "particle identifiers (e.g. n,U235_e0) when calling this constructor.\n\n"
    "Arguments:\n"
    "    self        the particle pair\n"
    "    particle    the incident or outgoing particle\n"
    "    residual    the target or residual nuclide"
  )
  .def_property_readonly(

    "particle",
    &Component::particle,
    "The first particle of the particle pair"
  )
  .def_property_readonly(

    "residual",
    &Component::residual,
    "The second particle of the particle pair"
  )
  .def_property_readonly(

    "pair_id",
    &Component::pairID,
    "The identifier of the particle pair"
  )
  .def_property_readonly(

    "mass_ratio",
    &Component::massRatio,
    "The mass ratio mb / ( ma + mb ) of the particle pair"
  )
  .def_property_readonly(

    "reduced_mass",
    [] ( const Component& self ) -> decltype(auto)
       { return removeUnit( self.reducedMass() ); },
    "The reduced mass mu = ma * mb / ( ma + mb ) of the particle pair"
  );
}

} // namespace rmatrix
