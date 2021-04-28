// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/ParticleChannel.hpp"
#include "conversion.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

template < typename ChannelType >
void wrapChannel( const std::string& channel, python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::Channel< ChannelType >;
  using ParticlePair = njoy::resonanceReconstruction::rmatrix::ParticlePair;
  using ChannelRadii = njoy::resonanceReconstruction::rmatrix::ChannelRadii;
  using ChannelQuantumNumbers = njoy::resonanceReconstruction::rmatrix::ChannelQuantumNumbers;
  using BoundaryCondition = njoy::resonanceReconstruction::rmatrix::BoundaryCondition;

  // create the component
  python::class_< Component > component(

    module,
    channel.c_str(),
    "A reaction channel\n\n"
    "The Channel class is used to used to store or calculate all channel specific\n"
    "data such as the incident particle pair, the outgoing particle pair, the\n"
    "associated Q-value, etc.\n\n"
    "The incident and outgoing particle pair identifiers determine the reaction\n"
    "to which the channel contributes. For example: 'n,Cl35->n,Cl35' for elastic\n"
    "scattering.\n\n"
    "The channel should always have a unique channel identifier associated to\n"
    "it. In most cases, such a unique identifier can be generated from the\n"
    "channel's particle pair identifier and (l,s,Jpi) set of quantum numbers,\n"
    "for example: 'n,Cl35{0,1,1+}'.\n\n"
    "Under special cases such as partial fission channels, a user defined\n"
    "identifier may be needed since the associated particle pair is the\n"
    "same (since they contribute to the same reaction), for example: \n"
    "'fission1{0,1,1+}' and 'fission2{0,1,1+}''"
  );

  // wrap the component
  component
  .def(

    python::init( [] ( const ParticlePair& incident, const ParticlePair& pair,
                       double q_value, const ChannelQuantumNumbers& numbers,
                       const ChannelRadii& radii,
                       const BoundaryCondition& boundary )
                     { return Component( incident, pair, toQValue( q_value ),
                                         numbers, radii, boundary ); } ),
    python::arg( "incident" ), python::arg( "pair" ), python::arg( "q_value" ),
    python::arg( "numbers" ), python::arg( "radii" ),
    python::arg( "boundary" ) = 0.0,
    "Initialise the channel\n\n"
    "Arguments:\n"
    "    self       the channel\n"
    "    incident   the current incident particle pair\n"
    "    pair       the particle pair associated to the channel\n"
    "    qValue     the Q value associated with the transition from\n"
    "               incident to the outgoing particle pair\n"
    "    numbers    the channel's quantum numbers\n"
    "    radii      the channel radii to be used for calculation of the\n"
    "               wave functions\n"
    "    boundary   the boundary condition (often the value of the\n"
    "               orbital angular momentum of the channel)"
  )
  .def(

    python::init( [] ( const ParticlePair& incident, const ParticlePair& pair,
                       const std::string& id, double q_value,
                       const ChannelQuantumNumbers& numbers,
                       const ChannelRadii& radii,
                       const BoundaryCondition& boundary )
                     { return Component( incident, pair, id, toQValue( q_value ),
                                         numbers, radii, boundary ); } ),
    python::arg( "incident" ), python::arg( "pair" ),
    python::arg( "id" ), python::arg( "qValue" ),
    python::arg( "numbers" ), python::arg( "radii" ),
    python::arg( "boundary" ) = 0.0,
    "Initialise the channel with a custom channel identifier\n\n"
    "In this case, the channel identifier is generated based on the channel's\n"
    "particle pair identifier and the channel's (l,s,Jpi) set of quantum numbers.\n"
    "For example: \'n,Cl35{0,1,1+}\'\n\n"
    "Arguments:\n"
    "    self       the channel\n"
    "    incident   the current incident particle pair\n"
    "    pair       the particle pair associated to the channel\n"
    "    id         the replacement particle pair identifier\n"
    "    qValue     the Q value associated with the transition from\n"
    "               incident to the outgoing particle pair\n"
    "    numbers    the channel's quantum numbers\n"
    "    radii      the channel radii to be used for calculation of the\n"
    "               wave functions\n"
    "    boundary   the boundary condition (often the value of the\n"
    "               orbital angular momentum of the channel)"
  )
  .def_property_readonly(

    "channel_id",
    &Component::channelID,
    "The unique channel identifier"
  )
  .def_property_readonly(

    "reaction_id",
    &Component::reactionID,
    "The identifier for the reaction to which this channel contributes"
  )
  .def_property_readonly(

    "incident_particle_pair",
    &Component::incidentParticlePair,
    "The current incident particle pair in the channel"
  )
  .def_property_readonly(

    "particle_pair",
    &Component::particlePair,
    "The particle pair in the channel"
  )
  .def_property_readonly(

    "quantum_numbers",
    &Component::quantumNumbers,
    "The (l,s,Jpi) quantum number set of the channel"
  )
  .def_property_readonly(

    "radii",
    &Component::radii,
    "The channel radii"
  )
  .def_property_readonly(

    "Q",
    [] ( const Component& self )
       { return removeUnit( self.Q() ); },
    "The Q value for going from the current incident channel to this channel"
  )
  .def_property_readonly(

    "boundary_condition",
    &Component::boundaryCondition,
    "The channel boundary condition"
  )
  .def_property_readonly(

    "is_incident_channel",
    &Component::isIncidentChannel,
    "Return whether or not the channel is an incident channel"
  )
  .def_property_readonly(

    "statistical_spin_factor",
    &Component::statisticalSpinFactor,
    "Return the statistical spin factor\n\n"
    "The statistical spin factor g of a channel is defined as follows:\n"
    "   g = ( 2 * J + 1 ) / ( 2 * ia + 1 ) / ( 2 * ib + 1 )\n"
    "in which J is the total angular momentum of the channel and ia and ib\n"
    "are the spins of the particles in the particle pair.\n\n"
    "For a particle pair involving a neutron (for which the particle spin is\n"
    "0.5), this reduces to:\n"
    "  g = ( 2 * J + 1 ) / ( 2 * I + 1 ) / 2\n"
    "in which I is the target nucleus spin value."
  )
  .def(

    "below_threshold",
    [] ( const Component& self, double energy )
       { return self.belowThreshold( toEnergy( energy ) ); },
    python::arg( "energy" ),
    "Return whether or not the energy is below the threshold for this channel\n\n"
    "The incident energy is below the threshold energy for the channel if\n"
    "    energy * ratio + q < 0.0\n"
    "where energy is the incident energy, ratio is the mass ratio M / ( m + M )\n"
    "for the incident particle pair and q is the Q value for this channel.\n\n"
    "Arguments:\n"
    "    self      the channel\n"
    "    energy    the energy to be tested\n"
  )
  .def(

    "sommerfeld_parameter",
    [] ( const Component& self, double energy )
       { return self.sommerfeldParameter( toEnergy( energy ) ); },
    python::arg( "energy" ),
    "Return the value of Sommerfeld parameter for the channel at a given energy\n\n"
    "The Sommerfeld parameter eta is an energy dependent quantity defined as\n"
    "follows:\n"
    "   eta = z * Z * mu / ( 4 * pi * epsilon0 * hbar^2 * k )\n"
    "in which z and Z are the electrical charge of the particles in the\n"
    "particle pair, mu is the reduced mass of the particle pair, hbar is the\n"
    "Planck constant, k is the wave number and epsilon0 is the vacuum\n"
    "permittivity.\n\n"
    "It is a dimensionless parameter.\n\n"
    "@param energy   the energy at which the sommerfeld parameter needs to be \n"
    "                evaluated\n"
  )
  .def(

    "wave_number",
    [] ( const Component& self, double energy )
       { return removeUnit( self.waveNumber( toEnergy( energy ) ) ); },
    python::arg( "energy" ),
    "Return the wave number of the channel at a given energy\n\n"
    "The wave number k is an energy dependent quantity defined as follows:\n"
    "   hbar^2 k^2 = 2 * mu * ( energy * ratio + q )\n"
    "in which mu is the reduced mass of the channel's particle pair and ratio\n"
    "is the mass ratio M / ( m + M ) for the incident particle pair, q is the\n"
    "Q value associated to the transition of the incident particle pair to the\n"
    "channel's particle pair and hbar is the Planck constant.\n\n"
    "@param energy   the energy at which the wave number needs to be \n"
    "                evaluated\n"
  )
  .def(

    "penetrability",
    [] ( const Component& self, double energy )
       { return self.penetrability( toEnergy( energy ) ); },
    python::arg( "energy" ),
    "Return the penetrability for this channel as a function of energy\n\n"
    "@param energy   the energy at which the penetrability needs to be \n"
    "                evaluated\n"
  )
  .def(

    "shift_factor",
    [] ( const Component& self, double energy )
       { return self.shiftFactor( toEnergy( energy ) ); },
    python::arg( "energy" ),
    "Return the shift factor for this channel as a function of energy\n\n"
    "@param energy   the energy at which the shift factor needs to be \n"
    "                evaluated\n"
  )
  .def(

    "phase_shift",
    [] ( const Component& self, double energy )
       { return self.phaseShift( toEnergy( energy ) ); },
    python::arg( "energy" ),
    "Return the phase shift for this channel as a function of energy\n\n"
    "@param energy   the energy at which the phase shift needs to be \n"
    "                evaluated\n"
  )
  .def(

    "coulomb_phase_shift",
    [] ( const Component& self, double energy )
       { return self.coulombPhaseShift( toEnergy( energy ) ); },
    python::arg( "energy" ),
    "Return the coulomb phase shift for this channel as a function of energy\n\n"
    "@param energy   the energy at which the coulomb phase shift needs to be \n"
    "                evaluated\n"
  );
}

void wrapParticleChannels( python::module& module, python::module& viewmodule ) {

  // channel types
  using Neutron = njoy::resonanceReconstruction::rmatrix::Neutron;
  using Photon = njoy::resonanceReconstruction::rmatrix::Photon;
  using ChargedParticle = njoy::resonanceReconstruction::rmatrix::ChargedParticle;
  using Fission = njoy::resonanceReconstruction::rmatrix::Fission;

  // wrap channels
  wrapChannel< Neutron >( "NeutronChannel", module, viewmodule );
  wrapChannel< Photon >( "PhotonChannel", module, viewmodule );
  wrapChannel< ChargedParticle >( "ChargedParticleChannel", module, viewmodule );
  wrapChannel< Fission >( "FissionChannel", module, viewmodule );
}

} // namespace rmatrix
