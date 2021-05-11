// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/ParticleChannelData.hpp"
#include "conversion.hpp"
#include "views.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapParticleChannelData( python::module& module, python::module& ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::ParticleChannelData;
  using ParticleChannel = njoy::resonanceReconstruction::rmatrix::ParticleChannel;

  // wrap views created by this component

  // create the component
  python::class_< Component > component(

    module,
    "ParticleChannelData",
    "Data for a single particle channel\n\n"
    "This class represent all data related to a single channel; the particle\n"
    "channel instance, the resonance energies and reduced widths. It provides\n"
    "an alternate method of constructing the SpinGroup and CompoundSystem for a\n"
    "set of channels and is most useful when the data is not presented by spin\n"
    "group to begin with, or if the spingroups are not unique."
  );

  // wrap the component
  component
  .def(

    python::init( [] ( const ParticleChannel& channel,
                       std::vector< double >&& energies,
                       std::vector< double >&& widths,
                       bool eliminated )
                     { return Component( channel, toEnergyArray( energies ),
                                         toReducedWidthArray( widths ),
                                         eliminated ); } ),
    python::arg( "channel" ), python::arg( "energies" ), python::arg( "widths" ),
    python::arg( "eliminated" ) = false,
    "Initialise the particle channel data\n\n"
    "Arguments:\n"
    "    self       the particle channel data\n"
    "    channel    the particle channel\n"
    "    energies   the resonanc energies\n"
    "    widths     the reduced resonance widths"
  )
  .def_property_readonly(

    "channel",
    &Component::channel,
    "The channel"
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

    "energies",
    [] ( const Component& self ) -> DoubleRange
       { return removeArrayUnit( self.energies() ); },
    "The energy values"
  )
  .def_property_readonly(

    "widths",
    [] ( const Component& self ) -> DoubleRange
       { return removeArrayUnit( self.widths() ); },
    "The reduced resonance widths"
  );
}

} // namespace rmatrix
