// system includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local includes
#include "resonanceReconstruction/rmatrix/ResonanceTable.hpp"
#include "conversion.hpp"
#include "views.hpp"

// namespace aliases
namespace python = pybind11;

namespace rmatrix {

void wrapResonanceTable( python::module& module, python::module& viewmodule ) {

  // type aliases
  using Component = njoy::resonanceReconstruction::rmatrix::ResonanceTable;
  using ChannelID = njoy::resonanceReconstruction::rmatrix::ChannelID;
  using Resonance = njoy::resonanceReconstruction::rmatrix::Resonance;
  using ChannelIDRange = RandomAccessAnyView< ChannelID >;
  using ResonanceRange = RandomAccessAnyView< Resonance >;

  // wrap views created by this component
  // none of these are supposed to be created directly by the user
  wrapRandomAccessAnyViewOf< ChannelID >(
      viewmodule,
      "any_view< ChannelID, random_access >" );
  wrapRandomAccessAnyViewOf< Resonance >(
      viewmodule,
      "any_view< Resonance, random_access >" );

  // create the component
  python::class_< Component > component(

    module,
    "ResonanceTable",
    "Resolved resonance parameters for a specific J,pi value"
  );

  // wrap the component
  component
  .def(

    python::init< std::vector< ChannelID >&&, std::vector< Resonance >&& >(),
    python::arg( "channels" ), python::arg( "resonances" ),
    "Initialise the resonance table\n\n"
    "The ResonanceTable class takes the reduced widths for all channels in a\n"
    "spin group. As the channels can differe from case to case, channels are\n"
    "identified using their channel IDs. Data can be extracted from the table\n"
    "using these IDs.\n\n"
    "Arguments:\n"
    "    self         the resonance table\n"
    "    channels     the channel identifiers\n"
    "    resonances   the resolved resonance parameters"
  )
  .def_property_readonly(

    "number_channels",
    &Component::numberChannels,
    "The number of channels"
  )
  .def_property_readonly(

    "number_resonances",
    &Component::numberResonances,
    "The number of resonances"
  )
  .def_property_readonly(

    "channels",
    [] ( const Component& self ) -> ChannelIDRange
       { return self.channels(); },
    "The channel identifiers"
  )
  .def_property_readonly(

    "resonances",
    [] ( const Component& self ) -> ResonanceRange
       { return self.resonances(); },
    "The resonances"
  )
  .def_property_readonly(

    "energies",
    [] ( const Component& self ) -> DoubleRange
       { return removeArrayUnit( self.energies() ); },
    "The resonance energies"
  );
}

} // namespace rmatrix
