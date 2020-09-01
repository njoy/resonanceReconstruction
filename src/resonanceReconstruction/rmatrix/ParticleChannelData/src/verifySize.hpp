static
void verifySize( const std::vector< Energy >& energies,
                 const std::vector< ReducedWidth >& widths ) {

  // verify that both are of the same size (no elements is allowed)

  if ( energies.size() != widths.size() ) {

    Log::error( "The number of energies and reduced resonance widths is not the"
                " same." );
    throw std::exception();
  }
}
