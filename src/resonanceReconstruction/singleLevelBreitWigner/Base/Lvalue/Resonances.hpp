/* Resonances is an alternative implementation of the storage and computation
 * of the resonance information designed to accomodate vectorized processing.
 * A scalar implementation was ultimately chosen over this vectorized 
 * implementation for reasons of expediency. 
 *
 * This implementation has been left as a starting point for potential future 
 * performance research.
 */
struct Resonances {
  std::vector< double > data;
  int rowLength;
  
  /* Rather than being stored contiguously, data subsets (energy, neutronWidths, 
   * etc) are organized into blocks along block boundaries, defined to the size
   * of a cacheline on commodity hardware at the time of writing, potentially
   * leaving gaps in the vector. 
   *
   * This extra step ensures data subsets are share an alignment. The block size
   * computation is done explicitly to provide a hint to the compiler regarding 
   * this property.
   */
  static constexpr int blockSize = 8;
  int rowBlocks;
  
  auto energies() const {
    return this->data
           | ranges::view::take_exactly( this->rowLength )
           | ranges::view::transform([](double d){ return d * electronVolts; });
  }
  
  auto neutronWidths() const {
    return this->data
           | ranges::view::drop_exactly( this->rowBlocks * blockSize )
           | ranges::view::take_exactly( this->rowLength )
           | ranges::view::transform([](double d){ return d * electronVolts; });

  }
  
  auto captureWidths() const {
    return this->data
           | ranges::view::drop_exactly( 2 * this->rowBlocks * blockSize )
           | ranges::view::take_exactly( this->rowLength )
           | ranges::view::transform([](double d){ return d * electronVolts; });

  }
  
  auto fissionWidths() const {
    return this->data
           | ranges::view::drop_exactly( 3 * this->rowBlocks * blockSize )
           | ranges::view::take_exactly( this->rowLength )
           | ranges::view::transform([](double d){ return d * electronVolts; });
  }
  
  auto weightedCompetitiveWidths() const {
    return this->data
           | ranges::view::drop_exactly( 4 * this->rowBlocks * blockSize )
           | ranges::view::take_exactly( this->rowLength )
           | ranges::view::transform([](double d){ return d * electronVolts; });
  }
  
  auto inversePenetrationFactors() const {
    return this->data
           | ranges::view::drop_exactly( 5 * this->rowBlocks * blockSize )
           | ranges::view::take_exactly( this->rowLength );
  }
  
  auto shiftFactors() const {
    return this->data
           | ranges::view::drop_exactly( 6 * this->rowBlocks * blockSize )
           | ranges::view::take_exactly( this->rowLength );
  }
  
  auto statisticalFactors() const {
    return this->data
           | ranges::view::drop_exactly( 7 * this->rowBlocks * blockSize )
           | ranges::view::take_exactly( this->rowLength );
  }

  #include "resonanceReconstruction/singleLevelBreitWigner/Base/Lvalue/Resonances/src/call.hpp"
};
