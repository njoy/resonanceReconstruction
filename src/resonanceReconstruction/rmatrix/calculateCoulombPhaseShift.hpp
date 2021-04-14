#ifndef NJOY_R2_RMATRIX_CALCULATEW
#define NJOY_R2_RMATRIX_CALCULATEW

// system includes

// other includes
#include "range/v3/numeric/accumulate.hpp"
#include "range/v3/view/indices.hpp"
#include "range/v3/view/transform.hpp"
#include "resonanceReconstruction/rmatrix/options.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @brief Default value for the Coulomb phase shift
 *
 *  For all channel types except for charged particle channels, the Coulomb
 *  phase shift is 0.0.
 */
template < typename ChannelType >
double calculateCoulombPhaseShift( const int, const double ) {
  return 0.0;
}

/**
 *  @brief The Coulomb phase shift for charged particle channels
 *
 *  For a charged particle channel, the Coulomb phase shift W is defined as:
 *     W = 0.0                                             (l == 0)
 *     W = \Sum_{n=1}^{l} \tan^{-1}( \frac{ \eta }{ n } )  (l != 0)
 *  in which l is the orbital angular momentum of the channel and eta is the
 *  parameter defined as follows:
 *     eta = za * zb * mu / hbar / k
 *  in which za and zb are the electrical charge of the particles in the
 *  particle pair, mu is the reduced mass of the particle pair, hbar is the
 *  Planck constant and k is the wave number.
 *
 *  @param[in] l     the oribital angular momentum
 *  @param[in] eta   the value fo the eta parameter to be used
 */
template <>
double calculateCoulombPhaseShift< ChargedParticle >( const int l,
                                                      const double eta ) {
  return ranges::accumulate(
             ranges::views::indices( 1, l + 1 )
               | ranges::views::transform(
                     [&] ( int n ) -> double
                         { return 1. / std::tan( eta / n ); } ), 0.0 );
}

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
