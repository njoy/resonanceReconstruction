#ifndef NJOY_R2_RMATRIX_CALCULATEP
#define NJOY_R2_RMATRIX_CALCULATEP

// system includes
#include <complex>

// other includes
#include "utility/horner.hpp"
#include "resonanceReconstruction/rmatrix/src/coh3-coulomb.hpp"
#include "resonanceReconstruction/rmatrix/ChannelType.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @brief Default value for the penetrability
 *
 *  For all channel types except for neutron and charged particle channels, the
 *  penetrability is 1.0.
 */
template < typename Type >
double calculatePenetrability( const unsigned int, const double, const double ) {

  return 1.0;
}

/**
 *  @brief The penetrability for neutron channels
 *
 *  For a neutron channel, the penetrability is defined using the neutron
 *  hardsphere functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 */
template <>
double calculatePenetrability< Neutron >( const unsigned int l,
                                          const double ratio,
                                          const double  ) {

  double squared = ratio * ratio;
  switch ( l ) {

    case 0 : return ratio;
    case 1 : return ratio * squared / ( 1. + squared );
    case 2 : {

      constexpr std::array< double, 3 > denominator = {{ 9., 3., 1. }};
      return ratio * squared * squared
             / utility::horner( std::rbegin( denominator ),
                                std::rend( denominator ),
                                squared );
    }
    case 3 : {

      constexpr std::array< double, 4 > denominator = {{ 225., 45., 6., 1. }};
      return ratio * squared * squared * squared
             / utility::horner( std::rbegin( denominator ),
                                std::rend( denominator ),
                                squared );
    }
    case 4 : {

      constexpr std::array< double, 5 > denominator = {{ 11025., 1575., 135., 10., 1. }};
      return ratio * squared * squared * squared * squared
             / utility::horner( std::rbegin( denominator ),
                                std::rend( denominator ),
                                squared );
    }
    default : throw std::exception();
  }
}

/**
 *  @brief The penetrability for charged particle channels
 *
 *  For a charged particle channel, the penetrability is defined using the
 *  Coulomb wave functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 *  @param[in] eta     the Sommerfeld parameter
 */
template <>
double calculatePenetrability< ChargedParticle >( const unsigned int l,
                                                  const double ratio,
                                                  const double eta ) {


  std::complex< double > gf, dgf;
  coulombWaveFunctions( l, ratio, eta, gf, dgf);

  double F = gf.imag();
  double G = gf.real();
  return ( F == 0. ) and ( G == 0. ) ? 0. : ratio / ( F * F + G * G );
}

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
