#ifndef NJOY_R2_RMATRIX_CALCULATEPHI
#define NJOY_R2_RMATRIX_CALCULATEPHI

// system includes
#include <complex>

// other includes
#include "utility/horner.hpp"
#include "resonanceReconstruction/rmatrix/src/coh3-coulomb.hpp"
#include "resonanceReconstruction/rmatrix/ChannelTypes.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @brief Default value for the phase shift
 *
 *  For all channel types except for neutron and charged particle channels, the
 *  shift factor is 0.0.
 */
template < typename Type >
double calculatePhaseShift( const unsigned int, const double, const double ) {

  return 0.0;
}

/**
 *  @brief The phase shift for neutron channels
 *
 *  For a neutron channel, the phase shift is defined using the neutron
 *  hardsphere functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 */
template <>
double calculatePhaseShift< Neutron >( const unsigned int l,
                                       const double ratio,
                                       const double ) {

  double squared = ratio * ratio;
  switch ( l ) {

    case 0 : return ratio;
    case 1 : return ratio - std::atan( ratio );
    case 2 : {

      double offset = 3. * ratio / ( 3. - squared );
      return ratio - std::atan( offset );
    }
    case 3 : {

      double offset = ratio * ( 15. - squared ) / ( 15. - 6. * squared );
      return ratio - std::atan( offset );
    }
    case 4 : {

      constexpr std::array< double, 3 > denominator = {{ 105., -45., 1. }};
      double offset = ratio * ( 105. - 10. * squared )
                      / utility::horner( std::rbegin( denominator ),
                                         std::rend( denominator ),
                                         squared );
      return ratio - std::atan( offset );
    }
    default : throw std::exception();
  }
}

/**
 *  @brief The phase shift for charged particle channels
 *
 *  For a charged particle channel, the phase shift is defined using the
 *  Coulomb wave functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 *  @param[in] eta     the Sommerfeld parameter
 */
template <>
double calculatePhaseShift< ChargedParticle >( const unsigned int l,
                                               const double ratio,
                                               const double eta ) {

  std::complex< double > gf, dgf;
  coulombWaveFunctions( l, ratio, eta, gf, dgf);

  double F = gf.imag();
  double G = gf.real();
  return ( F == 0. ) and ( G == 0. )
           ? 0.
           : std::acos( G / std::sqrt( F * F + G * G ) );
}

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
