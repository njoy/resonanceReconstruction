/**
 *  @class
 *  @brief Resolved resonance data for a specific l,J spin group using SLBW
 *
 *  This class contains the resolved resonance parameters and the associated
 *  incident channel for legacy SLBW ENDF data.
 */
template <>
class SpinGroup< SingleLevelBreitWigner > : protected SpinGroupBase< resolved::ResonanceTable > {

  /* aliases */
  using EnergySquared = decltype( std::declval< Energy >() *
                                  std::declval< Energy >() );

  /* fields */
  Energy qx_;

  std::vector< Width > elastic_;
  std::vector< Width > capture_;
  std::vector< Width > fission_;
  std::vector< Width > total_;
  std::vector< Energy > delta_;
  std::vector< EnergySquared > denominator_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/src/precompute.hpp"

protected:

  using SpinGroupBase::elasticID;
  using SpinGroupBase::captureID;
  using SpinGroupBase::fissionID;

  auto elastic() const { return ranges::view::all( this->elastic_ ); }
  auto capture() const { return ranges::view::all( this->capture_ ); }
  auto fission() const { return ranges::view::all( this->fission_ ); }
  auto total() const { return ranges::view::all( this->total_ ); }
  auto delta() const { return ranges::view::all( this->delta_ ); }
  auto denominator() const { return ranges::view::all( this->denominator_ ); }

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/src/ctor.hpp"

  /* methods */

  using SpinGroupBase::incidentChannel;
  using SpinGroupBase::orbitalAngularMomentum;
  using SpinGroupBase::totalAngularMomentum;
  using SpinGroupBase::resonanceTable;

  /**
   *  @brief Return competitive Q value
   */
  const Energy& QX() const { return this->qx_; }

  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/src/grid.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/src/evaluate.hpp"
};
