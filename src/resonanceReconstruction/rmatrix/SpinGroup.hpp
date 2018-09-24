/**
 *  @typedef
 *  @brief Channel types
 *
 *  Some of the data for a given channel will actually depend on the particle
 *  type of the channel. The penetrability, shift factor, phase shift and 
 *  coulomb phase shift will depend on whether the channel is a neutron, photon,
 *  charged particle or fission channel. This variant will allow us to capture
 *  that distinction.
 */
using ParticleChannel = std::variant< Channel< Neutron >,
                                      Channel< Photon >,
                                      Channel< ChargedParticle >,
                                      Channel< Fission > >;

/**
 *  @class
 *  @brief A spin group corresponding to a J,pi quantum number set
 */
template < typename Formalism, typename BoundaryOption >
class SpinGroup {

  /* fields */
  std::vector< ReactionID > reactions_;
  Matrix< std::complex< double > > matrix_;
  std::vector< std::complex< double > > diagonalLMatrix_;

  std::vector< unsigned int > incident_;
  std::vector< ParticleChannel > channels_;
  ResonanceTable parameters_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/determineIncidentChannels.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/makeReactionIdentifiers.hpp"

  #include "resonanceReconstruction/rmatrix/SpinGroup/src/penetrabilities.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/shiftFactors.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/phaseShifts.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/coulombShifts.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/boundaryConditions.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/calculateLDiagonal.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/calculateRLMatrix.hpp"

  #include "resonanceReconstruction/rmatrix/SpinGroup/src/verifyChannels.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/verifyIncidentChannels.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/verifyResonanceChannels.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/ctor.hpp"

  /**
   *  @brief Return the channels in the spin group
   */
  auto channels() const { return ranges::view::all( this->channels_ ); }

  /**
   *  @brief Return the current incident channels in the spin group
   */
  auto incidentChannels() const {
    return ranges::view::all( this->incident_ )
             | ranges::view::transform( [&] ( const unsigned int i )
                                            { return this->channels_[i]; } ); }

  /**
   *  @brief Return the current incident particle pair
   */
  const ParticlePair incidentPair() const {

    //! @todo for some reason debug fails when using const ParticlePair& as
    //!       the return type for this function and the lambda function inside
    //!       the visit function
    return std::visit(
               [&] ( const auto& channel ) -> const ParticlePair
                   { return channel.particlePair(); },
               this->incidentChannels().front() );
  };

  /**
   *  @brief Return the reaction identifiers associated to each channel
   *
   *  These reaction identifiers are not unique: a given reaction identifier
   *  may appear multiple times because multiple channels can contribute to the
   *  same reactions (e.g. multiple fission channels). The order in which these 
   *  are given equals the order of the channels in the spin group.
   */
  auto reactions() const { return ranges::view::all( this->reactions_ ); }

  /**
   *  @brief Return the resonance table
   */
  const ResonanceTable& resonanceTable() const { return this->parameters_; }

  #include "resonanceReconstruction/rmatrix/SpinGroup/src/switchIncidentPair.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/evaluate.hpp"
};

