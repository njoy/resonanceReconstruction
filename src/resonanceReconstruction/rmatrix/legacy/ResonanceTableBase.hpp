/**
 *  @class
 *  @brief Base interface for a table of resonances for a specific l,J value
 */
template < typename ResonanceType > class ResonanceTableBase {

  /* fields */
  std::vector< ResonanceType > resonances_;

public:

  /* constructor */

  /**
   *  @brief Constructor
   *
   *  @param[in] resonances   the resonances (ne values)
   */
  ResonanceTableBase( std::vector< ResonanceType >&& resonances ) :
      resonances_( resonances ) {}

  /* methods */

  /**
   *  @brief Return the number of resonances
   */
  unsigned int numberResonances() const { return this->resonances_.size(); }

  /**
   *  @brief Return the resonances
   */
  const std::vector< ResonanceType >& resonances() const {

    return this->resonances_;
  }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {

    return this->resonances()
             | ranges::view::transform( [] ( const auto& resonance )
                                           { return resonance.energy(); } );
  }
};
