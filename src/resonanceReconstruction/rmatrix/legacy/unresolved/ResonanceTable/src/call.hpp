/**
 *  @brief Return the unresolved resonance parameters at the given energy
 *
 *  @param[in] energy       the incident energy
 */
auto operator()( const Energy& energy ) const {

  return Resonance( energy,
                    this->level_spacing_table_( energy ),
                    this->elastic_table_( energy ),
                    this->capture_table_( energy ),
                    this->fission_table_( energy ),
                    this->competition_table_( energy ) );
}
