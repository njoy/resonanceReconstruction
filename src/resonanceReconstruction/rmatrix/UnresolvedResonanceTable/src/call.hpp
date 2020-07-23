/**
 *  @brief Return the unresolved resonance parameters at the given energy
 *
 *  @param[in] energy       the incident energy
 */
auto operator()( const Energy& energy ) const {

  Energy spacing = this->level_spacing_table_( energy );
  auto widths = this->width_tables_
                  | ranges::view::transform(
                        [energy] ( const auto& table )
                                 { return table( energy ); } );

  return UnresolvedResonance( energy, spacing, widths );
}
