// a basic interpolation table templated on x and y types and interpolation type
template < typename XType, typename YType, typename Interpolation >
using Table = interpolation::Table<
                interpolation::table::Type<
                  Interpolation,
                  interpolation::table::search::Binary,
                  interpolation::table::discontinuity::TakeLeft,
                  std::vector< XType >, std::vector< YType > >,
                interpolation::table::right::interval::Throws,
                interpolation::table::left::interval::Throws >;

// histogram interpolation table templated on x and y type
template < typename XType, typename YType >
using HistogramTable = Table< XType, YType, interpolation::Histogram >;

// lin-lin interpolation table templated on x and y type
template < typename XType, typename YType >
using LinLinTable = Table< XType, YType, interpolation::LinearLinear >;

// linear-log interpolation table templated on x and y type
template < typename XType, typename YType >
using LinLogTable = Table< XType, YType, interpolation::LinearLogarithmic >;

// log-linear interpolation table templated on x and y type
template < typename XType, typename YType >
using LogLinTable = Table< XType, YType, interpolation::LogarithmicLinear >;

// log-log interpolation table templated on x and y type
template < typename XType, typename YType >
using LogLogTable = Table< XType, YType, interpolation::LogarithmicLogarithmic >;

// variant composing the different interpolation types
template < typename XType, typename YType >
using Tablevariant =
  interpolation::Table<
    interpolation::table::Variant< HistogramTable< XType, YType >,
                                   LinLinTable< XType, YType >,
                                   LinLogTable< XType, YType >,
                                   LogLinTable< XType, YType >,
                                   LogLogTable< XType, YType > > >;

// a vector of tables
template < typename XType, typename YType >
using TableVector = interpolation::table::Vector< Tablevariant< XType, YType > >;

// interpolation table with multiple interpolation regions
template < typename XType, typename YType >
using MultiRegionTable = interpolation::Table< TableVector< XType, YType > >;
