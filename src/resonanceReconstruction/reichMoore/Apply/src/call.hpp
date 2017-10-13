template< typename Functor >
decltype(auto)
operator()( const ENDF::resolved::ReichMoore& rm, Functor&& callback ) const {
  switch( rm.NRO() ){
  case 0:
    switch( rm.NAPS() ){
    case 0:
      return callback( build( rm,
                              Scattering{},
                              channelRadius( rm.lValues().front().AWRI() ),
                              radius( rm.AP() ),
                              false ) );
    case 1:
      return callback( build( rm, Both{}, radius( rm.AP() ), true ) );
    }
  case 1:
    switch( rm.NAPS() ){
    case 0:
      return callback( build( rm,
                              Neither{},
                              channelRadius( rm.lValues().front().AWRI() ),
                              radius( rm.APE() ),
                              false ) );
    case 1:
      return callback( build( rm, Neither{}, radius( rm.APE() ), false ) );
    case 2:
      return callback( build( rm,
                              Channel{},
                              radius( rm.AP() ),
                              radius( rm.APE() ),
                              true ) );
    }
  }
}
