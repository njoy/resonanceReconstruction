template< typename Functor >
decltype(auto)
operator()( const ENDF::resolved::MLBW& mlbw, Functor&& callback ) const {
  switch( mlbw.NRO() ){
  case 0:
    switch( mlbw.NAPS() ){
    case 0:
      return callback( build( mlbw,
                              channelRadius( mlbw.lValues().front().AWRI() ),
                              radius( mlbw.AP() ) ) );
    case 1:
      return callback( build( mlbw, radius( mlbw.AP() ) ) );
    }
  case 1:
    switch( mlbw.NAPS() ){
    case 0:
      return callback( build( mlbw,
                              channelRadius( mlbw.lValues().front().AWRI() ),
                              radius( mlbw.APE() ) ) );
    case 1:
      return callback( build( mlbw, radius( mlbw.APE() ) ) );
    case 2:
      return callback( build( mlbw,
                              radius( mlbw.AP() ),
                              radius( mlbw.APE() ) ) );
    }
  }
}
