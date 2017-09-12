template< typename Functor >
decltype(auto)
operator()( const ENDF::resolved::SLBW& slbw, Functor&& callback ) const {
  switch( slbw.NRO() ){
  case 0:
    switch( slbw.NAPS() ){
    case 0:
      return callback( build( slbw,
                              channelRadius( slbw.LStates().front().AWRI() ),
                              radius( slbw.AP() ) ) );
    case 1:
      return callback( build( slbw, radius( slbw.AP() ) ) );
    }
  case 1:
    switch( slbw.NAPS() ){
    case 0:
      return callback( build( slbw,
                              channelRadius( slbw.LStates().front().AWRI() ),
                              radius( slbw.APE() ) ) );
    case 1:
      return callback( build( slbw, radius( slbw.APE() ) ) );
    case 2:
      return callback( build( slbw,
                              radius( slbw.AP() ),
                              radius( slbw.APE() ) ) );
    }
  } 
}
