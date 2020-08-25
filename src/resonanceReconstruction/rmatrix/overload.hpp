template < class... Types > struct overload : Types... {

  using Types::operator()...;
};
template < class... Types > overload( Types... ) -> overload< Types... >;
