class Apply : protected singleLevel::Apply {
protected:
  #include "resonanceReconstruction/breitWigner/multiLevel/Apply/src/build.hpp"

public:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
  #include "resonanceReconstruction/breitWigner/multiLevel/Apply/src/call.hpp"
#pragma GCC diagnostic pop
};
