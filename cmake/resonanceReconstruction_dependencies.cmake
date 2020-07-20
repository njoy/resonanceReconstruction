include( FetchContent )


########################################################################
# Forward declarations
########################################################################


########################################################################
# Declare project dependencies
########################################################################

FetchContent_Declare( ENDFtk
    GIT_REPOSITORY  http://github.com/njoy/ENDFtk
    GIT_TAG         origin/build/fetchcontent
    )

FetchContent_Declare( eigen-adapter
    GIT_REPOSITORY  http://github.com/njoy/eigen-adapter
    GIT_TAG         origin/master
    )

FetchContent_Declare( interpolation
    GIT_REPOSITORY  http://github.com/njoy/interpolation
    GIT_TAG         origin/build/fetchcontent
    )

FetchContent_Declare( dimwits
    GIT_REPOSITORY  http://github.com/njoy/DimensionalAnalysis
    GIT_TAG         origin/build/fetchcontent
    )

FetchContent_Declare( elementary
    GIT_REPOSITORY  http://github.com/njoy/elementary
    GIT_TAG         origin/master
    )


########################################################################
# Load dependencies
########################################################################

FetchContent_MakeAvailable(
    ENDFtk
    eigen-adapter
    interpolation
    dimwits
    elementary
    )
