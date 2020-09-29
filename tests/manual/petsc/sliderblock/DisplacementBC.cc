#include <portinfo>

#include "DisplacementBC.hh"

namespace _DisplacementBC {
    static const double velocity = 1.0;
}

// --------------------------------------------------------------------------------------------------
double
DisplacementBC::displacement(const double t) {
    return _DisplacementBC::velocity * t;
}


// End of file
