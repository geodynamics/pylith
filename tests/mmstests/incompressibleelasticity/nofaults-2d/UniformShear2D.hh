// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIncompressibleElasticity.hh" // USES TestIncompressibleElasticity_Data

namespace pylith {
    class UniformShear2D;
}

class pylith::UniformShear2D {
public:

    // Data factory methods
    static TestIncompressibleElasticity_Data* TriP1(void);

    static TestIncompressibleElasticity_Data* TriP2(void);

    static TestIncompressibleElasticity_Data* TriP3(void);

    static TestIncompressibleElasticity_Data* TriP4(void);

    static TestIncompressibleElasticity_Data* QuadQ1(void);

    static TestIncompressibleElasticity_Data* QuadQ1Distorted(void);

    static TestIncompressibleElasticity_Data* QuadQ2(void);

    static TestIncompressibleElasticity_Data* QuadQ3(void);

    static TestIncompressibleElasticity_Data* QuadQ4(void);

private:

    UniformShear2D(void); ///< Not implemented
}; // UniformShear2D

// End of file
