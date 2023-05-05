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

#include "TestLinearElasticity.hh" // USES TestLinearElasticity_Data

namespace pylith {
    class Gravity3D;
}

class pylith::Gravity3D {
public:

    // Data factory methods

    static TestLinearElasticity_Data* TetP2(void);

    static TestLinearElasticity_Data* TetP3(void);

    static TestLinearElasticity_Data* HexQ2(void);

    static TestLinearElasticity_Data* HexQ3(void);

private:

    Gravity3D(void); ///< Not implemented
}; // Gravity3D

// End of file
