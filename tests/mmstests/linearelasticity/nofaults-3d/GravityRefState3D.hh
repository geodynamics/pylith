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
    class GravityRefState3D;
}

class pylith::GravityRefState3D {
public:

    // Data factory methods

    static TestLinearElasticity_Data* TetP1(void);

    static TestLinearElasticity_Data* TetP2(void);

    static TestLinearElasticity_Data* TetP3(void);

    static TestLinearElasticity_Data* HexQ1(void);

    static TestLinearElasticity_Data* HexQ2(void);

    static TestLinearElasticity_Data* HexQ3(void);

private:

    GravityRefState3D(void); ///< Not implemented
}; // GravityRefState3D

// End of file
