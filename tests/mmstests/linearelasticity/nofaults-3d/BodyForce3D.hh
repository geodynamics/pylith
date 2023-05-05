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
    class BodyForce3D;
}

class pylith::BodyForce3D {
public:

    // Data factory methods

    static TestLinearElasticity_Data* TetP2(void);

    static TestLinearElasticity_Data* TetP3(void);

    static TestLinearElasticity_Data* HexQ2(void);

    static TestLinearElasticity_Data* HexQ3(void);

private:

    BodyForce3D(void); ///< Not implemented
}; // BodyForce3D

// End of file
