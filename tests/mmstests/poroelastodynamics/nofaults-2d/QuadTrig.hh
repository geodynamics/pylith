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

#include "TestLinearPoroElasticity.hh" // USES TestPoroLinearElasticity_Data

namespace pylith {
    class QuadTrig;
}

class pylith::QuadTrig {
public:

    // Data factory methods
    static TestLinearPoroElasticity_Data* TriP2(void);

    static TestLinearPoroElasticity_Data* TriP3(void);

    static TestLinearPoroElasticity_Data* TriP4(void);

    static TestLinearPoroElasticity_Data* TriP5(void);

    static TestLinearPoroElasticity_Data* QuadQ2(void);

    static TestLinearPoroElasticity_Data* QuadQ2Distorted(void);

    static TestLinearPoroElasticity_Data* QuadQ3(void);

    static TestLinearPoroElasticity_Data* QuadQ4(void);

    static TestLinearPoroElasticity_Data* QuadQ5(void);

private:

    QuadTrig(void);
}; // QuadTrig

// End of file
