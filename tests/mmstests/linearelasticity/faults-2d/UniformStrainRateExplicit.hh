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

#include "TestFaultKin.hh" // USES TestFaultKin_Data

namespace pylith {
    class UniformStrainRateExplicit;
}

class pylith::UniformStrainRateExplicit {
public:

    // Data factory methods
    static TestFaultKin_Data* TriP1(void);

    static TestFaultKin_Data* TriP2(void);

    static TestFaultKin_Data* TriP3(void);

    static TestFaultKin_Data* TriP4(void);

    static TestFaultKin_Data* QuadQ1(void);

    static TestFaultKin_Data* QuadQ2(void);

    static TestFaultKin_Data* QuadQ3(void);

    static TestFaultKin_Data* QuadQ4(void);

private:

    UniformStrainRateExplicit(void); ///< Not implemented
}; // UniformStrainRateExplicit

// End of file
