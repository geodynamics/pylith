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

#include "TestSelfGrav.hh" // USES TestSelfGrav_Data

namespace pylith {
    class SelfGrav3D;
}

class pylith::SelfGrav3D {
public:

    // Data factory methods

    static TestSelfGrav_Data *TetP2(void);

private:

    SelfGrav3D(void); ///< Not implemented
};                    // SelfGrav3D

// End of file
