// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2019 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include "pylith/testing/TestDriver.cc"

int
main(int argc,
     char* argv[]) {
    pylith::testing::TestDriver driver;
    return driver.run(argc, argv);
} // main


// End of file
