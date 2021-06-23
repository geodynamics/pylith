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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/mmstests/poroelasticity/TestIsotropicLinearPoroelasticity.hh
 *
 * @brief C++ TestIsotropicLinearPoroelasticity object
 *
 * C++ unit testing for IsotropicLinearPoroelasticity.
 */

#if !defined(pylith_mmstests_testisotropiclinearporoelasticity_hh)
#define pylith_mmstests_testisotropiclinearporoelasticity_hh

#include "TestPoroelasticity.hh" // ISA TestPoroelasticity

#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // HOLDSA IsotropicLinearPoroelasticity

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearPoroelasticity;
    } // tests/mmstests
} // pylith

/// C++ unit testing for IsotropicLinearPoroelasticity
class pylith::mmstests::TestIsotropicLinearPoroelasticity : public pylith::mmstests::TestPoroelasticity {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::materials::IsotropicLinearPoroelasticity* _rheology; ///< Rheology for testing.

}; // class TestIsotropicLinearPoroelasticity

#endif // pylith_mmstests_testisotropiclinearporoelasticity_hh

// End of file
