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
 * @file tests/mmstests/elasticity/TestIsotropicLinearElasticityLgDeform.hh
 *
 * @brief C++ TestIsotropicLinearElasticityLgDeform object
 *
 * C++ unit testing for IsotropicLinearElasticityLgDeform.
 */

#if !defined(pylith_mmstests_testisotropiclinearelasticityLgDeform_hh)
#define pylith_mmstests_testisotropiclinearelasticityLgDeform_hh

#include "TestElasticity.hh" // ISA TestElasticity

#include "pylith/materials/IsotropicLinearElasticityLgDeform.hh" // HOLDSA IsotropicLinearElasticityLgDeform

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearElasticityLgDeform;
    } // tests/mmstests
} // pylith

/// C++ unit testing for IsotropicLinearElasticityLgDeform
class pylith::mmstests::TestIsotropicLinearElasticityLgDeform : public pylith::mmstests::TestElasticity {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::materials::IsotropicLinearElasticityLgDeform* _rheology; ///< Rheology for testing.

}; // class TestIsotropicLinearElasticityLgDeform

#endif // pylith_mmstests_testisotropiclinearelasticityLgDeform_hh

// End of file
