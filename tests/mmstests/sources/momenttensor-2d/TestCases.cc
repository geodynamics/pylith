// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestMomentTensorSource
 *
 * Note: We only test discretization since the moment tensor source is a forcing
 * term that doesn't depend on the solution, so it has no Jacobian contribution.
 */

#include "TestMomentTensorSource.hh" // USES TestMomentTensorSource

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "MomentTensorSource2D.hh"

// TriP1
TEST_CASE("MomentTensorSource2D::TriP1::testDiscretization", "[MomentTensorSource2D][TriP1][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorSource2D::TriP1()).testDiscretization();
}

// TriP2
TEST_CASE("MomentTensorSource2D::TriP2::testDiscretization", "[MomentTensorSource2D][TriP2][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorSource2D::TriP2()).testDiscretization();
}

// QuadQ1
TEST_CASE("MomentTensorSource2D::QuadQ1::testDiscretization", "[MomentTensorSource2D][QuadQ1][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorSource2D::QuadQ1()).testDiscretization();
}

// QuadQ2
TEST_CASE("MomentTensorSource2D::QuadQ2::testDiscretization", "[MomentTensorSource2D][QuadQ2][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorSource2D::QuadQ2()).testDiscretization();
}

// ------------------------------------------------------------------------------------------------
#include "MomentTensorUniformStrain2D.hh"

// MomentTensorUniformStrain2D - Nonzero manufactured solution tests

// TriP1
TEST_CASE("MomentTensorUniformStrain2D::TriP1::testDiscretization", "[MomentTensorUniformStrain2D][TriP1][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorUniformStrain2D::TriP1()).testDiscretization();
}

// TriP2
TEST_CASE("MomentTensorUniformStrain2D::TriP2::testDiscretization", "[MomentTensorUniformStrain2D][TriP2][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorUniformStrain2D::TriP2()).testDiscretization();
}

// QuadQ1
TEST_CASE("MomentTensorUniformStrain2D::QuadQ1::testDiscretization", "[MomentTensorUniformStrain2D][QuadQ1][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorUniformStrain2D::QuadQ1()).testDiscretization();
}

// QuadQ2
TEST_CASE("MomentTensorUniformStrain2D::QuadQ2::testDiscretization", "[MomentTensorUniformStrain2D][QuadQ2][discretization]") {
    pylith::TestMomentTensorSource(pylith::MomentTensorUniformStrain2D::QuadQ2()).testDiscretization();
}


// End of file
