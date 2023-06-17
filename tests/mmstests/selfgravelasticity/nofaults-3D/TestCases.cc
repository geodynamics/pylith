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

/** Test cases for TestSelfGrav
 */

#include "TestSelfGrav.hh" // USES TestLineaerElasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "SelfGrav3D.hh"
// TetP2
TEST_CASE("SelfGrav3D::TetP2::testDiscretization", "[SelfGrav3D][TetP2][discretization]")
{
    pylith::TestSelfGrav(pylith::SelfGrav3D::TetP2()).testDiscretization();
}
TEST_CASE("SelfGrav3D::TetP2::testResidual", "[SelfGrav3D][TetP2][residual]")
{
    pylith::TestSelfGrav(pylith::SelfGrav3D::TetP2()).testResidual();
}
TEST_CASE("SelfGrav3D::TetP2::testJacobianTaylorSeries", "[SelfGrav3D][TetP2][Jacobian Taylor series]")
{
    pylith::TestSelfGrav(pylith::SelfGrav3D::TetP2()).testJacobianTaylorSeries();
}
TEST_CASE("SelfGrav3D::TetP2::testJacobianFiniteDiff", "[SelfGrav3D][TetP2][Jacobian finite difference]")
{
    pylith::TestSelfGrav(pylith::SelfGrav3D::TetP2()).testJacobianFiniteDiff();
}

// End of file
