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

/** Test cases for TestLinearPoroElasticity
 */

#include "TestLinearPoroElasticity.hh" // USES TestLineaerPoroElasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "QuadTrig.hh"
// TriP2
TEST_CASE("QuadTrig::TriP2::testDiscretization", "[QuadTrig][TriP2][discretization]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP2()).testDiscretization();
}
TEST_CASE("QuadTrig::TriP2::testResidual", "[QuadTrig][TriP2][residual]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP2()).testResidual();
}
TEST_CASE("QuadTrig::TriP2::testJacobianTaylorSeries", "[QuadTrig][TriP2][Jacobian Taylor series]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("QuadTrig::TriP2::testJacobianFiniteDiff", "[QuadTrig][TriP2][Jacobian finite difference]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("QuadTrig::TriP3::testDiscretization", "[QuadTrig][TriP3][discretization]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP3()).testDiscretization();
}
TEST_CASE("QuadTrig::TriP3::testResidual", "[QuadTrig][TriP3][residual]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP3()).testResidual();
}
TEST_CASE("QuadTrig::TriP3::testJacobianTaylorSeries", "[QuadTrig][TriP3][Jacobian Taylor series]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("QuadTrig::TriP3::testJacobianFiniteDiff", "[QuadTrig][TriP3][Jacobian finite difference]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP3()).testJacobianFiniteDiff();
}

// TriP4
TEST_CASE("QuadTrig::TriP4::testDiscretization", "[QuadTrig][TriP4][discretization]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP4()).testDiscretization();
}
TEST_CASE("QuadTrig::TriP4::testResidual", "[QuadTrig][TriP4][residual]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP4()).testResidual();
}
TEST_CASE("QuadTrig::TriP4::testJacobianTaylorSeries", "[QuadTrig][TriP4][Jacobian Taylor series]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP4()).testJacobianTaylorSeries();
}
TEST_CASE("QuadTrig::TriP4::testJacobianFiniteDiff", "[QuadTrig][TriP4][Jacobian finite difference]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP4()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("QuadTrig::QuadQ3::testDiscretization", "[QuadTrig][QuadQ3][discretization]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testDiscretization();
}
TEST_CASE("QuadTrig::QuadQ3::testResidual", "[QuadTrig][QuadQ3][residual]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testResidual();
}
TEST_CASE("QuadTrig::QuadQ3::testJacobianTaylorSeries", "[QuadTrig][QuadQ3][Jacobian Taylor series]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("QuadTrig::QuadQ3::testJacobianFiniteDiff", "[QuadTrig][QuadQ3][Jacobian finite difference]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("QuadTrig::QuadQ3::testDiscretization", "[QuadTrig][QuadQ3][discretization]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testDiscretization();
}
TEST_CASE("QuadTrig::QuadQ3::testResidual", "[QuadTrig][QuadQ3][residual]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testResidual();
}
TEST_CASE("QuadTrig::QuadQ3::testJacobianTaylorSeries", "[QuadTrig][QuadQ3][Jacobian Taylor series]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("QuadTrig::QuadQ3::testJacobianFiniteDiff", "[QuadTrig][QuadQ3][Jacobian finite difference]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testJacobianFiniteDiff();
}

// QuadQ4
TEST_CASE("QuadTrig::QuadQ4::testDiscretization", "[QuadTrig][QuadQ4][discretization]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ4()).testDiscretization();
}
TEST_CASE("QuadTrig::QuadQ4::testResidual", "[QuadTrig][QuadQ4][residual]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ4()).testResidual();
}
TEST_CASE("QuadTrig::QuadQ4::testJacobianTaylorSeries", "[QuadTrig][QuadQ4][Jacobian Taylor series]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ4()).testJacobianTaylorSeries();
}
TEST_CASE("QuadTrig::QuadQ4::testJacobianFiniteDiff", "[QuadTrig][QuadQ4][Jacobian finite difference]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ4()).testJacobianFiniteDiff();
}

// End of file
