// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestLinearPoroelasticity
 */

#include "TestLinearPoroelasticity.hh" // USES TestLineaerPoroelasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "PressureGradient.hh"

// TriP2P1P1
TEST_CASE("PressureGradient::TriP2P1P1::testDiscretization", "[PressureGradient][TriP2P1P1][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP2P1P1::testResidual", "[PressureGradient][TriP2P1P1][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testResidual();
}
TEST_CASE("PressureGradient::TriP2P1P1::testJacobianTaylorSeries", "[PressureGradient][TriP2P1P1][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP2P1P1::testJacobianFiniteDiff", "[PressureGradient][TriP2P1P1][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testJacobianFiniteDiff();
}

// TriP3P2P2
TEST_CASE("PressureGradient::TriP3P2P2::testDiscretization", "[PressureGradient][TriP3P2P2][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP3P2P2::testResidual", "[PressureGradient][TriP3P2P2][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testResidual();
}
TEST_CASE("PressureGradient::TriP3P2P2::testJacobianTaylorSeries", "[PressureGradient][TriP3P2P2][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP3P2P2::testJacobianFiniteDiff", "[PressureGradient][TriP3P2P2][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testDiscretization", "[PressureGradient][QuadQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testResidual", "[PressureGradient][QuadQ2Q1Q1][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testJacobianTaylorSeries", "[PressureGradient][QuadQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testJacobianFiniteDiff", "[PressureGradient][QuadQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testDiscretization", "[PressureGradient][QuadQ3Q2Q2][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testResidual", "[PressureGradient][QuadQ3Q2Q2][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testJacobianTaylorSeries", "[PressureGradient][QuadQ3Q2Q2][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testJacobianFiniteDiff", "[PressureGradient][QuadQ3Q2Q2][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testJacobianFiniteDiff();
}

// End of file
