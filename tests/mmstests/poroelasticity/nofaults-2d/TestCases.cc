// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
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

// TriP2P1P1 w/state variables
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testDiscretization", "[PressureGradient][TriP2P1P1_StateVars][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testResidual", "[PressureGradient][TriP2P1P1_StateVars][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testJacobianTaylorSeries", "[PressureGradient][TriP2P1P1_StateVars][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testJacobianFiniteDiff", "[PressureGradient][TriP2P1P1_StateVars][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testJacobianFiniteDiff();
}

// TriP3P2P2 with state variables
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testDiscretization", "[PressureGradient][TriP3P2P2_StateVars][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testResidual", "[PressureGradient][TriP3P2P2_StateVars][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testJacobianTaylorSeries", "[PressureGradient][TriP3P2P2_StateVars][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testJacobianFiniteDiff", "[PressureGradient][TriP3P2P2_StateVars][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1 with state variables
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testDiscretization", "[PressureGradient][QuadQ2Q1Q1_StateVars][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testResidual", "[PressureGradient][QuadQ2Q1Q1_StateVars][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testJacobianTaylorSeries", "[PressureGradient][QuadQ2Q1Q1_StateVars][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testJacobianFiniteDiff", "[PressureGradient][QuadQ2Q1Q1_StateVars][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2 with state variables
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testDiscretization", "[PressureGradient][QuadQ3Q2Q2][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testResidual", "[PressureGradient][QuadQ3Q2Q2][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testJacobianTaylorSeries", "[PressureGradient][QuadQ3Q2Q2][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testJacobianFiniteDiff", "[PressureGradient][QuadQ3Q2Q2][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testJacobianFiniteDiff();
}

// End of file
