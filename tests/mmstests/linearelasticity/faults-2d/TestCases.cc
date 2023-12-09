// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestFaultKin
 */

#include "TestFaultKin.hh" // USES TestLineaerElasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "TwoBlocksStatic.hh"
// TriP1
TEST_CASE("TwoBlocksStatic::TriP1::testDiscretization", "[TwoBlocksStatic][TriP1][discretization]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP1()).testDiscretization();
}
TEST_CASE("TwoBlocksStatic::TriP1::testResidual", "[TwoBlocksStatic][TriP1][residual]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP1()).testResidual();
}
TEST_CASE("TwoBlocksStatic::TriP1::testJacobianTaylorSeries", "[TwoBlocksStatic][TriP1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStatic::TriP1::testJacobianFiniteDiff", "[TwoBlocksStatic][TriP1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("TwoBlocksStatic::TriP2::testDiscretization", "[TwoBlocksStatic][TriP2][discretization]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP2()).testDiscretization();
}
TEST_CASE("TwoBlocksStatic::TriP2::testResidual", "[TwoBlocksStatic][TriP2][residual]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP2()).testResidual();
}
TEST_CASE("TwoBlocksStatic::TriP2::testJacobianTaylorSeries", "[TwoBlocksStatic][TriP2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStatic::TriP2::testJacobianFiniteDiff", "[TwoBlocksStatic][TriP2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("TwoBlocksStatic::TriP3::testDiscretization", "[TwoBlocksStatic][TriP3][discretization]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP3()).testDiscretization();
}
TEST_CASE("TwoBlocksStatic::TriP3::testResidual", "[TwoBlocksStatic][TriP3][residual]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP3()).testResidual();
}
TEST_CASE("TwoBlocksStatic::TriP3::testJacobianTaylorSeries", "[TwoBlocksStatic][TriP3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStatic::TriP3::testJacobianFiniteDiff", "[TwoBlocksStatic][TriP3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("TwoBlocksStatic::QuadQ1::testDiscretization", "[TwoBlocksStatic][QuadQ1][discretization]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ1()).testDiscretization();
}
TEST_CASE("TwoBlocksStatic::QuadQ1::testResidual", "[TwoBlocksStatic][QuadQ1][residual]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ1()).testResidual();
}
TEST_CASE("TwoBlocksStatic::QuadQ1::testJacobianTaylorSeries", "[TwoBlocksStatic][QuadQ1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStatic::QuadQ1::testJacobianFiniteDiff", "[TwoBlocksStatic][QuadQ1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("TwoBlocksStatic::QuadQ2::testDiscretization", "[TwoBlocksStatic][QuadQ2][discretization]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ2()).testDiscretization();
}
TEST_CASE("TwoBlocksStatic::QuadQ2::testResidual", "[TwoBlocksStatic][QuadQ2][residual]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ2()).testResidual();
}
TEST_CASE("TwoBlocksStatic::QuadQ2::testJacobianTaylorSeries", "[TwoBlocksStatic][QuadQ2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStatic::QuadQ2::testJacobianFiniteDiff", "[TwoBlocksStatic][QuadQ2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("TwoBlocksStatic::QuadQ3::testDiscretization", "[TwoBlocksStatic][QuadQ3][discretization]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ3()).testDiscretization();
}
TEST_CASE("TwoBlocksStatic::QuadQ3::testResidual", "[TwoBlocksStatic][QuadQ3][residual]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ3()).testResidual();
}
TEST_CASE("TwoBlocksStatic::QuadQ3::testJacobianTaylorSeries", "[TwoBlocksStatic][QuadQ3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStatic::QuadQ3::testJacobianFiniteDiff", "[TwoBlocksStatic][QuadQ3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoBlocksStatic::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "ThreeBlocksStatic.hh"
// TriP1
TEST_CASE("ThreeBlocksStatic::TriP1::testDiscretization", "[ThreeBlocksStatic][TriP1][discretization]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP1()).testDiscretization();
}
TEST_CASE("ThreeBlocksStatic::TriP1::testResidual", "[ThreeBlocksStatic][TriP1][residual]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP1()).testResidual();
}
TEST_CASE("ThreeBlocksStatic::TriP1::testJacobianTaylorSeries", "[ThreeBlocksStatic][TriP1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("ThreeBlocksStatic::TriP1::testJacobianFiniteDiff", "[ThreeBlocksStatic][TriP1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("ThreeBlocksStatic::TriP2::testDiscretization", "[ThreeBlocksStatic][TriP2][discretization]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP2()).testDiscretization();
}
TEST_CASE("ThreeBlocksStatic::TriP2::testResidual", "[ThreeBlocksStatic][TriP2][residual]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP2()).testResidual();
}
TEST_CASE("ThreeBlocksStatic::TriP2::testJacobianTaylorSeries", "[ThreeBlocksStatic][TriP2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("ThreeBlocksStatic::TriP2::testJacobianFiniteDiff", "[ThreeBlocksStatic][TriP2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("ThreeBlocksStatic::TriP3::testDiscretization", "[ThreeBlocksStatic][TriP3][discretization]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP3()).testDiscretization();
}
TEST_CASE("ThreeBlocksStatic::TriP3::testResidual", "[ThreeBlocksStatic][TriP3][residual]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP3()).testResidual();
}
TEST_CASE("ThreeBlocksStatic::TriP3::testJacobianTaylorSeries", "[ThreeBlocksStatic][TriP3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("ThreeBlocksStatic::TriP3::testJacobianFiniteDiff", "[ThreeBlocksStatic][TriP3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("ThreeBlocksStatic::QuadQ1::testDiscretization", "[ThreeBlocksStatic][QuadQ1][discretization]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ1()).testDiscretization();
}
TEST_CASE("ThreeBlocksStatic::QuadQ1::testResidual", "[ThreeBlocksStatic][QuadQ1][residual]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ1()).testResidual();
}
TEST_CASE("ThreeBlocksStatic::QuadQ1::testJacobianTaylorSeries", "[ThreeBlocksStatic][QuadQ1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("ThreeBlocksStatic::QuadQ1::testJacobianFiniteDiff", "[ThreeBlocksStatic][QuadQ1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("ThreeBlocksStatic::QuadQ2::testDiscretization", "[ThreeBlocksStatic][QuadQ2][discretization]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ2()).testDiscretization();
}
TEST_CASE("ThreeBlocksStatic::QuadQ2::testResidual", "[ThreeBlocksStatic][QuadQ2][residual]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ2()).testResidual();
}
TEST_CASE("ThreeBlocksStatic::QuadQ2::testJacobianTaylorSeries", "[ThreeBlocksStatic][QuadQ2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("ThreeBlocksStatic::QuadQ2::testJacobianFiniteDiff", "[ThreeBlocksStatic][QuadQ2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("ThreeBlocksStatic::QuadQ3::testDiscretization", "[ThreeBlocksStatic][QuadQ3][discretization]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ3()).testDiscretization();
}
TEST_CASE("ThreeBlocksStatic::QuadQ3::testResidual", "[ThreeBlocksStatic][QuadQ3][residual]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ3()).testResidual();
}
TEST_CASE("ThreeBlocksStatic::QuadQ3::testJacobianTaylorSeries", "[ThreeBlocksStatic][QuadQ3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("ThreeBlocksStatic::QuadQ3::testJacobianFiniteDiff", "[ThreeBlocksStatic][QuadQ3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::ThreeBlocksStatic::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "OneFaultShearNoSlip.hh"
// TriP1
TEST_CASE("OneFaultShearNoSlip::TriP1::testDiscretization", "[OneFaultShearNoSlip][TriP1][discretization]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP1()).testDiscretization();
}
TEST_CASE("OneFaultShearNoSlip::TriP1::testResidual", "[OneFaultShearNoSlip][TriP1][residual]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP1()).testResidual();
}
TEST_CASE("OneFaultShearNoSlip::TriP1::testJacobianTaylorSeries", "[OneFaultShearNoSlip][TriP1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("OneFaultShearNoSlip::TriP1::testJacobianFiniteDiff", "[OneFaultShearNoSlip][TriP1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("OneFaultShearNoSlip::TriP2::testDiscretization", "[OneFaultShearNoSlip][TriP2][discretization]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP2()).testDiscretization();
}
TEST_CASE("OneFaultShearNoSlip::TriP2::testResidual", "[OneFaultShearNoSlip][TriP2][residual]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP2()).testResidual();
}
TEST_CASE("OneFaultShearNoSlip::TriP2::testJacobianTaylorSeries", "[OneFaultShearNoSlip][TriP2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("OneFaultShearNoSlip::TriP2::testJacobianFiniteDiff", "[OneFaultShearNoSlip][TriP2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("OneFaultShearNoSlip::TriP3::testDiscretization", "[OneFaultShearNoSlip][TriP3][discretization]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP3()).testDiscretization();
}
TEST_CASE("OneFaultShearNoSlip::TriP3::testResidual", "[OneFaultShearNoSlip][TriP3][residual]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP3()).testResidual();
}
TEST_CASE("OneFaultShearNoSlip::TriP3::testJacobianTaylorSeries", "[OneFaultShearNoSlip][TriP3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("OneFaultShearNoSlip::TriP3::testJacobianFiniteDiff", "[OneFaultShearNoSlip][TriP3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("OneFaultShearNoSlip::QuadQ1::testDiscretization", "[OneFaultShearNoSlip][QuadQ1][discretization]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ1()).testDiscretization();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ1::testResidual", "[OneFaultShearNoSlip][QuadQ1][residual]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ1()).testResidual();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ1::testJacobianTaylorSeries", "[OneFaultShearNoSlip][QuadQ1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ1::testJacobianFiniteDiff", "[OneFaultShearNoSlip][QuadQ1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("OneFaultShearNoSlip::QuadQ2::testDiscretization", "[OneFaultShearNoSlip][QuadQ2][discretization]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ2()).testDiscretization();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ2::testResidual", "[OneFaultShearNoSlip][QuadQ2][residual]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ2()).testResidual();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ2::testJacobianTaylorSeries", "[OneFaultShearNoSlip][QuadQ2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ2::testJacobianFiniteDiff", "[OneFaultShearNoSlip][QuadQ2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("OneFaultShearNoSlip::QuadQ3::testDiscretization", "[OneFaultShearNoSlip][QuadQ3][discretization]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ3()).testDiscretization();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ3::testResidual", "[OneFaultShearNoSlip][QuadQ3][residual]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ3()).testResidual();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ3::testJacobianTaylorSeries", "[OneFaultShearNoSlip][QuadQ3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("OneFaultShearNoSlip::QuadQ3::testJacobianFiniteDiff", "[OneFaultShearNoSlip][QuadQ3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::OneFaultShearNoSlip::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "TwoFaultsShearNoSlip.hh"
// TriP1
TEST_CASE("TwoFaultsShearNoSlip::TriP1::testDiscretization", "[TwoFaultsShearNoSlip][TriP1][discretization]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP1()).testDiscretization();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP1::testResidual", "[TwoFaultsShearNoSlip][TriP1][residual]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP1()).testResidual();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP1::testJacobianTaylorSeries", "[TwoFaultsShearNoSlip][TriP1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP1::testJacobianFiniteDiff", "[TwoFaultsShearNoSlip][TriP1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("TwoFaultsShearNoSlip::TriP2::testDiscretization", "[TwoFaultsShearNoSlip][TriP2][discretization]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP2()).testDiscretization();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP2::testResidual", "[TwoFaultsShearNoSlip][TriP2][residual]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP2()).testResidual();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP2::testJacobianTaylorSeries", "[TwoFaultsShearNoSlip][TriP2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP2::testJacobianFiniteDiff", "[TwoFaultsShearNoSlip][TriP2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("TwoFaultsShearNoSlip::TriP3::testDiscretization", "[TwoFaultsShearNoSlip][TriP3][discretization]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP3()).testDiscretization();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP3::testResidual", "[TwoFaultsShearNoSlip][TriP3][residual]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP3()).testResidual();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP3::testJacobianTaylorSeries", "[TwoFaultsShearNoSlip][TriP3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("TwoFaultsShearNoSlip::TriP3::testJacobianFiniteDiff", "[TwoFaultsShearNoSlip][TriP3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("TwoFaultsShearNoSlip::QuadQ1::testDiscretization", "[TwoFaultsShearNoSlip][QuadQ1][discretization]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ1()).testDiscretization();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ1::testResidual", "[TwoFaultsShearNoSlip][QuadQ1][residual]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ1()).testResidual();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ1::testJacobianTaylorSeries", "[TwoFaultsShearNoSlip][QuadQ1][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ1::testJacobianFiniteDiff", "[TwoFaultsShearNoSlip][QuadQ1][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("TwoFaultsShearNoSlip::QuadQ2::testDiscretization", "[TwoFaultsShearNoSlip][QuadQ2][discretization]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ2()).testDiscretization();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ2::testResidual", "[TwoFaultsShearNoSlip][QuadQ2][residual]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ2()).testResidual();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ2::testJacobianTaylorSeries", "[TwoFaultsShearNoSlip][QuadQ2][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ2::testJacobianFiniteDiff", "[TwoFaultsShearNoSlip][QuadQ2][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("TwoFaultsShearNoSlip::QuadQ3::testDiscretization", "[TwoFaultsShearNoSlip][QuadQ3][discretization]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ3()).testDiscretization();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ3::testResidual", "[TwoFaultsShearNoSlip][QuadQ3][residual]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ3()).testResidual();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ3::testJacobianTaylorSeries", "[TwoFaultsShearNoSlip][QuadQ3][Jacobian Taylor series]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("TwoFaultsShearNoSlip::QuadQ3::testJacobianFiniteDiff", "[TwoFaultsShearNoSlip][QuadQ3][Jacobian finite difference]") {
    pylith::TestFaultKin(pylith::TwoFaultsShearNoSlip::QuadQ3()).testJacobianFiniteDiff();
}

// End of file
