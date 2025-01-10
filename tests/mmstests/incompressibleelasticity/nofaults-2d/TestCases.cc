// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestIncompressibleElasticity
 */

#include "TestIncompressibleElasticity.hh" // USES TestLineaerElasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "UniformPressure2D.hh"
// TriP1
TEST_CASE("UniformPressure2D::TriP1::testDiscretization", "[UniformPressure2D][TriP1][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP1()).testDiscretization();
}
TEST_CASE("UniformPressure2D::TriP1::testResidual", "[UniformPressure2D][TriP1][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP1()).testResidual();
}
TEST_CASE("UniformPressure2D::TriP1::testJacobianTaylorSeries", "[UniformPressure2D][TriP1][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformPressure2D::TriP1::testJacobianFiniteDiff", "[UniformPressure2D][TriP1][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("UniformPressure2D::TriP2::testDiscretization", "[UniformPressure2D][TriP2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP2()).testDiscretization();
}
TEST_CASE("UniformPressure2D::TriP2::testResidual", "[UniformPressure2D][TriP2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP2()).testResidual();
}
TEST_CASE("UniformPressure2D::TriP2::testJacobianTaylorSeries", "[UniformPressure2D][TriP2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformPressure2D::TriP2::testJacobianFiniteDiff", "[UniformPressure2D][TriP2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("UniformPressure2D::TriP3::testDiscretization", "[UniformPressure2D][TriP3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP3()).testDiscretization();
}
TEST_CASE("UniformPressure2D::TriP3::testResidual", "[UniformPressure2D][TriP3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP3()).testResidual();
}
TEST_CASE("UniformPressure2D::TriP3::testJacobianTaylorSeries", "[UniformPressure2D][TriP3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformPressure2D::TriP3::testJacobianFiniteDiff", "[UniformPressure2D][TriP3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("UniformPressure2D::QuadQ1::testDiscretization", "[UniformPressure2D][QuadQ1][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ1()).testDiscretization();
}
TEST_CASE("UniformPressure2D::QuadQ1::testResidual", "[UniformPressure2D][QuadQ1][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ1()).testResidual();
}
TEST_CASE("UniformPressure2D::QuadQ1::testJacobianTaylorSeries", "[UniformPressure2D][QuadQ1][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformPressure2D::QuadQ1::testJacobianFiniteDiff", "[UniformPressure2D][QuadQ1][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("UniformPressure2D::QuadQ2::testDiscretization", "[UniformPressure2D][QuadQ2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ2()).testDiscretization();
}
TEST_CASE("UniformPressure2D::QuadQ2::testResidual", "[UniformPressure2D][QuadQ2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ2()).testResidual();
}
TEST_CASE("UniformPressure2D::QuadQ2::testJacobianTaylorSeries", "[UniformPressure2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformPressure2D::QuadQ2::testJacobianFiniteDiff", "[UniformPressure2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("UniformPressure2D::QuadQ3::testDiscretization", "[UniformPressure2D][QuadQ3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ3()).testDiscretization();
}
TEST_CASE("UniformPressure2D::QuadQ3::testResidual", "[UniformPressure2D][QuadQ3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ3()).testResidual();
}
TEST_CASE("UniformPressure2D::QuadQ3::testJacobianTaylorSeries", "[UniformPressure2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformPressure2D::QuadQ3::testJacobianFiniteDiff", "[UniformPressure2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformPressure2D::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "UniformShear2D.hh"
// TriP1
TEST_CASE("UniformShear2D::TriP1::testDiscretization", "[UniformShear2D][TriP1][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP1()).testDiscretization();
}
TEST_CASE("UniformShear2D::TriP1::testResidual", "[UniformShear2D][TriP1][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP1()).testResidual();
}
TEST_CASE("UniformShear2D::TriP1::testJacobianTaylorSeries", "[UniformShear2D][TriP1][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP1()).testJacobianTaylorSeries();
}
// No finite-difference check; zero pivot.

// TriP2
TEST_CASE("UniformShear2D::TriP2::testDiscretization", "[UniformShear2D][TriP2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP2()).testDiscretization();
}
TEST_CASE("UniformShear2D::TriP2::testResidual", "[UniformShear2D][TriP2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP2()).testResidual();
}
TEST_CASE("UniformShear2D::TriP2::testJacobianTaylorSeries", "[UniformShear2D][TriP2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformShear2D::TriP2::testJacobianFiniteDiff", "[UniformShear2D][TriP2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("UniformShear2D::TriP3::testDiscretization", "[UniformShear2D][TriP3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP3()).testDiscretization();
}
TEST_CASE("UniformShear2D::TriP3::testResidual", "[UniformShear2D][TriP3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP3()).testResidual();
}
TEST_CASE("UniformShear2D::TriP3::testJacobianTaylorSeries", "[UniformShear2D][TriP3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformShear2D::TriP3::testJacobianFiniteDiff", "[UniformShear2D][TriP3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("UniformShear2D::QuadQ1::testDiscretization", "[UniformShear2D][QuadQ1][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ1()).testDiscretization();
}
TEST_CASE("UniformShear2D::QuadQ1::testResidual", "[UniformShear2D][QuadQ1][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ1()).testResidual();
}
TEST_CASE("UniformShear2D::QuadQ1::testJacobianTaylorSeries", "[UniformShear2D][QuadQ1][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ1()).testJacobianTaylorSeries();
}
// No finite-difference check; zero pivot.

// QuadQ2
TEST_CASE("UniformShear2D::QuadQ2::testDiscretization", "[UniformShear2D][QuadQ2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ2()).testDiscretization();
}
TEST_CASE("UniformShear2D::QuadQ2::testResidual", "[UniformShear2D][QuadQ2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ2()).testResidual();
}
TEST_CASE("UniformShear2D::QuadQ2::testJacobianTaylorSeries", "[UniformShear2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformShear2D::QuadQ2::testJacobianFiniteDiff", "[UniformShear2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("UniformShear2D::QuadQ3::testDiscretization", "[UniformShear2D][QuadQ3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ3()).testDiscretization();
}
TEST_CASE("UniformShear2D::QuadQ3::testResidual", "[UniformShear2D][QuadQ3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ3()).testResidual();
}
TEST_CASE("UniformShear2D::QuadQ3::testJacobianTaylorSeries", "[UniformShear2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformShear2D::QuadQ3::testJacobianFiniteDiff", "[UniformShear2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::UniformShear2D::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "Gravity2D.hh"
// TriP1
TEST_CASE("Gravity2D::TriP1::testDiscretization", "[Gravity2D][TriP1][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP1()).testDiscretization();
}
TEST_CASE("Gravity2D::TriP1::testResidual", "[Gravity2D][TriP1][residual]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP1()).testResidual();
}
TEST_CASE("Gravity2D::TriP1::testJacobianTaylorSeries", "[Gravity2D][TriP1][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::TriP1::testJacobianFiniteDiff", "[Gravity2D][TriP1][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("Gravity2D::TriP2::testDiscretization", "[Gravity2D][TriP2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP2()).testDiscretization();
}
TEST_CASE("Gravity2D::TriP2::testResidual", "[Gravity2D][TriP2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP2()).testResidual();
}
TEST_CASE("Gravity2D::TriP2::testJacobianTaylorSeries", "[Gravity2D][TriP2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::TriP2::testJacobianFiniteDiff", "[Gravity2D][TriP2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("Gravity2D::TriP3::testDiscretization", "[Gravity2D][TriP3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP3()).testDiscretization();
}
TEST_CASE("Gravity2D::TriP3::testResidual", "[Gravity2D][TriP3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP3()).testResidual();
}
TEST_CASE("Gravity2D::TriP3::testJacobianTaylorSeries", "[Gravity2D][TriP3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::TriP3::testJacobianFiniteDiff", "[Gravity2D][TriP3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("Gravity2D::QuadQ2::testDiscretization", "[Gravity2D][QuadQ2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ2()).testDiscretization();
}
TEST_CASE("Gravity2D::QuadQ2::testResidual", "[Gravity2D][QuadQ2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ2()).testResidual();
}
TEST_CASE("Gravity2D::QuadQ2::testJacobianTaylorSeries", "[Gravity2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::QuadQ2::testJacobianFiniteDiff", "[Gravity2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("Gravity2D::QuadQ3::testDiscretization", "[Gravity2D][QuadQ3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ3()).testDiscretization();
}
TEST_CASE("Gravity2D::QuadQ3::testResidual", "[Gravity2D][QuadQ3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ3()).testResidual();
}
TEST_CASE("Gravity2D::QuadQ3::testJacobianTaylorSeries", "[Gravity2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::QuadQ3::testJacobianFiniteDiff", "[Gravity2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::Gravity2D::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "BodyForce2D.hh"
// TriP1
TEST_CASE("BodyForce2D::TriP1::testDiscretization", "[BodyForce2D][TriP1][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP1()).testDiscretization();
}
TEST_CASE("BodyForce2D::TriP1::testResidual", "[BodyForce2D][TriP1][residual]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP1()).testResidual();
}
TEST_CASE("BodyForce2D::TriP1::testJacobianTaylorSeries", "[BodyForce2D][TriP1][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::TriP1::testJacobianFiniteDiff", "[BodyForce2D][TriP1][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("BodyForce2D::TriP2::testDiscretization", "[BodyForce2D][TriP2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP2()).testDiscretization();
}
TEST_CASE("BodyForce2D::TriP2::testResidual", "[BodyForce2D][TriP2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP2()).testResidual();
}
TEST_CASE("BodyForce2D::TriP2::testJacobianTaylorSeries", "[BodyForce2D][TriP2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::TriP2::testJacobianFiniteDiff", "[BodyForce2D][TriP2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("BodyForce2D::TriP3::testDiscretization", "[BodyForce2D][TriP3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP3()).testDiscretization();
}
TEST_CASE("BodyForce2D::TriP3::testResidual", "[BodyForce2D][TriP3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP3()).testResidual();
}
TEST_CASE("BodyForce2D::TriP3::testJacobianTaylorSeries", "[BodyForce2D][TriP3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::TriP3::testJacobianFiniteDiff", "[BodyForce2D][TriP3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("BodyForce2D::QuadQ2::testDiscretization", "[BodyForce2D][QuadQ2][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ2()).testDiscretization();
}
TEST_CASE("BodyForce2D::QuadQ2::testResidual", "[BodyForce2D][QuadQ2][residual]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ2()).testResidual();
}
TEST_CASE("BodyForce2D::QuadQ2::testJacobianTaylorSeries", "[BodyForce2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::QuadQ2::testJacobianFiniteDiff", "[BodyForce2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("BodyForce2D::QuadQ3::testDiscretization", "[BodyForce2D][QuadQ3][discretization]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ3()).testDiscretization();
}
TEST_CASE("BodyForce2D::QuadQ3::testResidual", "[BodyForce2D][QuadQ3][residual]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ3()).testResidual();
}
TEST_CASE("BodyForce2D::QuadQ3::testJacobianTaylorSeries", "[BodyForce2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::QuadQ3::testJacobianFiniteDiff", "[BodyForce2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestIncompressibleElasticity(pylith::BodyForce2D::QuadQ3()).testJacobianFiniteDiff();
}

// End of file
