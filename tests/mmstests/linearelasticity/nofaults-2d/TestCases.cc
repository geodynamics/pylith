// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestLinearElasticity
 */

#include "TestLinearElasticity.hh" // USES TestLineaerElasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "UniformStrain2D.hh"
// TriP1
TEST_CASE("UniformStrain2D::TriP1::testDiscretization", "[UniformStrain2D][TriP1][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP1()).testDiscretization();
}
TEST_CASE("UniformStrain2D::TriP1::testResidual", "[UniformStrain2D][TriP1][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP1()).testResidual();
}
TEST_CASE("UniformStrain2D::TriP1::testJacobianTaylorSeries", "[UniformStrain2D][TriP1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain2D::TriP1::testJacobianFiniteDiff", "[UniformStrain2D][TriP1][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("UniformStrain2D::TriP2::testDiscretization", "[UniformStrain2D][TriP2][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP2()).testDiscretization();
}
TEST_CASE("UniformStrain2D::TriP2::testResidual", "[UniformStrain2D][TriP2][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP2()).testResidual();
}
TEST_CASE("UniformStrain2D::TriP2::testJacobianTaylorSeries", "[UniformStrain2D][TriP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain2D::TriP2::testJacobianFiniteDiff", "[UniformStrain2D][TriP2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("UniformStrain2D::TriP3::testDiscretization", "[UniformStrain2D][TriP3][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP3()).testDiscretization();
}
TEST_CASE("UniformStrain2D::TriP3::testResidual", "[UniformStrain2D][TriP3][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP3()).testResidual();
}
TEST_CASE("UniformStrain2D::TriP3::testJacobianTaylorSeries", "[UniformStrain2D][TriP3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain2D::TriP3::testJacobianFiniteDiff", "[UniformStrain2D][TriP3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("UniformStrain2D::QuadQ1::testDiscretization", "[UniformStrain2D][QuadQ1][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ1()).testDiscretization();
}
TEST_CASE("UniformStrain2D::QuadQ1::testResidual", "[UniformStrain2D][QuadQ1][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ1()).testResidual();
}
TEST_CASE("UniformStrain2D::QuadQ1::testJacobianTaylorSeries", "[UniformStrain2D][QuadQ1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain2D::QuadQ1::testJacobianFiniteDiff", "[UniformStrain2D][QuadQ1][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("UniformStrain2D::QuadQ2::testDiscretization", "[UniformStrain2D][QuadQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ2()).testDiscretization();
}
TEST_CASE("UniformStrain2D::QuadQ2::testResidual", "[UniformStrain2D][QuadQ2][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ2()).testResidual();
}
TEST_CASE("UniformStrain2D::QuadQ2::testJacobianTaylorSeries", "[UniformStrain2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain2D::QuadQ2::testJacobianFiniteDiff", "[UniformStrain2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("UniformStrain2D::QuadQ3::testDiscretization", "[UniformStrain2D][QuadQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ3()).testDiscretization();
}
TEST_CASE("UniformStrain2D::QuadQ3::testResidual", "[UniformStrain2D][QuadQ3][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ3()).testResidual();
}
TEST_CASE("UniformStrain2D::QuadQ3::testJacobianTaylorSeries", "[UniformStrain2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain2D::QuadQ3::testJacobianFiniteDiff", "[UniformStrain2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain2D::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "Gravity2D.hh"
// TriP2
TEST_CASE("Gravity2D::TriP2::testDiscretization", "[Gravity2D][TriP2][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP2()).testDiscretization();
}
TEST_CASE("Gravity2D::TriP2::testResidual", "[Gravity2D][TriP2][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP2()).testResidual();
}
TEST_CASE("Gravity2D::TriP2::testJacobianTaylorSeries", "[Gravity2D][TriP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::TriP2::testJacobianFiniteDiff", "[Gravity2D][TriP2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("Gravity2D::TriP3::testDiscretization", "[Gravity2D][TriP3][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP3()).testDiscretization();
}
TEST_CASE("Gravity2D::TriP3::testResidual", "[Gravity2D][TriP3][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP3()).testResidual();
}
TEST_CASE("Gravity2D::TriP3::testJacobianTaylorSeries", "[Gravity2D][TriP3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::TriP3::testJacobianFiniteDiff", "[Gravity2D][TriP3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("Gravity2D::QuadQ2::testDiscretization", "[Gravity2D][QuadQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ2()).testDiscretization();
}
TEST_CASE("Gravity2D::QuadQ2::testResidual", "[Gravity2D][QuadQ2][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ2()).testResidual();
}
TEST_CASE("Gravity2D::QuadQ2::testJacobianTaylorSeries", "[Gravity2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::QuadQ2::testJacobianFiniteDiff", "[Gravity2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("Gravity2D::QuadQ3::testDiscretization", "[Gravity2D][QuadQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ3()).testDiscretization();
}
TEST_CASE("Gravity2D::QuadQ3::testResidual", "[Gravity2D][QuadQ3][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ3()).testResidual();
}
TEST_CASE("Gravity2D::QuadQ3::testJacobianTaylorSeries", "[Gravity2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity2D::QuadQ3::testJacobianFiniteDiff", "[Gravity2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity2D::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "GravityRefState2D.hh"
// No finite difference check because solution matches and no Jacobian is formed.
// TriP1
TEST_CASE("GravityRefState2D::TriP1::testDiscretization", "[GravityRefState2D][TriP1][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP1()).testDiscretization();
}
TEST_CASE("GravityRefState2D::TriP1::testResidual", "[GravityRefState2D][TriP1][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP1()).testResidual();
}
TEST_CASE("GravityRefState2D::TriP1::testJacobianTaylorSeries", "[GravityRefState2D][TriP1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP1()).testJacobianTaylorSeries();
}

// TriP2
TEST_CASE("GravityRefState2D::TriP2::testDiscretization", "[GravityRefState2D][TriP2][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP2()).testDiscretization();
}
TEST_CASE("GravityRefState2D::TriP2::testResidual", "[GravityRefState2D][TriP2][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP2()).testResidual();
}
TEST_CASE("GravityRefState2D::TriP2::testJacobianTaylorSeries", "[GravityRefState2D][TriP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP2()).testJacobianTaylorSeries();
}

// TriP3
TEST_CASE("GravityRefState2D::TriP3::testDiscretization", "[GravityRefState2D][TriP3][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP3()).testDiscretization();
}
TEST_CASE("GravityRefState2D::TriP3::testResidual", "[GravityRefState2D][TriP3][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP3()).testResidual();
}
TEST_CASE("GravityRefState2D::TriP3::testJacobianTaylorSeries", "[GravityRefState2D][TriP3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::TriP3()).testJacobianTaylorSeries();
}

// QuadQ1
TEST_CASE("GravityRefState2D::QuadQ1::testDiscretization", "[GravityRefState2D][QuadQ1][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ1()).testDiscretization();
}
TEST_CASE("GravityRefState2D::QuadQ1::testResidual", "[GravityRefState2D][QuadQ1][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ1()).testResidual();
}
TEST_CASE("GravityRefState2D::QuadQ1::testJacobianTaylorSeries", "[GravityRefState2D][QuadQ1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ1()).testJacobianTaylorSeries();
}

// QuadQ2
TEST_CASE("GravityRefState2D::QuadQ2::testDiscretization", "[GravityRefState2D][QuadQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ2()).testDiscretization();
}
TEST_CASE("GravityRefState2D::QuadQ2::testResidual", "[GravityRefState2D][QuadQ2][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ2()).testResidual();
}
TEST_CASE("GravityRefState2D::QuadQ2::testJacobianTaylorSeries", "[GravityRefState2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ2()).testJacobianTaylorSeries();
}

// QuadQ3
TEST_CASE("GravityRefState2D::QuadQ3::testDiscretization", "[GravityRefState2D][QuadQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ3()).testDiscretization();
}
TEST_CASE("GravityRefState2D::QuadQ3::testResidual", "[GravityRefState2D][QuadQ3][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ3()).testResidual();
}
TEST_CASE("GravityRefState2D::QuadQ3::testJacobianTaylorSeries", "[GravityRefState2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState2D::QuadQ3()).testJacobianTaylorSeries();
}

// ------------------------------------------------------------------------------------------------
#include "BodyForce2D.hh"
// TriP2
TEST_CASE("BodyForce2D::TriP2::testDiscretization", "[BodyForce2D][TriP2][discretization]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP2()).testDiscretization();
}
TEST_CASE("BodyForce2D::TriP2::testResidual", "[BodyForce2D][TriP2][residual]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP2()).testResidual();
}
TEST_CASE("BodyForce2D::TriP2::testJacobianTaylorSeries", "[BodyForce2D][TriP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::TriP2::testJacobianFiniteDiff", "[BodyForce2D][TriP2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("BodyForce2D::TriP3::testDiscretization", "[BodyForce2D][TriP3][discretization]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP3()).testDiscretization();
}
TEST_CASE("BodyForce2D::TriP3::testResidual", "[BodyForce2D][TriP3][residual]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP3()).testResidual();
}
TEST_CASE("BodyForce2D::TriP3::testJacobianTaylorSeries", "[BodyForce2D][TriP3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::TriP3::testJacobianFiniteDiff", "[BodyForce2D][TriP3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("BodyForce2D::QuadQ2::testDiscretization", "[BodyForce2D][QuadQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ2()).testDiscretization();
}
TEST_CASE("BodyForce2D::QuadQ2::testResidual", "[BodyForce2D][QuadQ2][residual]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ2()).testResidual();
}
TEST_CASE("BodyForce2D::QuadQ2::testJacobianTaylorSeries", "[BodyForce2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::QuadQ2::testJacobianFiniteDiff", "[BodyForce2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("BodyForce2D::QuadQ3::testDiscretization", "[BodyForce2D][QuadQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ3()).testDiscretization();
}
TEST_CASE("BodyForce2D::QuadQ3::testResidual", "[BodyForce2D][QuadQ3][residual]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ3()).testResidual();
}
TEST_CASE("BodyForce2D::QuadQ3::testJacobianTaylorSeries", "[BodyForce2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce2D::QuadQ3::testJacobianFiniteDiff", "[BodyForce2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::BodyForce2D::QuadQ3()).testJacobianFiniteDiff();
}

// ------------------------------------------------------------------------------------------------
#include "RigidBodyAcc2D.hh"
// TriP1
TEST_CASE("RigidBodyAcc2D::TriP1::testDiscretization", "[RigidBodyAcc2D][TriP1][discretization]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::TriP1()).testDiscretization();
}
TEST_CASE("RigidBodyAcc2D::TriP1::testResidual", "[RigidBodyAcc2D][TriP1][residual]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::TriP1()).testResidual();
}

// TriP2
TEST_CASE("RigidBodyAcc2D::TriP2::testDiscretization", "[RigidBodyAcc2D][TriP2][discretization]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::TriP2()).testDiscretization();
}
TEST_CASE("RigidBodyAcc2D::TriP2::testResidual", "[RigidBodyAcc2D][TriP2][residual]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::TriP2()).testResidual();
}

// TriP3
TEST_CASE("RigidBodyAcc2D::TriP3::testDiscretization", "[RigidBodyAcc2D][TriP3][discretization]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::TriP3()).testDiscretization();
}
TEST_CASE("RigidBodyAcc2D::TriP3::testResidual", "[RigidBodyAcc2D][TriP3][residual]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::TriP3()).testResidual();
}

// QuadQ1
TEST_CASE("RigidBodyAcc2D::QuadQ1::testDiscretization", "[RigidBodyAcc2D][QuadQ1][discretization]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::QuadQ1()).testDiscretization();
}
TEST_CASE("RigidBodyAcc2D::QuadQ1::testResidual", "[RigidBodyAcc2D][QuadQ1][residual]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::QuadQ1()).testResidual();
}

// QuadQ2
TEST_CASE("RigidBodyAcc2D::QuadQ2::testDiscretization", "[RigidBodyAcc2D][QuadQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::QuadQ2()).testDiscretization();
}
TEST_CASE("RigidBodyAcc2D::QuadQ2::testResidual", "[RigidBodyAcc2D][QuadQ2][residual]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::QuadQ2()).testResidual();
}

// QuadQ3
TEST_CASE("RigidBodyAcc2D::QuadQ3::testDiscretization", "[RigidBodyAcc2D][QuadQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::QuadQ3()).testDiscretization();
}
TEST_CASE("RigidBodyAcc2D::QuadQ3::testResidual", "[RigidBodyAcc2D][QuadQ3][residual]") {
    pylith::TestLinearElasticity(pylith::RigidBodyAcc2D::QuadQ3()).testResidual();
}

// ------------------------------------------------------------------------------------------------
#include "PlanePWave2D.hh"
// TriP1
TEST_CASE("PlanePWave2D::TriP1::testDiscretization", "[PlanePWave2D][TriP1][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP1()).testDiscretization();
}
TEST_CASE("PlanePWave2D::TriP1::testResidual", "[PlanePWave2D][TriP1][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP1()).testResidual();
}

// TriP2
TEST_CASE("PlanePWave2D::TriP2::testDiscretization", "[PlanePWave2D][TriP2][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP2()).testDiscretization();
}
TEST_CASE("PlanePWave2D::TriP2::testResidual", "[PlanePWave2D][TriP2][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP2()).testResidual();
}

// TriP3
TEST_CASE("PlanePWave2D::TriP3::testDiscretization", "[PlanePWave2D][TriP3][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP3()).testDiscretization();
}
TEST_CASE("PlanePWave2D::TriP3::testResidual", "[PlanePWave2D][TriP3][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP3()).testResidual();
}

// TriP4
TEST_CASE("PlanePWave2D::TriP4::testDiscretization", "[PlanePWave2D][TriP4][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP4()).testDiscretization();
}
TEST_CASE("PlanePWave2D::TriP4::testResidual", "[PlanePWave2D][TriP4][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::TriP4()).testResidual();
}

// QuadQ1
TEST_CASE("PlanePWave2D::QuadQ1::testDiscretization", "[PlanePWave2D][QuadQ1][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ1()).testDiscretization();
}
TEST_CASE("PlanePWave2D::QuadQ1::testResidual", "[PlanePWave2D][QuadQ1][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ1()).testResidual();
}

// QuadQ2
TEST_CASE("PlanePWave2D::QuadQ2::testDiscretization", "[PlanePWave2D][QuadQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ2()).testDiscretization();
}
TEST_CASE("PlanePWave2D::QuadQ2::testResidual", "[PlanePWave2D][QuadQ2][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ2()).testResidual();
}

// QuadQ3
TEST_CASE("PlanePWave2D::QuadQ3::testDiscretization", "[PlanePWave2D][QuadQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ3()).testDiscretization();
}
TEST_CASE("PlanePWave2D::QuadQ3::testResidual", "[PlanePWave2D][QuadQ3][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ3()).testResidual();
}

// QuadQ4
TEST_CASE("PlanePWave2D::QuadQ4::testDiscretization", "[PlanePWave2D][QuadQ4][discretization]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ4()).testDiscretization();
}
TEST_CASE("PlanePWave2D::QuadQ4::testResidual", "[PlanePWave2D][QuadQ4][residual]") {
    pylith::TestLinearElasticity(pylith::PlanePWave2D::QuadQ4()).testResidual();
}

// End of file
