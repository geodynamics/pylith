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

/** Test cases for TestLinearElasticity
 */

#include "TestLinearElasticity.hh" // USES TestLineaerElasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "UniformStrain3D.hh"
// TetP1
TEST_CASE("UniformStrain3D::TetP1::testDiscretization", "[UniformStrain3D][TetP1][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP1()).testDiscretization();
}
TEST_CASE("UniformStrain3D::TetP1::testResidual", "[UniformStrain3D][TetP1][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP1()).testResidual();
}
TEST_CASE("UniformStrain3D::TetP1::testJacobianTaylorSeries", "[UniformStrain3D][TetP1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain3D::TetP1::testJacobianFiniteDiff", "[UniformStrain3D][TetP1][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP1()).testJacobianFiniteDiff();
}

// TetP2
TEST_CASE("UniformStrain3D::TetP2::testDiscretization", "[UniformStrain3D][TetP2][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP2()).testDiscretization();
}
TEST_CASE("UniformStrain3D::TetP2::testResidual", "[UniformStrain3D][TetP2][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP2()).testResidual();
}
TEST_CASE("UniformStrain3D::TetP2::testJacobianTaylorSeries", "[UniformStrain3D][TetP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain3D::TetP2::testJacobianFiniteDiff", "[UniformStrain3D][TetP2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP2()).testJacobianFiniteDiff();
}

// TetP3
TEST_CASE("UniformStrain3D::TetP3::testDiscretization", "[UniformStrain3D][TetP3][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP3()).testDiscretization();
}
TEST_CASE("UniformStrain3D::TetP3::testResidual", "[UniformStrain3D][TetP3][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP3()).testResidual();
}
#if 0
TEST_CASE("UniformStrain3D::TetP3::testJacobianTaylorSeries", "[UniformStrain3D][TetP3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain3D::TetP3::testJacobianFiniteDiff", "[UniformStrain3D][TetP3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::TetP3()).testJacobianFiniteDiff();
}
#endif

// HexQ1
TEST_CASE("UniformStrain3D::HexQ1::testDiscretization", "[UniformStrain3D][HexQ1][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ1()).testDiscretization();
}
TEST_CASE("UniformStrain3D::HexQ1::testResidual", "[UniformStrain3D][HexQ1][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ1()).testResidual();
}
TEST_CASE("UniformStrain3D::HexQ1::testJacobianTaylorSeries", "[UniformStrain3D][HexQ1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain3D::HexQ1::testJacobianFiniteDiff", "[UniformStrain3D][HexQ1][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ1()).testJacobianFiniteDiff();
}

// HexQ2
TEST_CASE("UniformStrain3D::HexQ2::testDiscretization", "[UniformStrain3D][HexQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ2()).testDiscretization();
}
TEST_CASE("UniformStrain3D::HexQ2::testResidual", "[UniformStrain3D][HexQ2][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ2()).testResidual();
}
TEST_CASE("UniformStrain3D::HexQ2::testJacobianTaylorSeries", "[UniformStrain3D][HexQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain3D::HexQ2::testJacobianFiniteDiff", "[UniformStrain3D][HexQ2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ2()).testJacobianFiniteDiff();
}

// HexQ3
TEST_CASE("UniformStrain3D::HexQ3::testDiscretization", "[UniformStrain3D][HexQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ3()).testDiscretization();
}
TEST_CASE("UniformStrain3D::HexQ3::testResidual", "[UniformStrain3D][HexQ3][residual]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ3()).testResidual();
}
#if 0
TEST_CASE("UniformStrain3D::HexQ3::testJacobianTaylorSeries", "[UniformStrain3D][HexQ3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformStrain3D::HexQ3::testJacobianFiniteDiff", "[UniformStrain3D][HexQ3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::UniformStrain3D::HexQ3()).testJacobianFiniteDiff();
}
#endif

// ------------------------------------------------------------------------------------------------
#include "Gravity3D.hh"
// TetP2
TEST_CASE("Gravity3D::TetP2::testDiscretization", "[Gravity3D][TetP2][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP2()).testDiscretization();
}
TEST_CASE("Gravity3D::TetP2::testResidual", "[Gravity3D][TetP2][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP2()).testResidual();
}
TEST_CASE("Gravity3D::TetP2::testJacobianTaylorSeries", "[Gravity3D][TetP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP2()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity3D::TetP2::testJacobianFiniteDiff", "[Gravity3D][TetP2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP2()).testJacobianFiniteDiff();
}

// TetP3
TEST_CASE("Gravity3D::TetP3::testDiscretization", "[Gravity3D][TetP3][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP3()).testDiscretization();
}
TEST_CASE("Gravity3D::TetP3::testResidual", "[Gravity3D][TetP3][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP3()).testResidual();
}
#if 0
TEST_CASE("Gravity3D::TetP3::testJacobianTaylorSeries", "[Gravity3D][TetP3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP3()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity3D::TetP3::testJacobianFiniteDiff", "[Gravity3D][TetP3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::TetP3()).testJacobianFiniteDiff();
}
#endif

// HexQ2
TEST_CASE("Gravity3D::HexQ2::testDiscretization", "[Gravity3D][HexQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ2()).testDiscretization();
}
TEST_CASE("Gravity3D::HexQ2::testResidual", "[Gravity3D][HexQ2][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ2()).testResidual();
}
TEST_CASE("Gravity3D::HexQ2::testJacobianTaylorSeries", "[Gravity3D][HexQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ2()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity3D::HexQ2::testJacobianFiniteDiff", "[Gravity3D][HexQ2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ2()).testJacobianFiniteDiff();
}

// HexQ3
TEST_CASE("Gravity3D::HexQ3::testDiscretization", "[Gravity3D][HexQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ3()).testDiscretization();
}
TEST_CASE("Gravity3D::HexQ3::testResidual", "[Gravity3D][HexQ3][residual]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ3()).testResidual();
}
#if 0
TEST_CASE("Gravity3D::HexQ3::testJacobianTaylorSeries", "[Gravity3D][HexQ3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ3()).testJacobianTaylorSeries();
}
TEST_CASE("Gravity3D::HexQ3::testJacobianFiniteDiff", "[Gravity3D][HexQ3][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::Gravity3D::HexQ3()).testJacobianFiniteDiff();
}
#endif

// ------------------------------------------------------------------------------------------------
#include "GravityRefState3D.hh"
// No finite difference check because solution matches and no Jacobian is formed.
// TetP1
TEST_CASE("GravityRefState3D::TetP1::testDiscretization", "[GravityRefState3D][TetP1][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP1()).testDiscretization();
}
TEST_CASE("GravityRefState3D::TetP1::testResidual", "[GravityRefState3D][TetP1][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP1()).testResidual();
}
TEST_CASE("GravityRefState3D::TetP1::testJacobianTaylorSeries", "[GravityRefState3D][TetP1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP1()).testJacobianTaylorSeries();
}

// TetP2
TEST_CASE("GravityRefState3D::TetP2::testDiscretization", "[GravityRefState3D][TetP2][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP2()).testDiscretization();
}
TEST_CASE("GravityRefState3D::TetP2::testResidual", "[GravityRefState3D][TetP2][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP2()).testResidual();
}
TEST_CASE("GravityRefState3D::TetP2::testJacobianTaylorSeries", "[GravityRefState3D][TetP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP2()).testJacobianTaylorSeries();
}

// TetP3
TEST_CASE("GravityRefState3D::TetP3::testDiscretization", "[GravityRefState3D][TetP3][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP3()).testDiscretization();
}
TEST_CASE("GravityRefState3D::TetP3::testResidual", "[GravityRefState3D][TetP3][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP3()).testResidual();
}
#if 0
TEST_CASE("GravityRefState3D::TetP3::testJacobianTaylorSeries", "[GravityRefState3D][TetP3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::TetP3()).testJacobianTaylorSeries();
}
#endif

// HexQ1
TEST_CASE("GravityRefState3D::HexQ1::testDiscretization", "[GravityRefState3D][HexQ1][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ1()).testDiscretization();
}
TEST_CASE("GravityRefState3D::HexQ1::testResidual", "[GravityRefState3D][HexQ1][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ1()).testResidual();
}
TEST_CASE("GravityRefState3D::HexQ1::testJacobianTaylorSeries", "[GravityRefState3D][HexQ1][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ1()).testJacobianTaylorSeries();
}

// HexQ2
TEST_CASE("GravityRefState3D::HexQ2::testDiscretization", "[GravityRefState3D][HexQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ2()).testDiscretization();
}
TEST_CASE("GravityRefState3D::HexQ2::testResidual", "[GravityRefState3D][HexQ2][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ2()).testResidual();
}
TEST_CASE("GravityRefState3D::HexQ2::testJacobianTaylorSeries", "[GravityRefState3D][HexQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ2()).testJacobianTaylorSeries();
}

// HexQ3
TEST_CASE("GravityRefState3D::HexQ3::testDiscretization", "[GravityRefState3D][HexQ3][discretization]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ3()).testDiscretization();
}
TEST_CASE("GravityRefState3D::HexQ3::testResidual", "[GravityRefState3D][HexQ3][residual]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ3()).testResidual();
}
#if 0
TEST_CASE("GravityRefState3D::HexQ3::testJacobianTaylorSeries", "[GravityRefState3D][HexQ3][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::GravityRefState3D::HexQ3()).testJacobianTaylorSeries();
}
#endif

// ------------------------------------------------------------------------------------------------
#include "BodyForce3D.hh"
// TetP2
TEST_CASE("BodyForce3D::TetP2::testDiscretization", "[BodyForce3D][TetP2][discretization]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::TetP2()).testDiscretization();
}
TEST_CASE("BodyForce3D::TetP2::testResidual", "[BodyForce3D][TetP2][residual]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::TetP2()).testResidual();
}
TEST_CASE("BodyForce3D::TetP2::testJacobianTaylorSeries", "[BodyForce3D][TetP2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::TetP2()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce3D::TetP2::testJacobianFiniteDiff", "[BodyForce3D][TetP2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::TetP2()).testJacobianFiniteDiff();
}

// HexQ2
TEST_CASE("BodyForce3D::HexQ2::testDiscretization", "[BodyForce3D][HexQ2][discretization]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::HexQ2()).testDiscretization();
}
TEST_CASE("BodyForce3D::HexQ2::testResidual", "[BodyForce3D][HexQ2][residual]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::HexQ2()).testResidual();
}
TEST_CASE("BodyForce3D::HexQ2::testJacobianTaylorSeries", "[BodyForce3D][HexQ2][Jacobian Taylor series]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::HexQ2()).testJacobianTaylorSeries();
}
TEST_CASE("BodyForce3D::HexQ2::testJacobianFiniteDiff", "[BodyForce3D][HexQ2][Jacobian finite difference]") {
    pylith::TestLinearElasticity(pylith::BodyForce3D::HexQ2()).testJacobianFiniteDiff();
}
