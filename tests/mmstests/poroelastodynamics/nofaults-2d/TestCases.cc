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
// TEST_CASE("QuadTrig::TriP2::testDiscretization", "[QuadTrig][TriP2][discretization]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP2()).testDiscretization();
// }
// TEST_CASE("QuadTrig::TriP2::testResidual", "[QuadTrig][TriP2][residual]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP2()).testResidual();
// }

// TriP3
// TEST_CASE("QuadTrig::TriP3::testDiscretization", "[QuadTrig][TriP3][discretization]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP3()).testDiscretization();
// }
// TEST_CASE("QuadTrig::TriP3::testResidual", "[QuadTrig][TriP3][residual]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP3()).testResidual();
// }

// TriP4
// TEST_CASE("QuadTrig::TriP4::testDiscretization", "[QuadTrig][TriP4][discretization]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP4()).testDiscretization();
// }
// TEST_CASE("QuadTrig::TriP4::testResidual", "[QuadTrig][TriP4][residual]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::TriP4()).testResidual();
// }

// QuadQ2
TEST_CASE("QuadTrig::QuadQ2::testDiscretization", "[QuadTrig][QuadQ2][discretization]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ2()).testDiscretization();
}
TEST_CASE("QuadTrig::QuadQ2::testResidual", "[QuadTrig][QuadQ2][residual]") {
    pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ2()).testResidual();
}

// QuadQ3
// TEST_CASE("QuadTrig::QuadQ3::testDiscretization", "[QuadTrig][QuadQ3][discretization]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testDiscretization();
// }
// TEST_CASE("QuadTrig::QuadQ3::testResidual", "[QuadTrig][QuadQ3][residual]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ3()).testResidual();
// }

// QuadQ4
// TEST_CASE("QuadTrig::QuadQ4::testDiscretization", "[QuadTrig][QuadQ4][discretization]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ4()).testDiscretization();
// }
// TEST_CASE("QuadTrig::QuadQ4::testResidual", "[QuadTrig][QuadQ4][residual]") {
//     pylith::TestLinearPoroElasticity(pylith::QuadTrig::QuadQ4()).testResidual();
// }

// End of file
