// =================================================================================================
// This code is part of SpatialData, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/scales/Scales.hh" // USES Scales

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath> // USES fabs()
#include <valarray> // USES std::valarray

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace scales {
        class TestScales;
    } // scales
} // pylith

class pylith::scales::TestScales {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructors.
    static
    void testConstructors(void);

    /// Test accessors.
    static
    void testAccessors(void);

    /// Test nondimensionalize() and dimensionalize().
    static
    void testNondimensionalize(void);

    /// Test Scalesie() and dimensionalize() with arrays.
    static
    void testNondimensionalizeArray(void);

}; // class TestScales

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestScales::testConstructors", "[TestScales]") {
    pylith::scales::TestScales::testConstructors();
}
TEST_CASE("TestScales::testAccessors", "[TestScales]") {
    pylith::scales::TestScales::testAccessors();
}
TEST_CASE("TestScales::testNondimensionalize", "[TestScales]") {
    pylith::scales::TestScales::testNondimensionalize();
}
TEST_CASE("TestScales::testNondimensionalizeArray", "[TestScales]") {
    pylith::scales::TestScales::testNondimensionalizeArray();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::scales::TestScales::testConstructors(void) {
    const double defaultLength = 1.0;
    const double defaultDisplacement = 1.0;
    const double defaultPressure = 1.0;
    const double defaultTime = 1.0;
    const double defaultTemperature = 1.0;

    Scales scales;
    CHECK(defaultLength == scales._length);
    CHECK(defaultDisplacement == scales._displacement);
    CHECK(defaultPressure == scales._rigidity);
    CHECK(defaultTime == scales._time);
    CHECK(defaultTemperature == scales._temperature);

    scales._length = 2.0;
    scales._displacement = 2.1;
    scales._rigidity = 3.0;
    scales._time = 4.0;
    scales._temperature = 5.0;
    Scales scalesCopy(scales);
    CHECK(scales._length == scalesCopy._length);
    CHECK(scales._displacement == scalesCopy._displacement);
    CHECK(scales._rigidity == scalesCopy._rigidity);
    CHECK(scales._time == scalesCopy._time);
    CHECK(scales._temperature == scalesCopy._temperature);

    Scales scalesAssign;
    scalesAssign = scales;
    CHECK(scales._length == scalesAssign._length);
    CHECK(scales._displacement == scalesAssign._displacement);
    CHECK(scales._rigidity == scalesAssign._rigidity);
    CHECK(scales._time == scalesAssign._time);
    CHECK(scales._temperature == scalesAssign._temperature);
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test accessors.
void
pylith::scales::TestScales::testAccessors(void) {
    const double length = 4.0;
    const double displacement = 4.1;
    const double pressure = 5.0;
    const double time = 6.0;
    const double temperature = 8.0;

    Scales scales;

    // Length scale
    INFO("Testing setting length scale.");
    scales.setLengthScale(length);
    CHECK(length == scales.getLengthScale());
    CHECK(1.0 == scales.getDisplacementScale());
    CHECK(1.0 == scales.getRigidityScale());
    CHECK(1.0 == scales.getTimeScale());
    CHECK(1.0 == scales.getTemperatureScale());

    // Displacement scale
    INFO("Testing setting displacement scale.");
    scales.setDisplacementScale(displacement);
    CHECK(length == scales.getLengthScale());
    CHECK(displacement == scales.getDisplacementScale());
    CHECK(1.0 == scales.getRigidityScale());
    CHECK(1.0 == scales.getTimeScale());
    CHECK(1.0 == scales.getTemperatureScale());

    // Pressure scale
    INFO("Testing setting pressure scale.");
    scales.setRigidityScale(pressure);
    CHECK(length == scales.getLengthScale());
    CHECK(displacement == scales.getDisplacementScale());
    CHECK(pressure == scales.getRigidityScale());
    CHECK(1.0 == scales.getTimeScale());
    CHECK(1.0 == scales.getTemperatureScale());

    // Time scale
    INFO("Testing setting time scale.");
    scales.setTimeScale(time);
    CHECK(length == scales.getLengthScale());
    CHECK(displacement == scales.getDisplacementScale());
    CHECK(pressure == scales.getRigidityScale());
    CHECK(time == scales.getTimeScale());
    CHECK(1.0 == scales.getTemperatureScale());

    // Temperature scale
    INFO("Testing setting temperature scale.");
    scales.setTemperatureScale(temperature);
    CHECK(length == scales.getLengthScale());
    CHECK(displacement == scales.getDisplacementScale());
    CHECK(pressure == scales.getRigidityScale());
    CHECK(time == scales.getTimeScale());
    CHECK(temperature == scales.getTemperatureScale());
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test nondimensionalize() and dimensionalize().
void
pylith::scales::TestScales::testNondimensionalize(void) {
    const double scale = 4.0;
    const double value = 3.0;
    const double valueE = 0.75;

    Scales scales;
    const double tolerance = 1.0e-6;
    CHECK_THAT(scales.nondimensionalize(value, scale), Catch::Matchers::WithinAbs(valueE, tolerance));
    CHECK_THAT(scales.dimensionalize(valueE, scale), Catch::Matchers::WithinAbs(value, tolerance));
} // testNondimensionalize


// ------------------------------------------------------------------------------------------------
// Test nondimensionalize() and dimensionalize() with arrays.
void
pylith::scales::TestScales::testNondimensionalizeArray(void) {
    const double scale = 10.0;
    const size_t nvalues = 3;
    const double values[nvalues] = { 2.0, 5.0, 7.0 };
    const double valuesE[nvalues] = { 0.2, 0.5, 0.7 };

    Scales scales;

    std::valarray<double> v(values, nvalues);
    scales.nondimensionalize(&v[0], nvalues, scale);
    const double tolerance = 1.0e-6;
    for (size_t i = 0; i < nvalues; ++i) {
        CHECK_THAT(v[i], Catch::Matchers::WithinAbs(valuesE[i], tolerance));
    } // for

    v = std::valarray<double>(valuesE, nvalues);
    scales.dimensionalize(&v[0], nvalues, scale);
    for (size_t i = 0; i < nvalues; ++i) {
        CHECK_THAT(v[i], Catch::Matchers::WithinAbs(values[i], tolerance));
    } // for
} // testNondimensionalizeArray


// End of file
