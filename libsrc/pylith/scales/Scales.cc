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

#include "Scales.hh" // implementation of class methods

#include "pylith/utils/journals.hh" // USES journal macros
#include "pylith/utils/Exceptions.hh" // USES Exception

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor
pylith::scales::Scales::Scales(void) :
    _length(1.0),
    _displacement(1.0),
    _rigidity(1.0),
    _time(1.0),
    _temperature(1.0) {}


// ----------------------------------------------------------------------
// Default destructor
pylith::scales::Scales::~Scales(void) {}


// ----------------------------------------------------------------------
// Copy constructor.
pylith::scales::Scales::Scales(const Scales& dim) :
    _length(dim._length),
    _displacement(dim._displacement),
    _rigidity(dim._rigidity),
    _time(dim._time),
    _temperature(dim._temperature) {}


// ----------------------------------------------------------------------
// Assignment operator.
const pylith::scales::Scales&
pylith::scales::Scales::operator=(const Scales& dim) {
    if (this != &dim) {
        _length = dim._length;
        _displacement = dim._displacement;
        _rigidity = dim._rigidity;
        _time = dim._time;
        _temperature = dim._temperature;
    } // if

    return *this;
} // operator=


// ----------------------------------------------------------------------
// Set value to nondimensionalize position.
void
pylith::scales::Scales::setLengthScale(const double value) {
    if (value <= 0.0) {
        PYLITH_ERROR(pylith::ValueError, pylith::journal::user_input,
                     "Length scale (" << value << ") must be positive.");
    } // if
    _length = value;
} // setLengthScale


// ----------------------------------------------------------------------
// Set value to nondimensionalize displacement.
void
pylith::scales::Scales::setDisplacementScale(const double value) {
    if (value <= 0.0) {
        PYLITH_ERROR(pylith::ValueError, pylith::journal::user_input,
                     "Displacement scale (" << value << ") must be positive.");
    } // if
    _displacement = value;
} // setDisplacementScale


// ----------------------------------------------------------------------
// Set value to nondimensionalize rigidity (elastic moduli).
void
pylith::scales::Scales::setRigidityScale(const double value) {
    if (value <= 0.0) {
        PYLITH_ERROR(pylith::ValueError, pylith::journal::user_input,
                     "Rigidity scale (" << value << ") must be positive.");
    } // if
    _rigidity = value;
} // setRigidityScale


// ----------------------------------------------------------------------
// Set value to nondimensionalize time scale in seconds (SI units).
void
pylith::scales::Scales::setTimeScale(const double value) {
    if (value <= 0.0) {
        PYLITH_ERROR(pylith::ValueError, pylith::journal::user_input,
                     "Time scale (" << value << ") must be positive.");
    } // if
    _time = value;
} // setTimeScale


// ----------------------------------------------------------------------
// Set value to nondimensionalize temperature scale in Kelvin (SI units).
void
pylith::scales::Scales::setTemperatureScale(const double value) {
    if (value <= 0.0) {
        PYLITH_ERROR(pylith::ValueError, pylith::journal::user_input,
                     "Temperature scale (" << value << ") must be positive.");
    } // if
    _temperature = value;
} // setTemperatureScale


// End of file
