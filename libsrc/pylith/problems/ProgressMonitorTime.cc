// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ProgressMonitorTime.hh" // implementation of class methods

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/units/Parser.hh" // USES Parser

#include <ctime> // USES C time_t
#include <iomanip> // USES std::setw, etcm
#include <iostream> // USES std::ofstream
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::ProgressMonitorTime::ProgressMonitorTime(void) :
    _baseTime(1.0),
    _baseUnit("second")
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ProgressMonitorTime::~ProgressMonitorTime(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set units for simulation time in output.
void
pylith::problems::ProgressMonitorTime::setTimeUnit(const char* value) {
    _baseUnit = value;

    spatialdata::units::Parser parser;
    _baseTime = parser.parse(_baseUnit.c_str());
} // setTimeUnit


// ---------------------------------------------------------------------------------------------------------------------
// Set units for simulation time in output.
const char*
pylith::problems::ProgressMonitorTime::getTimeUnit(void) const {
    return _baseUnit.c_str();
} // getTimeUnit


// ---------------------------------------------------------------------------------------------------------------------
// Open progress monitor.
void
pylith::problems::ProgressMonitorTime::_open(void) {
    _sout.open(getFilename());
    _sout << "Timestamp                     Simulation t   % complete   Est. completion" << std::endl;
} // _open


// ---------------------------------------------------------------------------------------------------------------------
// Close progress monitor.
void
pylith::problems::ProgressMonitorTime::_close(void) {
    if (_sout.is_open()) {
        _sout.close();
    } // if
} // _close


// ---------------------------------------------------------------------------------------------------------------------
// Update progress.
void
pylith::problems::ProgressMonitorTime::_update(const double current,
                                               const time_t& now,
                                               const double percentComplete,
                                               const char* finished) {
    assert(_sout.is_open());
    const double tSimNorm = current / _baseTime;
    std::tm* now_tm = localtime(&now);
    _sout << asctime(now_tm) << "   "
          << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(8) << tSimNorm << "*" << _baseUnit
          << std::setprecision(0) << std::setw(10) << percentComplete
          << "   " << finished
          << std::endl;
} // _update


// End of file
