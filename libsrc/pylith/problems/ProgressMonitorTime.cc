// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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
// Deallocate PETSc and local data structures.
void
pylith::problems::ProgressMonitorTime::deallocate(void) {
    ProgressMonitor::deallocate();

    if (_sout.is_open()) {
        _sout.close();
    } // if
} // deallocate


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
    _sout.setf(std::ios::fixed);
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
    std::string now_str = asctime(now_tm);
    now_str = now_str.erase(now_str.find_last_not_of('\n')+1);
    _sout << now_str << "   "
          << std::setprecision(2) << std::setw(8) << tSimNorm
          << "*" << std::left << std::setw(6) << _baseUnit << std::right
          << std::setprecision(0) << std::setw(13) << percentComplete
          << "   " << finished
          << std::endl;
} // _update


// End of file
