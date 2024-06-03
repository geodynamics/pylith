// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/problems/ProgressMonitor.hh" // implementation of class methods

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <mpi.h> // USES MPI_Comm_rank, MPI_COMM_WORLD

#include <fstream> // HASA std::ofstream
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::ProgressMonitor::ProgressMonitor(void) :
    _updatePercent(5.0),
    _filename("progress.txt"),
    _iUpdate(-1),
    _isMaster(true) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ProgressMonitor::~ProgressMonitor(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::ProgressMonitor::deallocate(void) {
}


// ------------------------------------------------------------------------------------------------
// Set how often to report status.
void
pylith::problems::ProgressMonitor::setUpdatePercent(const double value) {
    if (value <= 0.0) {
        std::ostringstream msg;
        msg << "Update percentage value (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    } // if

    _updatePercent = value;
} // setUpdatePercent


// ------------------------------------------------------------------------------------------------
// Get how often to report status.
double
pylith::problems::ProgressMonitor::getUpdatePercent(void) const {
    return _updatePercent;
} // getUpdatePercent


// ------------------------------------------------------------------------------------------------
// Set filename for output.
void
pylith::problems::ProgressMonitor::setFilename(const char* filename) {
    if (!strlen(filename)) {
        throw std::runtime_error("Progress monitor output filename set to empty string.");
    } // if
    _filename = filename;
} // setFilename


// ------------------------------------------------------------------------------------------------
// Get filename for output.
const char*
pylith::problems::ProgressMonitor::getFilename(void) const {
    return _filename.c_str();
} // getFilename


// ------------------------------------------------------------------------------------------------
// Open progress monitor.
void
pylith::problems::ProgressMonitor::open(void) {
    _iUpdate = -1;
    _startTime = time(NULL);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    _isMaster = 0 == rank;

    if (_isMaster) {
        _open();
    } // if
} // open


// ------------------------------------------------------------------------------------------------
// Close progress monitor.
void
pylith::problems::ProgressMonitor::close(void) {
    if (_isMaster) {
        _close();
    } // if
} // close


// ------------------------------------------------------------------------------------------------
// Compute finish time.
std::string
pylith::problems::ProgressMonitor::_calcFinishTime(const double percentComplete,
                                                   const time_t& now,
                                                   const time_t& startTime) {
    std::string finished;
    if (percentComplete > 0.0) {
        const double durationSec = 100.0 / percentComplete * difftime(now, startTime);

        struct std::tm start_tm = *localtime(&startTime);
        struct std::tm finished_tm = start_tm;
        finished_tm.tm_sec += durationSec;
        mktime(&finished_tm);
        finished = asctime(&finished_tm);
        finished = finished.erase(finished.find_last_not_of('\n')+1);
    } else {
        finished = "TBD";
    } // if/else

    return finished;
} // calcFinishTime


// End of file
