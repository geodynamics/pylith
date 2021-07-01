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

#include "DataWriterVTK.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES StratumIS
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfieldIS

#include "petscdmplex.h"

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

extern
PetscErrorCode DMPlexVTKWriteAll(PetscObject odm,
                                 PetscViewer viewer);

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterVTK::DataWriterVTK(void) :
    _timeConstant(1.0),
    _precision(6),
    _filename("output.vtk"),
    _timeFormat("%f"),
    _viewer(NULL),
    _dm(NULL),
    _isOpenTimeStep(false),
    _wroteVertexHeader(false),
    _wroteCellHeader(false) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterVTK::~DataWriterVTK(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriterVTK::deallocate(void) { // deallocate
    PYLITH_METHOD_BEGIN;

    closeTimeStep(); // Insure time step is closed.
    close(); // Insure clean up.
    DataWriter::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriterVTK::DataWriterVTK(const DataWriterVTK& w) :
    DataWriter(w),
    _timeConstant(w._timeConstant),
    _precision(w._precision),
    _filename(w._filename),
    _timeFormat(w._timeFormat),
    _viewer(NULL),
    _dm(NULL),
    _isOpenTimeStep(w._isOpenTimeStep),
    _wroteVertexHeader(w._wroteVertexHeader),
    _wroteCellHeader(w._wroteCellHeader) { // copy constructor
} // copy constructor


// ------------------------------------------------------------------------------------------------
// Set value used to normalize time stamp in name of VTK file.
void
pylith::meshio::DataWriterVTK::timeConstant(const PylithScalar value) {
    PYLITH_METHOD_BEGIN;

    if (value <= 0.0) {
        std::ostringstream msg;
        msg << "Time used to normalize time stamp in VTK data files must be "
            << "positive.\nCurrent value is " << value << ".";
        throw std::runtime_error(msg.str());
    } // if
    _timeConstant = value;

    PYLITH_METHOD_END;
} // timeConstant


// ------------------------------------------------------------------------------------------------
// Set precision of floating point values in output.
void
pylith::meshio::DataWriterVTK::precision(const int value) {
    PYLITH_METHOD_BEGIN;

    if (value <= 0) {
        std::ostringstream msg;
        msg << "Floating point precision (" << value << ") must be positive.";
        throw std::runtime_error(msg.str());
    } // if

    _precision = value;

    PYLITH_METHOD_END;
} // precision


// ------------------------------------------------------------------------------------------------
// Prepare for writing files.
void
pylith::meshio::DataWriterVTK::open(const pylith::topology::Mesh& mesh,
                                    const bool isInfo) {
    PYLITH_METHOD_BEGIN;

    DataWriter::open(mesh, isInfo);

    // Save handle for actions required in closeTimeStep() and close();
    PetscErrorCode err = 0;
    err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
    _dm = mesh.getDM();assert(_dm);
    err = PetscObjectReference((PetscObject) _dm);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // open


// ------------------------------------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriterVTK::close(void) {
    PYLITH_METHOD_BEGIN;

    if (_isOpen) {
        assert(_dm);
        PetscErrorCode err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
    } // if

    DataWriter::close();

    PYLITH_METHOD_END;
} // close


// ------------------------------------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriterVTK::openTimeStep(const PylithScalar t,
                                            const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    assert(_dm && _dm == mesh.getDM());
    assert(_isOpen && !_isOpenTimeStep);

    const std::string& filename = _vtkFilename(t);
    PetscErrorCode err = PetscViewerVTKOpen(mesh.getComm(), filename.c_str(), FILE_MODE_WRITE,
                                            &_viewer);PYLITH_CHECK_ERROR(err);

    _isOpenTimeStep = true;

    PYLITH_METHOD_END;
} // openTimeStep


// ------------------------------------------------------------------------------------------------
/// Cleanup after writing data for a time step.
void
pylith::meshio::DataWriterVTK::closeTimeStep(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;

    // Destroy the viewer (which also writes the file).
    err = PetscViewerDestroy(&_viewer);PYLITH_CHECK_ERROR(err);

    // Remove label
    if (_isOpenTimeStep) {
        assert(_dm);
        PetscBool hasLabel = PETSC_FALSE;
        err = DMHasLabel(_dm, "vtk", &hasLabel);PYLITH_CHECK_ERROR(err);
        if (hasLabel) {
            err = DMClearLabelStratum(_dm, "vtk", 1);PYLITH_CHECK_ERROR(err);
            err = DMClearLabelStratum(_dm, "vtk", 2);PYLITH_CHECK_ERROR(err);
        } // if
    } // if

    _isOpenTimeStep = false;
    _wroteVertexHeader = false;
    _wroteCellHeader = false;

    PYLITH_METHOD_END;
} // closeTimeStep


// ------------------------------------------------------------------------------------------------
// Write field over vertices to file.
void
pylith::meshio::DataWriterVTK::writeVertexField(const PylithScalar t,
                                                const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;
    assert(_isOpen && _isOpenTimeStep);

    PetscErrorCode err;
    PetscViewerVTKFieldType ft = subfield.getDescription().vectorFieldType == pylith::topology::FieldBase::VECTOR ?
                                 PETSC_VTK_POINT_VECTOR_FIELD : PETSC_VTK_POINT_FIELD;
    err = PetscViewerVTKAddField(_viewer, (PetscObject)_dm, DMPlexVTKWriteAll, PETSC_DEFAULT, ft,
                                 PETSC_TRUE, (PetscObject)subfield.getVector());PYLITH_CHECK_ERROR(err);
    err = PetscObjectReference((PetscObject) subfield.getVector());PYLITH_CHECK_ERROR(err); // Viewer destroys Vec

    _wroteVertexHeader = true;

    PYLITH_METHOD_END;
} // writeVertexField


// ------------------------------------------------------------------------------------------------
// Write field over cells to file.
void
pylith::meshio::DataWriterVTK::writeCellField(const PylithScalar t,
                                              const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;

    assert(_isOpen && _isOpenTimeStep);

    PetscErrorCode err;
    PetscViewerVTKFieldType ft = subfield.getDescription().vectorFieldType == pylith::topology::FieldBase::VECTOR ?
                                 PETSC_VTK_CELL_VECTOR_FIELD : PETSC_VTK_CELL_FIELD;
    err = PetscViewerVTKAddField(_viewer, (PetscObject)_dm, DMPlexVTKWriteAll, PETSC_DEFAULT, ft,
                                 PETSC_TRUE, (PetscObject)subfield.getVector());PYLITH_CHECK_ERROR(err);
    err = PetscObjectReference((PetscObject) subfield.getVector());PYLITH_CHECK_ERROR(err); // Viewer destroys Vec

    _wroteCellHeader = true;

    PYLITH_METHOD_END;
} // writeCellField


// ------------------------------------------------------------------------------------------------
// Generate filename for VTK file.
std::string
pylith::meshio::DataWriterVTK::_vtkFilename(const PylithScalar t) const {
    PYLITH_METHOD_BEGIN;

    std::ostringstream filename;
    const int indexExt = _filename.find(".vtk");
    if (!DataWriter::_isInfo) {
        // If data with multiple time steps, then add time stamp to filename
        char sbuffer[256];
        sprintf(sbuffer, _timeFormat.c_str(), t * _timeScale / _timeConstant);
        std::string timestamp(sbuffer);
        const size_t pos = timestamp.find(".");
        if (pos != std::string::npos) {
            timestamp.erase(pos, 1);
        } // if
        filename << std::string(_filename, 0, indexExt) << "_t" << timestamp << ".vtk";
    } else {
        filename << std::string(_filename, 0, indexExt) << "_info.vtk";
    } // if/else

    PYLITH_METHOD_RETURN(std::string(filename.str()));
} // _vtkFilename


// End of file
