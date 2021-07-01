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

#include "DataWriterHDF5.hh" // Implementation of class methods

#include "HDF5.hh" // USES HDF5
#include "Xdmf.hh" // USES Xdmf

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/MeshOps.hh" // USES isCohesiveCell()
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "petscviewerhdf5.h"
#include <mpi.h> // USES MPI routines

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

extern "C" {
    extern PetscErrorCode VecView_Seq(Vec,
                                      PetscViewer);

    extern PetscErrorCode VecView_MPI(Vec,
                                      PetscViewer);

}

#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8
#define PYLITH_HDF5_USE_API_18
#endif

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterHDF5::DataWriterHDF5(void) :
    _filename("output.h5"),
    _viewer(0),
    _tstamp(0),
    _tstampIndex(0) {
    PyreComponent::setName("datawriterhdf5");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterHDF5::~DataWriterHDF5(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriterHDF5::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    DataWriter::deallocate();

    PetscErrorCode err = 0;
    err = PetscViewerDestroy(&_viewer);PYLITH_CHECK_ERROR(err);assert(!_viewer);
    err = VecDestroy(&_tstamp);PYLITH_CHECK_ERROR(err);assert(!_tstamp);

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriterHDF5::DataWriterHDF5(const DataWriterHDF5& w) :
    DataWriter(w),
    _filename(w._filename),
    _viewer(0),
    _tstamp(0),
    _tstampIndex(0) {}


// ---------------------------------------------------------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriterHDF5::open(const pylith::topology::Mesh& mesh,
                                     const bool isInfo) {
    PYLITH_METHOD_BEGIN;

    DataWriter::open(mesh, isInfo);

    try {
        PetscErrorCode err = 0;

        deallocate();

        const std::string& filename = hdf5Filename();

        _timesteps.clear();
        _tstampIndex = 0;
        PetscMPIInt commRank;
        err = MPI_Comm_rank(mesh.getComm(), &commRank);PYLITH_CHECK_ERROR(err);
        const int localSize = (!commRank) ? 1 : 0;
        err = VecCreateMPI(mesh.getComm(), localSize, 1, &_tstamp);PYLITH_CHECK_ERROR(err);assert(_tstamp);
        err = VecSetBlockSize(_tstamp, 1);PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject) _tstamp, "time");PYLITH_CHECK_ERROR(err);

        err = PetscViewerHDF5Open(mesh.getComm(), filename.c_str(), FILE_MODE_WRITE, &_viewer);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5SetBaseDimension2(_viewer, PETSC_TRUE);PYLITH_CHECK_ERROR(err);

        err = PetscViewerHDF5PushGroup(_viewer, "/geometry");PYLITH_CHECK_ERROR(err);
        PetscVec coordsGlobalVec = NULL;
        DataWriter::getCoordsGlobalVec(&coordsGlobalVec, mesh);
        PetscBool isseq;
        err = PetscObjectTypeCompare((PetscObject) coordsGlobalVec, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
        if (isseq) {
            err = VecView_Seq(coordsGlobalVec, _viewer);PYLITH_CHECK_ERROR(err);
        } else {
            err = VecView_MPI(coordsGlobalVec, _viewer);PYLITH_CHECK_ERROR(err);
        } // if/else
        err = VecDestroy(&coordsGlobalVec);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);

        PetscInt vStart, vEnd, cellHeight, cStart, cEnd, conesSize = 0, numCorners, numCornersLocal = -1;

        PetscDM dmMesh = mesh.getDM();assert(dmMesh);
        err = DMPlexGetVTKCellHeight(dmMesh, &cellHeight);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
        for (PetscInt cell = cStart; cell < cEnd; ++cell) {
            PetscInt *closure = NULL;
            PetscInt closureSize, v;

            if (pylith::topology::MeshOps::isCohesiveCell(dmMesh, cell)) { continue; }
            err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            PetscInt numCornersCell = 0;
            for (v = 0; v < closureSize*2; v += 2) {
                if ((closure[v] >= vStart) && (closure[v] < vEnd)) {
                    ++numCornersCell;
                } // if
            } // for
            if (-1 == numCornersLocal) {
                numCornersLocal = numCornersCell;
            } else {
                assert(numCornersCell == numCornersLocal); // All cells in output must have the same number of corners.
            } // if/else
            conesSize += numCornersCell;
            err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        } // for
        err = MPI_Allreduce(&numCornersLocal, &numCorners, 1, MPIU_INT, MPI_MAX, mesh.getComm());PYLITH_CHECK_ERROR(err);

        PetscIS globalVertexNumbers = NULL;
        const PetscInt *gvertex = NULL;
        PetscVec cellVec = NULL;
        PetscScalar *vertices = NULL;

        err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
        err = VecCreate(mesh.getComm(), &cellVec);PYLITH_CHECK_ERROR(err);
        err = VecSetSizes(cellVec, conesSize, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
        err = VecSetBlockSize(cellVec, numCorners);PYLITH_CHECK_ERROR(err);
        err = VecSetFromOptions(cellVec);PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject) cellVec, "cells");PYLITH_CHECK_ERROR(err);
        err = VecGetArray(cellVec, &vertices);PYLITH_CHECK_ERROR(err);
        for (PetscInt cell = cStart, v = 0; cell < cEnd; ++cell) {
            PetscInt *closure = NULL;
            PetscInt closureSize, nC = 0, p;

            if (pylith::topology::MeshOps::isCohesiveCell(dmMesh, cell)) { continue; }
            err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            for (p = 0; p < closureSize*2; p += 2) {
                if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
                    closure[nC++] = closure[p];
                } // if
            } // for
            DMPolytopeType ct;
            err = DMPlexGetCellType(dmMesh, cell, &ct);PYLITH_CHECK_ERROR(err);
            err = DMPlexInvertCell(ct, closure);PYLITH_CHECK_ERROR(err);
            for (p = 0; p < nC; ++p) {
                const PetscInt gv = gvertex[closure[p] - vStart];
                vertices[v++] = gv < 0 ? -(gv+1) : gv;
            }
            err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            // assert(v == (cell-cStart+1)*numCorners); Would be true without the label check
        } // for
        err = VecRestoreArray(cellVec, &vertices);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5PushGroup(_viewer, "/topology");PYLITH_CHECK_ERROR(err);
        err = VecView(cellVec, _viewer);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);
        err = VecDestroy(&cellVec);PYLITH_CHECK_ERROR(err);
        err = ISRestoreIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);

        hid_t h5 = -1;
        err = PetscViewerHDF5GetFileId(_viewer, &h5);PYLITH_CHECK_ERROR(err);
        assert(h5 >= 0);
        const int cellDim = mesh.getDimension();
        HDF5::writeAttribute(h5, "/topology/cells", "cell_dim", (void*)&cellDim, H5T_NATIVE_INT);
    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error while opening HDF5 file " << hdf5Filename() << ".\n" << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Unknown error while opening HDF5 file " << hdf5Filename() << ".";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // open


// ---------------------------------------------------------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriterHDF5::close(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;
    err = PetscViewerDestroy(&_viewer);PYLITH_CHECK_ERROR(err);assert(!_viewer);
    err = VecDestroy(&_tstamp);PYLITH_CHECK_ERROR(err);assert(!_tstamp);

    _timesteps.clear();
    _tstampIndex = 0;

    if (isOpen()) {
        // Write Xdmf file on process 0
        PetscMPIInt commRank;
        err = MPI_Comm_rank(PETSC_COMM_WORLD, &commRank);PYLITH_CHECK_ERROR(err);
        if (!commRank) {
            try {
                Xdmf::write(hdf5Filename().c_str());
            } catch (const std::exception& err) {
                pythia::journal::error_t error("datawriter");
                error << err.what() << pythia::journal::endl;
            } // catch
        } // if
    } // if

    DataWriter::close();

    PYLITH_METHOD_END;
} // close


// ---------------------------------------------------------------------------------------------------------------------
// Write field over vertices to file.
void
pylith::meshio::DataWriterHDF5::writeVertexField(const PylithScalar t,
                                                 const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;
    assert(_viewer);

    const char* name = subfield.getDescription().label.c_str();
    try {
        PetscErrorCode err;

        if (_timesteps.find(name) == _timesteps.end()) {
            _timesteps[name] = 0;
        } else {
            _timesteps[name] += 1;
        }
        const int istep = _timesteps[name];
        // Add time stamp to "/time" if necessary.
        MPI_Comm comm;
        err = PetscObjectGetComm((PetscObject)subfield.getDM(), &comm);PYLITH_CHECK_ERROR(err);
        PetscMPIInt commRank;
        err = MPI_Comm_rank(comm, &commRank);PYLITH_CHECK_ERROR(err);
        if (_tstampIndex == istep) {
            _writeTimeStamp(t, commRank);
        } // if

        err = PetscViewerHDF5PushGroup(_viewer, "/vertex_fields");PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5PushTimestepping(_viewer);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5SetTimestep(_viewer, istep);PYLITH_CHECK_ERROR(err);

        PetscVec vector = subfield.getVector();assert(vector);
        PetscBool isseq;
        err = PetscObjectTypeCompare((PetscObject) vector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
        if (isseq) {
            err = VecView_Seq(vector, _viewer);PYLITH_CHECK_ERROR(err);
        } else {
            err = VecView_MPI(vector, _viewer);PYLITH_CHECK_ERROR(err);
        }
        err = PetscViewerHDF5PopTimestepping(_viewer);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);

        if (0 == istep) {
            hid_t h5 = -1;
            err = PetscViewerHDF5GetFileId(_viewer, &h5);PYLITH_CHECK_ERROR(err);
            assert(h5 >= 0);
            std::string fullName = std::string("/vertex_fields/") + std::string(name);
            const char* sattr = pylith::topology::FieldBase::vectorFieldString(subfield.getDescription().vectorFieldType);
            HDF5::writeAttribute(h5, fullName.c_str(), "vector_field_type", sattr);
        } // if

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error while writing field '" << name << "' at time "
            << t << " to HDF5 file '" << hdf5Filename() << "'.\n" << err.what();
        throw std::runtime_error(msg.str());

    } catch (...) {
        std::ostringstream msg;
        msg << "Error while writing field '" << name << "' at time "
            << t << " to HDF5 file '" << hdf5Filename() << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writeVertexField


// ---------------------------------------------------------------------------------------------------------------------
// Write field over cells to file.
void
pylith::meshio::DataWriterHDF5::writeCellField(const PylithScalar t,
                                               const pylith::meshio::OutputSubfield& subfield) {
    PYLITH_METHOD_BEGIN;

    assert(_viewer);

    const char* name = subfield.getDescription().label.c_str();
    try {
        PetscErrorCode err;

        if (_timesteps.find(name) == _timesteps.end()) {
            _timesteps[name] = 0;
        } else {
            _timesteps[name] += 1;
        }
        const int istep = _timesteps[name];
        // Add time stamp to "/time" if necessary.
        MPI_Comm comm;
        err = PetscObjectGetComm((PetscObject)subfield.getDM(), &comm);PYLITH_CHECK_ERROR(err);
        PetscMPIInt commRank;
        err = MPI_Comm_rank(comm, &commRank);PYLITH_CHECK_ERROR(err);
        if (_tstampIndex == istep) {
            _writeTimeStamp(t, commRank);
        } // if

        err = PetscViewerHDF5PushGroup(_viewer, "/cell_fields");PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5PushTimestepping(_viewer);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5SetTimestep(_viewer, istep);PYLITH_CHECK_ERROR(err);

        PetscVec vector = subfield.getVector();assert(vector);
        PetscBool isseq;
        err = PetscObjectTypeCompare((PetscObject) vector, VECSEQ, &isseq);PYLITH_CHECK_ERROR(err);
        if (isseq) {
            err = VecView_Seq(vector, _viewer);PYLITH_CHECK_ERROR(err);
        } else {
            err = VecView_MPI(vector, _viewer);PYLITH_CHECK_ERROR(err);
        } // if/else
        err = PetscViewerHDF5PopTimestepping(_viewer);PYLITH_CHECK_ERROR(err);
        err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);

        if (0 == istep) {
            hid_t h5 = -1;
            err = PetscViewerHDF5GetFileId(_viewer, &h5);PYLITH_CHECK_ERROR(err);
            assert(h5 >= 0);
            std::string fullName = std::string("/cell_fields/") + std::string(name);
            const char* sattr = pylith::topology::FieldBase::vectorFieldString(subfield.getDescription().vectorFieldType);
            HDF5::writeAttribute(h5, fullName.c_str(), "vector_field_type", sattr);
        } // if
    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error while writing field '" << name << "' at time "
            << t << " to HDF5 file '" << hdf5Filename() << "'.\n" << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        std::ostringstream msg;
        msg << "Error while writing field '" << name << "' at time "
            << t << " to HDF5 file '" << hdf5Filename() << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writeCellField


// ---------------------------------------------------------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::DataWriterHDF5::writePointNames(const pylith::string_vector& names,
                                                const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    assert(_viewer);

    char* namesFixedLength = NULL;
    try {
        // Put station names into array of fixed length strings
        // (numNames*maxStringLegnth) on each process, and then write
        // collectively.
        int mpierr;
        MPI_Comm comm = mesh.getComm();
        const int commRank = mesh.getCommRank();
        int nprocs = 0;
        mpierr = MPI_Comm_size(comm, &nprocs);assert(MPI_SUCCESS == mpierr);

        // Number of names on each process.
        const int numNamesLocal = names.size();
        int_array numNamesArray(nprocs);
        // Use void* for compatibility with OpenMPI 1.3 on Travis-CI
        mpierr = MPI_Allgather((void*)&numNamesLocal, 1, MPI_INT, &numNamesArray[0], 1, MPI_INT, comm);assert(MPI_SUCCESS == mpierr);
        const int numNames = numNamesArray.sum();

        // Get maximum string length.
        int maxStringLengthLocal = 0;
        int maxStringLength = 0;
        for (int i = 0; i < numNamesLocal; ++i) {
            maxStringLengthLocal = std::max(maxStringLengthLocal, int(names[i].length()));
        } // for
        maxStringLengthLocal += 1; // add space for null terminator.
        mpierr = MPI_Allreduce(&maxStringLengthLocal, &maxStringLength, 1, MPI_INT, MPI_MAX, comm);assert(MPI_SUCCESS == mpierr);

        namesFixedLength = (numNamesLocal*maxStringLength > 0) ? new char[numNamesLocal*maxStringLength] : NULL;
        for (int i = 0; i < numNamesLocal; ++i) {
            const int index = i*maxStringLength;
            strncpy(&namesFixedLength[index], names[i].c_str(), maxStringLength-1);
            namesFixedLength[index+maxStringLength-1] = '\0';
            // Fill remaining portion of string with null characters.
            for (int j = names[i].length(); j < maxStringLength; ++j) {
                namesFixedLength[index+j] = '\0';
            } // for
        } // for

        hid_t h5 = -1;
        PetscErrorCode petscerr = PetscViewerHDF5GetFileId(_viewer, &h5);PYLITH_CHECK_ERROR(petscerr);
        assert(h5 >= 0);
        const char* parent = "/";
        const char* name = "stations";

        // Open group
#if defined(PYLITH_HDF5_USE_API_18)
        hid_t group = H5Gopen2(h5, parent, H5P_DEFAULT);
#else
        hid_t group = H5Gopen(h5, parent);
#endif
        if (group < 0) { throw std::runtime_error("Could not open group.");}

        hid_t datatype = H5Tcopy(H5T_C_S1);
        if (datatype < 0) { throw std::runtime_error("Could not create datatype.");}
        herr_t err = H5Tset_size(datatype, maxStringLength);
        if (err < 0) { throw std::runtime_error("Could not set size of datatype.");}

        // Create the filespace
        const int ndims = 1;
        hsize_t dims[ndims];
        dims[0] = numNames;
        hid_t filespace = H5Screate_simple(ndims, dims, NULL);
        if (filespace < 0) { throw std::runtime_error("Could not create filespace.");}

#if defined(PYLITH_HDF5_USE_API_18)
        hid_t dataset = H5Dcreate2(group, name, datatype, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
        hid_t dataset = H5Dcreate(group, name, datatype, filespace, H5P_DEFAULT);
#endif
        if (dataset < 0) { throw std::runtime_error("Could not create dataset.");}
        err = H5Sclose(filespace);
        if (err < 0) { throw std::runtime_error("Could not close filespace.");}

        // Create the memspace
        dims[0] = numNamesLocal;
        hid_t memspace = H5Screate_simple(ndims, dims, NULL);
        if (memspace < 0) {throw std::runtime_error("Could not create memspace.");}

        hid_t dataspace = H5Dget_space(dataset);
        if (dataspace < 0) {throw std::runtime_error("Could not get dataspace.");}
        hsize_t offset[1] = {0};
        for (int i = 0; i < commRank; ++i) {
            offset[0] += numNamesArray[i];
        } // for
        hsize_t count[1];
        count[0] = numNamesLocal;
        err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
        if (err < 0) {throw std::runtime_error("Could not select hyperslab.");}

        hid_t property = H5Pcreate(H5P_DATASET_XFER);
        if (property < 0) {
            throw std::runtime_error("Could not create property.");
        }
        H5Pset_dxpl_mpio(property, H5FD_MPIO_COLLECTIVE);

        err = H5Dwrite(dataset, datatype, memspace, dataspace, property, namesFixedLength);
        if (err < 0) {throw std::runtime_error("Could not write dataset.");}

        err = H5Sclose(memspace);
        if (err < 0) {throw std::runtime_error("Could not close memspace.");}
        err = H5Sclose(dataspace);
        if (err < 0) {throw std::runtime_error("Could not close dataspace.");}
        err = H5Pclose(property);
        if (err < 0) {throw std::runtime_error("Could not close property.");}
        err = H5Dclose(dataset);
        if (err < 0) {throw std::runtime_error("Could not close dataset.");}
        err = H5Tclose(datatype);
        if (err < 0) {throw std::runtime_error("Could not close datatype.");}

        err = H5Gclose(group);
        if (err < 0) {throw std::runtime_error("Could not close group.");}

        delete[] namesFixedLength;namesFixedLength = NULL;
    } catch (const std::exception& err) {
        delete[] namesFixedLength;namesFixedLength = NULL;

        std::ostringstream msg;
        msg << "Error while writing stations to HDF5 file '" << hdf5Filename() << "'.\n" << err.what();
        throw std::runtime_error(msg.str());
    } catch (...) {
        delete[] namesFixedLength;namesFixedLength = NULL;

        std::ostringstream msg;
        msg << "Error while writing stations to HDF5 file '" << hdf5Filename() << "'.";
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // writePointNames


// ---------------------------------------------------------------------------------------------------------------------
// Generate filename for HDF5 file.
std::string
pylith::meshio::DataWriterHDF5::hdf5Filename(void) const {
    PYLITH_METHOD_BEGIN;

    std::ostringstream filename;
    const int indexExt = _filename.find(".h5");
    if (DataWriter::_isInfo) {
        filename << std::string(_filename, 0, indexExt) << "_info.h5";
    } else {
        filename << _filename;
    } // if/else

    PYLITH_METHOD_RETURN(std::string(filename.str()));
} // hdf5Filename


// ---------------------------------------------------------------------------------------------------------------------
// Write time stamp to file.
void
pylith::meshio::DataWriterHDF5::_writeTimeStamp(const PylithScalar t,
                                                const int commRank) {
    assert(_tstamp);
    PetscErrorCode err = 0;

    if (!commRank) {
        const PylithScalar tDim = t * DataWriter::_timeScale;
        err = VecSetValue(_tstamp, 0, tDim, INSERT_VALUES);PYLITH_CHECK_ERROR(err);
    } // if
    err = VecAssemblyBegin(_tstamp);PYLITH_CHECK_ERROR(err);
    err = VecAssemblyEnd(_tstamp);PYLITH_CHECK_ERROR(err);

    err = PetscViewerHDF5PushGroup(_viewer, "/");PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5PushTimestepping(_viewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5SetTimestep(_viewer, _tstampIndex);PYLITH_CHECK_ERROR(err);
    err = VecView(_tstamp, _viewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5PopTimestepping(_viewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerHDF5PopGroup(_viewer);PYLITH_CHECK_ERROR(err);

    _tstampIndex++;
} // _writeTimeStamp


// End of file
