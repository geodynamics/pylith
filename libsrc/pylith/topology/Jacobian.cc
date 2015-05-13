// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
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

#include "Jacobian.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh
#include "Field.hh" // USES Field

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR
#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Jacobian::Jacobian(const Field& field,
                                     const char* matrixType,
                                     const bool blockOkay) :
  _matrix(0),
  _valuesChanged(true)
{ // constructor
  PYLITH_METHOD_BEGIN;

  PetscDM dmMesh = field.dmMesh();assert(dmMesh);

  const char* msg = "Could not create PETSc sparse matrix associated with system Jacobian.";
  PetscErrorCode err = DMCreateMatrix(dmMesh, &_matrix);PYLITH_CHECK_ERROR_MSG(err, msg);

  _type = matrixType;

  PYLITH_METHOD_END;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Jacobian::~Jacobian(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Jacobian::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = MatDestroy(&_matrix);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Get PETSc matrix.
const PetscMat
pylith::topology::Jacobian::matrix(void) const
{ // matrix
  return _matrix;
} // matrix

// ----------------------------------------------------------------------
// Get PETSc matrix.
PetscMat
pylith::topology::Jacobian::matrix(void)
{ // matrix
  return _matrix;
} // matrix

// ----------------------------------------------------------------------
// Get matrix type.
const char*
pylith::topology::Jacobian::matrixType(void) const
{ // matrixType
  return _type.c_str();
} // matrixType

// ----------------------------------------------------------------------
// Assemble matrix.
void
pylith::topology::Jacobian::assemble(const char* mode)
{ // assemble
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = 0;
  if (0 == strcmp(mode, "final_assembly")) {
    err = MatAssemblyBegin(_matrix, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(_matrix, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);

#if 0 // DEBUGGING
    // Check for empty row
    const PetscInt *cols;
    PetscInt rStart, rEnd, ncols;

    err = MatGetOwnershipRange(_matrix, &rStart, &rEnd);PYLITH_CHECK_ERROR(err);
    for(PetscInt r = rStart; r < rEnd; ++r) {
      PetscInt c;

      err = MatGetRow(_matrix,r, &ncols, &cols, PETSC_NULL);PYLITH_CHECK_ERROR(err);
      if (!ncols) {
        std::ostringstream msg;
        msg << "ERROR: Matrix empty row " << r << " in ["<<rStart<<","<<rEnd<<")" << std::endl;
        throw std::runtime_error(msg.str());
      }
      for(c = 0; c < ncols; ++c) {
        if (cols[c] == r) break;
      }
      if (c == ncols) {
        std::ostringstream msg;
        msg << "ERROR: Matrix row " << r << " in ["<<rStart<<","<<rEnd<<") is missing diagonal element" << std::endl;
        throw std::runtime_error(msg.str());
      }
      err = MatRestoreRow(_matrix,r, &ncols, &cols, PETSC_NULL);PYLITH_CHECK_ERROR(err);
    }
#endif

  } else if (0 == strcmp(mode, "flush_assembly")) {
    err = MatAssemblyBegin(_matrix, MAT_FLUSH_ASSEMBLY);PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(_matrix, MAT_FLUSH_ASSEMBLY);PYLITH_CHECK_ERROR(err);
  } else
    throw std::runtime_error("Unknown mode for assembly of sparse matrix "
			     "associated with system Jacobian.");

  _valuesChanged = true;

  PYLITH_METHOD_END;
} // assemble

// ----------------------------------------------------------------------
// Set entries in matrix to zero (retain structure).
void
pylith::topology::Jacobian::zero(void)
{ // zero
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = MatZeroEntries(_matrix);PYLITH_CHECK_ERROR(err);
  _valuesChanged = true;

  PYLITH_METHOD_END;
} // zero

// ----------------------------------------------------------------------
// View matrix to stdout.
void
pylith::topology::Jacobian::view(void) const
{ // view
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = MatView(_matrix, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // view

// ----------------------------------------------------------------------
// Write matrix to binary file.
void
pylith::topology::Jacobian::write(const char* filename,
                                  const MPI_Comm comm)
{ // write
  PYLITH_METHOD_BEGIN;

  PetscViewer viewer;
  PetscErrorCode err = PetscViewerBinaryOpen(comm, filename, FILE_MODE_WRITE, &viewer);PYLITH_CHECK_ERROR(err);

  err = MatView(_matrix, viewer); PYLITH_CHECK_ERROR(err);
  err = PetscViewerDestroy(&viewer); PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // write

// ----------------------------------------------------------------------
// Verify symmetry of matrix.
void
pylith::topology::Jacobian::verifySymmetry(void) const
{ // verifySymmetry
  PYLITH_METHOD_BEGIN;

  const PetscMat matSparse = _matrix;
  PetscErrorCode err;

  int nrows = 0;
  int ncols = 0;
  err = MatGetSize(matSparse, &nrows, &ncols);PYLITH_CHECK_ERROR(err);

  PetscMat matDense;
  PetscMat matSparseAIJ;
  err = MatConvert(matSparse, MATSEQAIJ, MAT_INITIAL_MATRIX, &matSparseAIJ);PYLITH_CHECK_ERROR(err);
  err = MatConvert(matSparseAIJ, MATSEQDENSE, MAT_INITIAL_MATRIX, &matDense);PYLITH_CHECK_ERROR(err);

  scalar_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  err = MatGetValues(matDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);PYLITH_CHECK_ERROR(err);
  const PylithScalar tolerance = 1.0e-06;
  bool isSymmetric = true;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int indexIJ = ncols*iRow+iCol;
      const int indexJI = nrows*iCol+iRow;
      const PylithScalar valIJ = vals[indexIJ];
      const PylithScalar valJI = vals[indexJI];
      if (fabs(valIJ) > 1.0) {
        if (fabs(1.0 - valJI/valIJ) > tolerance) {
          std::cerr << "Mismatch: " 
                    << "(" << iRow << ", " << iCol << ") = " << valIJ
                    << ", (" << iCol << ", " << iRow << ") = " << valJI
                    << std::endl;
          isSymmetric = false;
        } // if
      } else {
        if (fabs(valJI - valIJ) > tolerance) {
          std::cerr << "Mismatch: " 
                    << "(" << iRow << ", " << iCol << ") = " << valIJ
                    << ", (" << iCol << ", " << iRow << ") = " << valJI
                    << std::endl;
          isSymmetric = false;
        } // if
      } // if/else
    } // for
  err = MatDestroy(&matDense);PYLITH_CHECK_ERROR(err);
  err = MatDestroy(&matSparseAIJ);PYLITH_CHECK_ERROR(err);
  if (!isSymmetric)
    throw std::runtime_error("Jacobian matrix is not symmetric.");

  PYLITH_METHOD_END;
} // verifySymmetry

// ----------------------------------------------------------------------
// Get flag indicating if sparse matrix values have been
// updated.
bool
pylith::topology::Jacobian::valuesChanged(void) const
{ // valuesChanged
  return _valuesChanged;
} // valuesChanged

// ----------------------------------------------------------------------
// Reset flag indicating if sparse matrix values have been updated.
void
pylith::topology::Jacobian::resetValuesChanged(void)
{ // resetValuesChanged
  _valuesChanged = false;
} // resteValuesChanged


// End of file 
