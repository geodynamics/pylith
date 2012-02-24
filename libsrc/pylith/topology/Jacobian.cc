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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Jacobian.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh
#include "SubMesh.hh" // USES SubMesh
#include "Field.hh" // USES Field

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Jacobian::Jacobian(const Field<Mesh>& field,
                                     const char* matrixType,
                                     const bool blockOkay) :
  _matrix(0),
  _valuesChanged(true)
{ // constructor
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = field.mesh().sieveMesh();
  const ALE::Obj<Mesh::RealSection>& fieldSection = field.section();
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Jacobian");

  // Set blockFlag to -1 if okay to set block size equal to fiber
  // dimension, otherwise use a block size of 1.
  const int blockFlag = (blockOkay) ? -1 : 1;

  PetscErrorCode err = DMMeshCreateMatrix(sieveMesh, fieldSection,
					matrixType, &_matrix, blockFlag);
  CHECK_PETSC_ERROR_MSG(err, "Could not create PETSc sparse matrix "
			"associated with system Jacobian.");
  logger.stagePop();

  _type = matrixType;
} // constructor

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Jacobian::Jacobian(const Field<SubMesh>& field,
                                     const char* matrixType,
                                     const bool blockOkay) :
  _matrix(0),
  _valuesChanged(true)
{ // constructor
  const ALE::Obj<SubMesh::SieveMesh>& sieveMesh = field.mesh().sieveMesh();
  const ALE::Obj<SubMesh::RealSection>& fieldSection = field.section();
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Jacobian");

  // Set blockFlag to -1 if okay to set block size equal to fiber
  // dimension, otherwise use a block size of 1.
  const int blockFlag = (blockOkay) ? -1 : 1;

  PetscErrorCode err = DMMeshCreateMatrix(sieveMesh, fieldSection,
					  matrixType, &_matrix, blockFlag);
  CHECK_PETSC_ERROR_MSG(err, "Could not create PETSc sparse matrix "
      "associated with subsystem Jacobian.");
  logger.stagePop();

  _type = matrixType;
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
  PetscErrorCode err = MatDestroy(&_matrix);CHECK_PETSC_ERROR(err);
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
  PetscErrorCode err = 0;
  if (0 == strcmp(mode, "final_assembly")) {
    err = MatAssemblyBegin(_matrix, MAT_FINAL_ASSEMBLY); 
    CHECK_PETSC_ERROR(err);
    err = MatAssemblyEnd(_matrix, MAT_FINAL_ASSEMBLY);
    CHECK_PETSC_ERROR(err);
  } else if (0 == strcmp(mode, "flush_assembly")) {
    err = MatAssemblyBegin(_matrix, MAT_FLUSH_ASSEMBLY);
    CHECK_PETSC_ERROR(err);
    err = MatAssemblyEnd(_matrix, MAT_FLUSH_ASSEMBLY);
    CHECK_PETSC_ERROR(err);
  } else
    throw std::runtime_error("Unknown mode for assembly of sparse matrix "
			     "associated with system Jacobian.");

  _valuesChanged = true;
} // assemble

// ----------------------------------------------------------------------
// Set entries in matrix to zero (retain structure).
void
pylith::topology::Jacobian::zero(void)
{ // zero
  PetscErrorCode err = MatZeroEntries(_matrix);
  CHECK_PETSC_ERROR(err);
  _valuesChanged = true;
} // zero

// ----------------------------------------------------------------------
// View matrix to stdout.
void
pylith::topology::Jacobian::view(void) const
{ // view
  PetscErrorCode err = MatView(_matrix, PETSC_VIEWER_STDOUT_WORLD);
  CHECK_PETSC_ERROR(err);
} // view

// ----------------------------------------------------------------------
// Write matrix to binary file.
void
pylith::topology::Jacobian::write(const char* filename,
                                  const MPI_Comm comm)
{ // write
  PetscViewer viewer;

  PetscErrorCode err = 
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_WRITE, &viewer);
  CHECK_PETSC_ERROR(err);

  err = MatView(_matrix, viewer); CHECK_PETSC_ERROR(err);
  err = PetscViewerDestroy(&viewer); CHECK_PETSC_ERROR(err);
} // write

// ----------------------------------------------------------------------
// Verify symmetry of matrix.
void
pylith::topology::Jacobian::verifySymmetry(void) const
{ // verifySymmetry
  const PetscMat matSparse = _matrix;

  int nrows = 0;
  int ncols = 0;
  MatGetSize(matSparse, &nrows, &ncols);

  PetscMat matDense;
  PetscMat matSparseAIJ;
  MatConvert(matSparse, MATSEQAIJ, MAT_INITIAL_MATRIX, &matSparseAIJ);
  MatConvert(matSparseAIJ, MATSEQDENSE, MAT_INITIAL_MATRIX, &matDense);

  double_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatGetValues(matDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);
  const double tolerance = 1.0e-06;
  bool isSymmetric = true;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int indexIJ = ncols*iRow+iCol;
      const int indexJI = nrows*iCol+iRow;
      const double valIJ = vals[indexIJ];
      const double valJI = vals[indexJI];
      if (fabs(valIJ) > 1.0)
	if (fabs(1.0 - valJI/valIJ) > tolerance) {
	  std::cerr << "Mismatch: " 
		    << "(" << iRow << ", " << iCol << ") = " << valIJ
		    << ", (" << iCol << ", " << iRow << ") = " << valJI
		    << std::endl;
	  isSymmetric = false;
	} // if
      else
	if (fabs(valJI - valIJ) > tolerance) {
	  std::cerr << "Mismatch: " 
		    << "(" << iRow << ", " << iCol << ") = " << valIJ
		    << ", (" << iCol << ", " << iRow << ") = " << valJI
		    << std::endl;
	  isSymmetric = false;
	} // if
    } // for
  MatDestroy(&matDense);
  MatDestroy(&matSparseAIJ);
  if (!isSymmetric)
    throw std::runtime_error("Jacobian matrix is not symmetric.");
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
