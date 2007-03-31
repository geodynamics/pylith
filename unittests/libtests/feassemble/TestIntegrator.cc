// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIntegrator.hh" // Implementation of class methods

#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D

#include <petscmat.h>

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegrator );

// ----------------------------------------------------------------------
// Test clone().
void
pylith::feassemble::TestIntegrator::testCopy(void)
{ // testCopy
  // Test copy constructor by testing value of minJacobian value in quadrature

  Quadrature1D quadrature;
  const double minJacobian = 4.0;
  quadrature.minJacobian(minJacobian);
  
  Integrator iOrig;
  iOrig.quadrature(&quadrature);

  Integrator iCopy = Integrator(iOrig);
  CPPUNIT_ASSERT_EQUAL(minJacobian, iCopy._quadrature->minJacobian());
} // testCopy

// ----------------------------------------------------------------------
// Test quadrature().
void
pylith::feassemble::TestIntegrator::testQuadrature(void)
{ // testQuadrature
  // Since quadrature is cloned, test setting quadrature by testing
  // value of minJacobian

  Quadrature1D quadrature;
  const double minJacobian = 4.0;
  quadrature.minJacobian(minJacobian);
  
  Integrator integrator;
  integrator.quadrature(&quadrature);

  CPPUNIT_ASSERT_EQUAL(minJacobian, integrator._quadrature->minJacobian());
} // testQuadrature


#if 0
// ----------------------------------------------------------------------
namespace pylith {
  namespace feassemble {
    class _TestIntegrator;
  } // feassemble
} // pylith

/// Helper class for TestIntegrator
class pylith::feassemble::_TestIntegrator {

public :
  /** Setup mesh.
   *
   * @param data Integrator data
   */
  static 
  ALE::Obj<ALE::Field::Mesh>
  _setupMesh(const IntegratorData& data);
}; // _TestIntegrator


// ----------------------------------------------------------------------
// Test integrateAction()
void
pylith::feassemble::TestIntegrator::_testIntegrateAction(Integrator* integrator,
					   const IntegratorData& data) const
{ // _testIntegrateAction
  CPPUNIT_ASSERT(false);

  typedef ALE::Field::Mesh        Mesh;
  typedef Mesh::real_section_type real_section_type;

  ALE::Obj<Mesh> mesh = _TestIntegrator::_setupMesh(data);

  // Fiber dimension (number of values in field per vertex) for fields
  const int fiberDim = data.fiberDim;

  // Setup input field for action
  const ALE::Obj<real_section_type>& fieldIn =
    mesh->getRealSection("fieldIn");
  fieldIn->setName("fieldIn");
  fieldIn->setFiberDimension(mesh->depthStratum(0), fiberDim);
  mesh->allocate(fieldIn);
  int iVertex = 0;
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  for (Mesh::label_sequence::iterator vIter=vertices->begin();
       vIter != verticesEnd;
       ++vIter, ++iVertex)
    fieldIn->updatePoint(*vIter, &data.fieldIn[iVertex*fiberDim]);

  // Setup field for action result
  const ALE::Obj<real_section_type>& fieldOut =
    mesh->getRealSection("fieldOut");
  fieldOut->setName("fieldOut");
  fieldOut->setFiberDimension(mesh->depthStratum(0), fiberDim);
  mesh->allocate(fieldOut);

  // Integrate action
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  integrator->integrateAction(fieldOut, fieldIn, coordinates);
  //fieldOut->view("field out");
  
  // Check values in output field
  iVertex = 0;
  const double tolerance = 1.0e-06;
  for (Mesh::label_sequence::iterator vIter=vertices->begin();
       vIter != verticesEnd;
       ++vIter, ++iVertex) {
    const real_section_type::value_type* vals = 
      fieldOut->restrictPoint(*vIter);
    const double* valsE = &data.valsAction[iVertex*fiberDim];
    const int dim = fieldOut->getFiberDimension(*vIter);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dim);
    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[iDim]/valsE[iDim], tolerance);
  } // for
} // _testIntegrateAction

// ----------------------------------------------------------------------
// Test integrate()
void
pylith::feassemble::TestIntegrator::_testIntegrate(Integrator* integrator,
					  const IntegratorData& data) const
{ // _testIntegrate
  CPPUNIT_ASSERT(false);

  typedef ALE::Field::Mesh        Mesh;
  typedef Mesh::real_section_type real_section_type;

  journal::debug_t debug("TestIntegrator");

  try {
    ALE::Obj<Mesh> mesh = _TestIntegrator::_setupMesh(data);

    // Fiber dimension (number of values in field per vertex) for fields
    const int fiberDim = data.fiberDim;

    // Setup input field for action
    const ALE::Obj<real_section_type>& fieldIn =
      mesh->getRealSection("fieldIn");
    fieldIn->setName("fieldIn");
    fieldIn->setFiberDimension(mesh->depthStratum(0), fiberDim);
    mesh->allocate(fieldIn);
    int iVertex = 0;
    const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
    const Mesh::label_sequence::iterator verticesEnd = vertices->end();
    for (topology_type::label_sequence::iterator vIter=vertices->begin();
	 vIter != verticesEnd;
	 ++vIter, ++iVertex)
      fieldIn->updatePoint(*vIter, &data.fieldIn[iVertex*fiberDim]);
    
    // Integrate
    PetscMat mat;
    const ALE::Obj<real_section_type>& coordinates = 
      mesh->getRealSection("coordinates");
    integrator->integrate(&mat, fieldIn, coordinates);

    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    
    // Check matrix size
    int nRows = 0;
    int nCols = 0;
    MatGetSize(mat, &nRows, &nCols);
    debug
      << journal::at(__HERE__)
      << "# rows: " << nRows
      << ", # cols: " << nCols
      << journal::endl;
    const int nRowsE = data.numVertices * data.fiberDim;
    const int nColsE = data.numVertices * data.fiberDim;
    CPPUNIT_ASSERT_EQUAL(nRowsE, nRows);
    CPPUNIT_ASSERT_EQUAL(nColsE, nCols);

    // Create dense matrix
    PetscMat matDense;
    PetscMat matSparseSeq;
    MatConvert(mat, MATSEQAIJ, MAT_INITIAL_MATRIX, &matSparseSeq);
    MatConvert(matSparseSeq, MATSEQDENSE, MAT_INITIAL_MATRIX, &matDense);
    MatDestroy(matSparseSeq);
    
    // Get values associated with dense matrix
    double* vals = 0;
    int* rows = 0;
    int* cols = 0;
    const int size = nRows*nCols;
    if (size > 0) {
      vals = new double[size];
      rows = new int[nRows];
      for (int iRow=0; iRow < nRows; ++iRow)
	rows[iRow] = iRow;
      cols = new int[nCols];
      for (int iCol=0; iCol < nCols; ++iCol)
	cols[iCol] = iCol;
    } // if
    MatGetValues(matDense, nRows, rows, nCols, cols, vals);
    delete[] rows; rows = 0;
    delete[] cols; cols = 0;

    // Check values from dense matrix
    const double tolerance = 1.0e-06;
    for (int iRow=0; iRow < nRows; ++iRow)
      for (int iCol=0; iCol < nCols; ++iCol) {
	const int index = iRow * nCols + iCol;
	debug
	  << journal::at(__HERE__)
	  << "index: " << index
	  << ", valsE: " << data.valsMatrix[index]
	  << ", vals: " << vals[index]
	  << journal::endl;
	if (fabs(data.valsMatrix[index]) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[index]/data.valsMatrix[index],
				       tolerance);
	else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(data.valsMatrix[index], vals[index],
				       tolerance);
      } // for
    delete[] vals; vals = 0;
    MatDestroy(matDense);
  } catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
    throw;
  } catch (ALE::Exception& err) {
    std::cerr << err << std::endl;
    throw;
  } // try/catch
} // _testIntegrate

// ----------------------------------------------------------------------
// Setup mesh.
ALE::Obj<ALE::Mesh>
pylith::feassemble::_TestIntegrator::_setupMesh(const IntegratorData& data)
{ // _setupMesh
  typedef ALE::Field::Mesh Mesh;
  typedef Meh::sieve_type  sieve_type;

  const int cellDim = data.cellDim;
  const int numCorners = data.numCorners;
  const int spaceDim = data.spaceDim;
  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const double* vertCoords = data.vertices;
  const int* cells = data.cells;
  CPPUNIT_ASSERT(0 != vertCoords);
  CPPUNIT_ASSERT(0 != cells);

  ALE::Obj<Mesh> mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

  const bool interpolate = false;
  ALE::New::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
  mesh->setSieve(sieve);
  mesh->stratify();
  ALE::New::SieveBuilder<Mesh>::buildCoordinatesNew(mesh, spaceDim, vertCoords);

  return mesh;
} // _setupMesh
#endif


// End of file 
