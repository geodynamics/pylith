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

#include "TestCellFilterAvg.hh" // Implementation of class methods

#include "pylith/meshio/CellFilterAvg.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature2D.hh" // USES Quadrature

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestCellFilterAvg );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestCellFilterAvg::testConstructor(void)
{ // testConstructor
  CellFilterAvg filter;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter()
void
pylith::meshio::TestCellFilterAvg::testFilter(void)
{ // testFilter
  const char* filename = "data/quad4.mesh";
  const int fiberDim = 2;
  const int ncells = 2;
  const char* fieldName = "field data";
  const VectorFieldEnum fieldType = OTHER_FIELD;
  const double fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
  };
  const int cellDim = 2;
  const int numBasis = 4;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const double basis[] = {
    1.0, 1.0,
    1.0, 1.0,
    1.0, 1.0,
    1.0, 1.0,
  };
  const double basisDerivRef[] = {
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
  };
  const double quadPtsRef[] = {
    1.0, 0.0,
   -1.0, 0.0,};
  const double quadWts[] = { 1.5, 0.5 };
  const double minJacobian = 1.0;



  const VectorFieldEnum fieldTypeE = OTHER_FIELD;
  const int fiberDimE = 1;
  const double fieldValuesE[] = {
    (1.5*1.1 + 0.5*1.2)/2.0,
    (1.5*2.1 + 0.5*2.2)/2.0,
  };

  MeshIOAscii iohandler;
  topology::Mesh(PETSC_COMM_WORLD, cellDim);
  iohandler.filename(filename);
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  // Set cell field
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  ALE::Obj<SieveMesh::real_section_type> field = 
    new SieveMesh::real_section_type(mesh.comm(), mesh.debug());
  field->setChart(SieveMesh::real_section_type::chart_type(0, cells->size()));
  field->setFiberDimension(cells, fiberDim);
  sieveMesh->allocate(field);

  CPPUNIT_ASSERT_EQUAL(ncells, int(cells->size()));

  int ipt = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++ipt) {
    const double* values = &fieldValues[ipt*fiberDim];
    field->updatePoint(*c_iter, values);
  } // for

  feassemble::Quadrature2D quadrature;
  quadrature.initialize(basis, basisDerivRef, quadPtsRef, quadWts,
			cellDim, numBasis, numQuadPts, spaceDim);

  CellFilterAvg filter;
  filter.quadrature(&quadrature);

  VectorFieldEnum fieldTypeF = SCALAR_FIELD;
  const ALE::Obj<real_section_type>& fieldF =
    filter.filter(&fieldTypeF, field, meshMesh);

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldTypeF);
  ipt = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++ipt) {
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldF->getFiberDimension(*c_iter));
    const double* values = fieldF->restrictPoint(*c_iter);
    CPPUNIT_ASSERT(0 != values);
    const double tolerance = 1.0e-06;
    for (int i=0; i < fiberDimE; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				   values[i]/fieldValuesE[ipt*fiberDimE+i],
				   tolerance);
  } // for
} // testFilter


// End of file 
