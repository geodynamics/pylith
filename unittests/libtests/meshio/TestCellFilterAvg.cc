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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/CellFilterAvg.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestCellFilterAvg );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestCellFilterAvg::testConstructor(void)
{ // testConstructor
  CellFilterAvg<topology::Mesh> filter;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter()
void
pylith::meshio::TestCellFilterAvg::testFilter(void)
{ // testFilter
  typedef pylith::topology::Mesh::SieveMesh SieveMesh;
  typedef pylith::topology::Mesh::RealSection RealSection;

  const char* filename = "data/quad4.mesh";
  const int fiberDim = 2;
  const int ncells = 2;
  const std::string label = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = 
    topology::FieldBase::MULTI_SCALAR;
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

  const topology::FieldBase::VectorFieldEnum fieldTypeE = 
    topology::FieldBase::SCALAR;
  const int fiberDimE = 1;
  const double fieldValuesE[] = {
    (1.5*1.1 + 0.5*1.2)/2.0,
    (1.5*2.1 + 0.5*2.2)/2.0,
  };

  MeshIOAscii iohandler;
  topology::Mesh mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  // Set cell field
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  field.allocate();
  field.vectorFieldType(fieldType);
  field.label(label.c_str());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  CPPUNIT_ASSERT_EQUAL(ncells, int(cells->size()));
  int ipt = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++ipt) {
    const double* values = &fieldValues[ipt*fiberDim];
    section->updatePoint(*c_iter, values);
  } // for

  feassemble::Quadrature<topology::Mesh> quadrature;
  quadrature.initialize(basis, numQuadPts, numBasis,
			basisDerivRef, numQuadPts, numBasis, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  CellFilterAvg<topology::Mesh> filter;
  filter.quadrature(&quadrature);

  const topology::Field<topology::Mesh>& fieldF = filter.filter(field);
  const ALE::Obj<RealSection>& sectionF = fieldF.section();
  CPPUNIT_ASSERT(!sectionF.isNull());

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldF.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(label, std::string(fieldF.label()));

  ipt = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++ipt) {
    CPPUNIT_ASSERT_EQUAL(fiberDimE, sectionF->getFiberDimension(*c_iter));
    const double* values = sectionF->restrictPoint(*c_iter);
    CPPUNIT_ASSERT(0 != values);
    const double tolerance = 1.0e-06;
    for (int i=0; i < fiberDimE; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				   values[i]/fieldValuesE[ipt*fiberDimE+i],
				   tolerance);
  } // for
} // testFilter


// End of file 
