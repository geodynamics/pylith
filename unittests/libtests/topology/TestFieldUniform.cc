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

#include "TestFieldUniform.hh" // Implementation of class methods
#include "pylith/topology/FieldUniform.hh" // USES FieldUniform

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldUniform );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestFieldUniform {
      const int spaceDim = 2;
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 1;
      const int ncorners = 4;
      const int cells[] = { 0, 1, 2, 3 };
      const double coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
    } // _TestFieldUniform
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldUniform::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  mesh.createSieveMesh();
  const int fiberDim = 4;
  FieldUniform field(mesh.sieveMesh(), fiberDim);

  CPPUNIT_ASSERT(!field._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(fiberDim, field._fiberDim);
} // testConstructor

// ----------------------------------------------------------------------
// Test createSection() with points.
void
pylith::topology::TestFieldUniform::testCreateSectionPoints(void)
{ // testCreateSectionPoints
  const int fiberDim = 3;
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  FieldUniform field(sieveMesh, fiberDim);
  field.createSection(vertices);

  const ALE::Obj<SieveRealSection>& section = field.section();
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testCreateSectionPoints

// ----------------------------------------------------------------------
// Test createSection() with chart.
void
pylith::topology::TestFieldUniform::testCreateSectionChart(void)
{ // testCreateSectionChart
  const int fiberDim = 2;
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  // Create field with to use to create new field
  FieldUniform fieldSrc(sieveMesh, fiberDim);
  { // Setup source field
    fieldSrc.newSection();
    const ALE::Obj<SieveRealSection>& section = fieldSrc.section();
    section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
    section->setFiberDimension(vertices, fiberDim);
    sieveMesh->allocate(section);
  } // Setup source field

  FieldUniform field(sieveMesh, fiberDim);
  field.createSection(fieldSrc.section()->getChart());

  const ALE::Obj<SieveRealSection>& section = field.section();
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testCreateSectionChart

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldUniform::_buildMesh(Mesh* mesh)
{ // _buildMesh
  assert(0 != mesh);

  mesh->createSieveMesh(_TestFieldUniform::cellDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();

  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(sieveMesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  ALE::Obj<ALE::Mesh::sieve_type> s = 
    new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());
  
  const int cellDim = _TestFieldUniform::cellDim;
  const int ncells = _TestFieldUniform::ncells;
  const int* cells = _TestFieldUniform::cells;
  const int nvertices = _TestFieldUniform::nvertices;
  const int ncorners = _TestFieldUniform::ncorners;
  const int spaceDim = _TestFieldUniform::spaceDim;
  const double* coordinates = _TestFieldUniform::coordinates;
  const bool interpolate = false;
  ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, ncells, (int*) cells,
					      nvertices, interpolate, 
					      ncorners);
  std::map<SieveMesh::point_type,SieveMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 coordinates);

} // _buildMesh

// End of file 
