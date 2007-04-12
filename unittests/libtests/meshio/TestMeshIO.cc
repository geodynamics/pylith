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

#include "TestMeshIO.hh" // Implementation of class methods

#include "data/MeshData.hh"

// ----------------------------------------------------------------------
// Get simple mesh for testing I/O.
ALE::Obj<ALE::Mesh>*
pylith::meshio::TestMeshIO::createMesh(const MeshData& data)
{ // createMesh
  typedef ALE::Mesh Mesh;
  typedef Mesh::sieve_type sieve_type;
  typedef Mesh::label_type label_type;

  // buildTopology() requires zero based index
  CPPUNIT_ASSERT(true == data.useIndexZero);

  const int cellDim = data.cellDim;
  const int numCorners = data.numCorners;
  const int spaceDim = data.spaceDim;
  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const double* vertCoords = data.vertices;
  const int* cells = data.cells;
  const int* materialIds = data.materialIds;
  CPPUNIT_ASSERT(0 != vertCoords);
  CPPUNIT_ASSERT(0 != cells);
  CPPUNIT_ASSERT(0 != materialIds);

  ALE::Obj<Mesh>* meshHandle = new ALE::Obj<Mesh>;
  *meshHandle = new Mesh(PETSC_COMM_WORLD, cellDim);
  ALE::Obj<Mesh> mesh = *meshHandle;
  mesh.addRef();
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

  const bool interpolate = false;
  ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
					 const_cast<int*>(cells), numVertices,
					 interpolate, numCorners);
  mesh->setSieve(sieve);
  mesh->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);

  const ALE::Obj<Mesh::label_sequence>& cellsMesh = mesh->heightStratum(0);

  const ALE::Obj<label_type>& labelMaterials = mesh->createLabel("material-id");
  
  int i = 0;
  for(Mesh::label_sequence::iterator e_iter = 
	cellsMesh->begin();
      e_iter != cellsMesh->end();
      ++e_iter)
    mesh->setValue(labelMaterials, *e_iter, materialIds[i++]);

  return meshHandle;
} // createMesh

// ----------------------------------------------------------------------
// Check values in mesh against data.
void
pylith::meshio::TestMeshIO::checkVals(const ALE::Obj<Mesh>& mesh,
				      const MeshData& data)
{ // checkVals
  typedef ALE::Mesh::label_type label_type;

  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh->getDimension());

  // Check vertices
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<Mesh::real_section_type>& coordsField =
    mesh->getRealSection("coordinates");
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, 
		       coordsField->getFiberDimension(*vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for(Mesh::label_sequence::iterator v_iter = 
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const Mesh::real_section_type::value_type *vertexCoords = 
      coordsField->restrictPoint(*v_iter);
    const double tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      if (data.vertices[i] < 1.0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[i++], vertexCoords[iDim],
				   tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vertexCoords[iDim]/data.vertices[i++],
				   tolerance);
      }
  } // for

  // check cells
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);
  const int numCorners = sieve->nCone(*cells->begin(), 
				      mesh->depth())->size();
  CPPUNIT_ASSERT_EQUAL(data.numCorners, numCorners);

  const ALE::Obj<Mesh::numbering_type>& vNumbering = 
    mesh->getFactory()->getLocalNumbering(mesh, 0);

  const int offset = (data.useIndexZero) ? 0 : 1;
  i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter) {
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*e_iter);
    for(sieve_type::traits::coneSequence::iterator c_iter = cone->begin();
	c_iter != cone->end();
	++c_iter)
      CPPUNIT_ASSERT_EQUAL(data.cells[i++], 
			   vNumbering->getIndex(*c_iter) + offset);
  } // for

  // check materials
  const int size = numCells;
  int* materialIds = (size > 0) ? new int[size] : 0;
  const ALE::Obj<label_type>& labelMaterials = 
    mesh->getLabel("material-id");
  const int idDefault = -999;

  i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    materialIds[i++] = mesh->getValue(labelMaterials, *e_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);
  delete[] materialIds; materialIds = 0;

  // :TODO: Check groups of vertices
} // checkVals

// End of file 
