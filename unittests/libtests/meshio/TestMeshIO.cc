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

#include <Mesh.hh>

#include "data/MeshData.hh"

// ----------------------------------------------------------------------
// Get simple mesh for testing I/O.
ALE::Obj<ALE::Mesh>*
pylith::meshio::TestMeshIO::createMesh(const MeshData& data)
{ // createMesh
  typedef ALE::Mesh::topology_type topology_type;
  typedef topology_type::sieve_type sieve_type;
  typedef ALE::Sifter<int, sieve_type::point_type, int> patch_label_type;

  // buildTopology() requires zero based index
  assert(true == data.useIndexZero);

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

  ALE::Obj<ALE::Mesh>* meshHandle = new ALE::Obj<ALE::Mesh>;
  *meshHandle = new ALE::Mesh(PETSC_COMM_WORLD, cellDim);
  ALE::Obj<ALE::Mesh> mesh = *meshHandle;
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
  ALE::Obj<topology_type> topology = new topology_type(mesh->comm());

  const bool interpolate = false;
  ALE::New::SieveBuilder<sieve_type>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
  sieve->stratify();
  topology->setPatch(0, sieve);
  topology->stratify();
  mesh->setTopology(topology);
  ALE::New::SieveBuilder<sieve_type>::buildCoordinates(
		    mesh->getRealSection("coordinates"), spaceDim, vertCoords);

  const Mesh::real_section_type::patch_type patch = 0;
  const ALE::Obj<Mesh::topology_type::label_sequence>& cellsMesh = 
    topology->heightStratum(patch, 0);

  const ALE::Obj<patch_label_type>& labelMaterials = 
    topology->createLabel(patch, "material-id");
  
  int i = 0;
  for(Mesh::topology_type::label_sequence::iterator e_iter = 
	cellsMesh->begin();
      e_iter != cellsMesh->end();
      ++e_iter)
    topology->setValue(labelMaterials, *e_iter, materialIds[i++]);

  return meshHandle;
} // createMesh

// ----------------------------------------------------------------------
// Check values in mesh against data.
void
pylith::meshio::TestMeshIO::checkVals(const ALE::Obj<ALE::Mesh>& mesh,
				      const MeshData& data)
{ // checkVals
  typedef ALE::Sifter<int, sieve_type::point_type, int> patch_label_type;

  const Mesh::real_section_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = mesh->getTopology();

  // Check mesh dimension
  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh->getDimension());

  // Check vertices
  const ALE::Obj<Mesh::topology_type::label_sequence>& vertices = 
    topology->depthStratum(patch, 0);
  const ALE::Obj<Mesh::real_section_type>& coordsField =
    mesh->getRealSection("coordinates");
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, 
		       coordsField->getFiberDimension(patch, 
						      *vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for(Mesh::topology_type::label_sequence::iterator v_iter = 
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const Mesh::real_section_type::value_type *vertexCoords = 
      coordsField->restrict(patch, *v_iter);
    const double tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vertexCoords[iDim]/data.vertices[i++],
				   tolerance);
  } // for

  // check cells
  const ALE::Obj<sieve_type>& sieve = topology->getPatch(patch);
  const ALE::Obj<Mesh::topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);
  const int numCorners = sieve->nCone(*cells->begin(), 
				      topology->depth())->size();
  CPPUNIT_ASSERT_EQUAL(data.numCorners, numCorners);

  const ALE::Obj<Mesh::numbering_type>& vNumbering = 
    mesh->getFactory()->getLocalNumbering(topology, patch, 0);

  const int offset = (data.useIndexZero) ? 0 : 1;
  i = 0;
  for(Mesh::topology_type::label_sequence::iterator e_iter = cells->begin();
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
  const ALE::Obj<patch_label_type>& labelMaterials = 
    topology->getLabel(patch, "material-id");
  const int idDefault = -999;

  i = 0;
  for(Mesh::topology_type::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    materialIds[i++] = topology->getValue(labelMaterials, *e_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);
  delete[] materialIds; materialIds = 0;

  // :TODO: Check groups of vertices
} // checkVals

// End of file 
