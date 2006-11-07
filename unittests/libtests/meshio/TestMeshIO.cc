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
// Check values in mesh against data.
void
pylith::meshio::TestMeshIO::_checkVals(const ALE::Obj<ALE::Mesh>& mesh,
				       const MeshData& data)
{ // _checkVals
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

  // :TODO: Check groups of vertices
} // _checkVals

// End of file 
