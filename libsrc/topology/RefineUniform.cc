// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "RefineUniform.hh" // implementation of class methods

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::topology::RefineUniform::RefineUniform(void)
{ // constructor
} // constructor
 
// ----------------------------------------------------------------------
// Destructor
pylith::topology::RefineUniform::~RefineUniform(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Refine mesh.
void
pylith::topology::RefineUniform::refine(ALE::Obj<Mesh>* const newMesh,
					const ALE::Obj<Mesh>& mesh,
					const int levels)
{ // refine
  typedef Mesh::point_type                 point_type;
  typedef std::pair<point_type,point_type> edge_type;

  *newMesh = new Mesh(mesh->comm(), mesh->getDimension(), mesh->debug());
  Obj<Mesh::sieve_type> newSieve = new Mesh::sieve_type(mesh->comm(), mesh->debug());
  std::map<edge_type, point_type> edge2vertex;

  (*newMesh)->setSieve(newSieve);
  ALE::MeshBuilder<Mesh>::refineTetrahedra(*mesh, *(*newMesh), edge2vertex);

  // Fix material ids
  const int                         numCells     = mesh->heightStratum(0)->size();
  const ALE::Obj<Mesh::label_type>& materials    = mesh->getLabel("material-id");
  const ALE::Obj<Mesh::label_type>& newMaterials = (*newMesh)->createLabel("material-id");

  for(int c = 0; c < numCells; ++c) {
    const int material = mesh->getValue(materials, c);

    for(int i = 0; i < 8; ++i) {(*newMesh)->setValue(newMaterials, c*8+i, material);}
  }
  // Fix groups, assuming vertex groups
  const int                               numNewVertices = (*newMesh)->depthStratum(0)->size();
  const int                               numNewCells    = (*newMesh)->heightStratum(0)->size();
  const ALE::Obj<std::set<std::string> >& sectionNames   = mesh->getIntSections();

  for(std::set<std::string>::const_iterator name = sectionNames->begin(); name != sectionNames->end(); ++name) {
    const ALE::Obj<Mesh::int_section_type>&   group    = mesh->getIntSection(*name);
    const ALE::Obj<Mesh::int_section_type>&   newGroup = (*newMesh)->getIntSection(*name);
    const Mesh::int_section_type::chart_type& chart    = group->getChart();

    group->setChart(Mesh::int_section_type::chart_type(numNewCells, numNewCells + numNewVertices));
    for(int p = chart.min(); p < chart.max(); ++p) {
      if (group->getFiberDimension(p)) {
        newGroup->setFiberDimension(p, 1);
      }
    }
    for(std::map<edge_type, point_type>::const_iterator e_iter = edge2vertex.begin(); e_iter != edge2vertex.end(); ++e_iter) {
      const point_type vertexA = e_iter->first.first;
      const point_type vertexB = e_iter->first.second;

      if (group->getFiberDimension(vertexA) && group->getFiberDimension(vertexB)) {
        if (group->restrictPoint(vertexA)[0] == group->restrictPoint(vertexB)[0]) {
          newGroup->setFiberDimension(e_iter->second, 1);
        }
      }
    }
    newGroup->allocatePoint();
    for(int p = chart.min(); p < chart.max(); ++p) {
      if (group->getFiberDimension(p)) {
        newGroup->updatePoint(p, group->restrictPoint(p));
      }
    }
    for(std::map<edge_type, point_type>::const_iterator e_iter = edge2vertex.begin(); e_iter != edge2vertex.end(); ++e_iter) {
      const point_type vertexA = e_iter->first.first;
      const point_type vertexB = e_iter->first.second;

      if (group->getFiberDimension(vertexA) && group->getFiberDimension(vertexB)) {
        if (group->restrictPoint(vertexA)[0] == group->restrictPoint(vertexB)[0]) {
          newGroup->updatePoint(e_iter->second, group->restrictPoint(vertexA));
        }
      }
    }
  }
} // refine


// End of file 
