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

#include "Mesh.hh" // USES Mesh

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
template<typename Point>
class Edge : public std::pair<Point, Point> {
public:
  Edge() : std::pair<Point, Point>() {};
  Edge(const Point l) : std::pair<Point, Point>(l, l) {};
  Edge(const Point l, const Point r) : std::pair<Point, Point>(l, r) {};
  ~Edge() {};
  friend std::ostream& operator<<(std::ostream& stream, const Edge& edge) {
    stream << "(" << edge.first << ", " << edge.second << ")";
    return stream;
  };
};

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
pylith::topology::RefineUniform::refine(Mesh* const newMesh,
					const Mesh& mesh,
					const int levels,
					const SolutionFields* fields)
{ // refine
  assert(0 != newMesh);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  assert(!cells.isNull());

  newMesh->debug(mesh.debug());

  // Assume number of corners per cell is the same for the entire mesh
  assert(cells->size() > 0);
  const int cellNumCorners = sieveMesh->getNumCellCorners(*cells->begin());

  if (3 == mesh.dimension() && 4 == cellNumCorners)
    _refineTet4(newMesh, mesh, levels);

  // TODO: Add other refinement cases here

  else {
    std::ostringstream msg;
    msg << "Unknown case for uniform global refinement.\n"
	<< "mesh dimension: " << mesh.dimension()
	<< ", number of corners in cell: " << cellNumCorners;
    throw std::runtime_error(msg.str());
  } // else
} // refine
    

// ----------------------------------------------------------------------
// Refine tet4 mesh.
void
pylith::topology::RefineUniform::_refineTet4(Mesh* const newMesh,
					     const Mesh& mesh,
					     const int levels)
{ // _refineTet4
  assert(0 != newMesh);

  typedef SieveMesh::point_type point_type;
  typedef Edge<point_type> edge_type;

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh>& newSieveMesh = newMesh->sieveMesh();
  assert(!newSieveMesh.isNull());

  ALE::Obj<SieveMesh::sieve_type> newSieve =
    new SieveMesh::sieve_type(mesh.comm(), mesh.debug());

  std::map<edge_type, point_type> edge2vertex;
    
  newSieveMesh->setSieve(newSieve);
  ALE::MeshBuilder<Mesh>::refineTetrahedra(*mesh.sieveMesh(), *newSieveMesh,
					   edge2vertex);

  // Fix material ids
  const int numCells = sieveMesh->heightStratum(0)->size();
  const ALE::Obj<SieveMesh::label_type>& materials =
    sieveMesh->getLabel("material-id");
  const ALE::Obj<SieveMesh::label_type>& newMaterials =
    newSieveMesh->createLabel("material-id");
  
  const int numNewCellsPerCell = 8; // :KLUDGE: depends on levels
  for(int icell=0; icell < numCells; ++icell) {
    const int material = sieveMesh->getValue(materials, icell);
    
    for(int i=0; i < numNewCellsPerCell; ++i)
      newSieveMesh->setValue(newMaterials, icell*8+i, material);
  } // for
  
  // Fix groups, assuming vertex groups
  const int numNewVertices = newSieveMesh->depthStratum(0)->size();
  const int numNewCells = newSieveMesh->heightStratum(0)->size();
  const ALE::Obj<std::set<std::string> >& sectionNames =
    sieveMesh->getIntSections();
  
  const std::set<std::string>::const_iterator namesBegin = 
    sectionNames->begin();
  const std::set<std::string>::const_iterator namesEnd = 
    sectionNames->end();
  for (std::set<std::string>::const_iterator name=namesBegin;
      name != namesEnd;
      ++name) {
    const ALE::Obj<Mesh::IntSection>& group = sieveMesh->getIntSection(*name);
    const ALE::Obj<Mesh::IntSection>& newGroup =
      newSieveMesh->getIntSection(*name);
    const Mesh::IntSection::chart_type& chart = group->getChart();
      
    newGroup->setChart(Mesh::IntSection::chart_type(numNewCells, 
						    numNewCells + numNewVertices));
    const Mesh::IntSection::chart_type& newChart = newGroup->getChart();
      
    const int chartMax = chart.max();
    for (int p = chart.min(), pNew = newChart.min(); p < chartMax; ++p, ++pNew) {
      if (group->getFiberDimension(p))
	newGroup->setFiberDimension(pNew, 1);
    } // for
    const std::map<edge_type, point_type>::const_iterator edge2VertexEnd =
      edge2vertex.end();
    for (std::map<edge_type, point_type>::const_iterator e_iter=edge2vertex.begin();
	 e_iter != edge2VertexEnd;
	 ++e_iter) {
      const point_type vertexA = e_iter->first.first;
      const point_type vertexB = e_iter->first.second;
      
      if (group->getFiberDimension(vertexA) && group->getFiberDimension(vertexB))
	if (group->restrictPoint(vertexA)[0] == group->restrictPoint(vertexB)[0])
	  newGroup->setFiberDimension(e_iter->second, 1);
    } // for

    newGroup->allocatePoint();
    for (int p=chart.min(), pNew = newChart.min(); p < chartMax; ++p, ++pNew) {
      if (group->getFiberDimension(p))
	newGroup->updatePoint(pNew, group->restrictPoint(p));
    } // for
    for (std::map<edge_type, point_type>::const_iterator e_iter=edge2vertex.begin();
	 e_iter != edge2VertexEnd;
	 ++e_iter) {
      const point_type vertexA = e_iter->first.first;
      const point_type vertexB = e_iter->first.second;
	
      if (group->getFiberDimension(vertexA) && group->getFiberDimension(vertexB))
	if (group->restrictPoint(vertexA)[0] == group->restrictPoint(vertexB)[0])
	  newGroup->updatePoint(e_iter->second, group->restrictPoint(vertexA));
    } // for
  } // for
} // _refineTet4


// End of file 
