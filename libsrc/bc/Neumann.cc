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

#include "Neumann.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA Field
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cstring> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealSection SubRealSection;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann::Neumann(void) :
  _boundaryMesh(0),
  _tractions(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void)
{ // destructor
  delete _boundaryMesh; _boundaryMesh = 0;
  delete _tractions; _tractions = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::Neumann::initialize(const topology::Mesh& mesh,
				const double upDir[3])
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _db);

  delete _boundaryMesh; _boundaryMesh = 0;
  delete _tractions; _tractions = 0;

  _boundaryMesh = new topology::SubMesh(mesh, _label.c_str());
  assert(0 != _boundaryMesh);

  double_array up(upDir, 3);

  // check compatibility of quadrature and boundary mesh
  if (_quadrature->cellDim() != _boundaryMesh->dimension()) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for Neumann traction "
	<< "boundary condition '" << _label << "'.\n"
	<< "Dimension of boundary mesh: " << _boundaryMesh->dimension()
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& submesh = _boundaryMesh->sieveMesh();
  assert(!submesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    submesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const int boundaryDepth = submesh->depth()-1; // depth of bndry cells

  // Make sure surface cells are compatible with quadrature.
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = 
      submesh->getNumCellCorners(*c_iter, boundaryDepth);
    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell for Neumann traction "
	  << "boundary condition '" << _label << "'.\n"
	  << "Cell " << *c_iter << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  // Create section for traction vector in global coordinates
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const int fiberDim = spaceDim * numQuadPts;
  
  _tractions = new topology::Field<topology::SubMesh>(*_boundaryMesh);
  assert(0 != _tractions);
  _tractions->newSection(cells, fiberDim);
  _tractions->allocate();

  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int jacobianSize = spaceDim * cellDim;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientation(orientationSize);
  double_array cellVertices(numCorners*spaceDim);

  // Set names based on dimension of problem.
  // 1-D problem = {'normal-traction'}
  // 2-D problem = {'shear-traction', 'normal-traction'}
  // 3-D problem = {'horiz-shear-traction', 'vert-shear-traction',
  //                'normal-traction'}
  _db->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* valueNames[] = {"normal-traction"};
      _db->queryVals(valueNames, 1);
      break;
    } // case 1
    case 2 : {
      const char* valueNames[] = {"shear-traction", "normal-traction"};
      _db->queryVals(valueNames, 2);
      break;
    } // case 2
    case 3 : {
      const char* valueNames[] = {"horiz-shear-traction",
				  "vert-shear-traction",
				  "normal-traction"};
      _db->queryVals(valueNames, 3);
      break;
    } // case 3
    default :
      assert(0);
    } // switch

  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  double_array tractionDataLocal(spaceDim);
  double_array quadPtRef(cellDim);
  double_array quadPtsGlobal(numQuadPts*spaceDim);
  const double_array& quadPtsRef = _quadrature->quadPtsRef();

  // Container for cell tractions rotated to global coordinates.
  double_array cellTractionsGlobal(fiberDim);

  // Get sections.
  const ALE::Obj<RealSection>& coordinates =
    submesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						numCorners*spaceDim);

  const ALE::Obj<SubRealSection>& tractSection = _tractions->section();
  assert(!tractSection.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double pressureScale = _normalizer->pressureScale();

  // Compute quadrature information
  _quadrature->computeGeometry(*_boundaryMesh, cells);

  // Loop over cells in boundary mesh, compute orientations, and then
  // compute corresponding traction vector in global coordinates
  // (store values in _tractionGlobal).
  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cellsEnd;
      ++c_iter) {
    _quadrature->retrieveGeometry(*c_iter);
    const double_array& quadPtsNondim = _quadrature->quadPts();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    coordsVisitor.clear();
    submesh->restrictClosure(*c_iter, coordsVisitor);
    
    cellTractionsGlobal = 0.0;
    for(int iQuad=0, iRef=0, iSpace=0; iQuad < numQuadPts;
	++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      // Get traction vector in local coordinate system at quadrature point
      const int err = _db->query(&tractionDataLocal[0], spaceDim,
				 &quadPtsGlobal[iSpace], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find traction values at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << " " << quadPtsGlobal[i+iSpace];
	msg << ") for traction boundary condition " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      _normalizer->nondimensionalize(&tractionDataLocal[0], spaceDim,
				     pressureScale);

      // Compute Jacobian and determinant at quadrature point, then get
      // orientation.
      memcpy(&quadPtRef[0], &quadPtsRef[iRef], cellDim*sizeof(double));
      memcpy(&cellVertices[0], coordsVisitor.getValues(), 
	     cellVertices.size()*sizeof(double));
      cellGeometry.jacobian(&jacobian, &jacobianDet, cellVertices, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      orientation /= jacobianDet;

      // Rotate traction vector from local coordinate system to global
      // coordinate system
      for(int iDim = 0; iDim < spaceDim; ++iDim) {
	for(int jDim = 0; jDim < spaceDim; ++jDim)
	  cellTractionsGlobal[iDim+iSpace] +=
	    orientation[jDim*spaceDim+iDim] * tractionDataLocal[jDim];
      } // for
    } // for

      // Update tractionsGlobal
    tractSection->updatePoint(*c_iter, &cellTractionsGlobal[0]);
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::Neumann::integrateResidual(const topology::Field<topology::Mesh>& residual,
				       const double t,
				       topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _tractions);

  PetscErrorCode err = 0;

  // Get cell information
  const ALE::Obj<SieveSubMesh>& submesh = _boundaryMesh->sieveMesh();
  assert(!submesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    submesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<SubRealSection>& tractSection = _tractions->section();
  assert(!tractSection.isNull());
  const ALE::Obj<RealSection>& residualSection = residual.section();

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();
  const int cellVecSize = numBasis*spaceDim;
  double_array tractionsCell(numQuadPts*spaceDim);

  // Loop over faces and integrate contribution from each face
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Get geometry information for current cell
    _quadrature->retrieveGeometry(*c_iter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict tractions to cell
    tractSection->restrictPoint(*c_iter, 
				&tractionsCell[0], tractionsCell.size());

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for traction bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += 
	      tractionsCell[iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for
    // Assemble cell contribution into field
    // :TODO: Use an UpdateAddVisitor?
    submesh->updateAdd(residualSection, *c_iter, &_cellVector[0]);

    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::Neumann::integrateJacobian(PetscMat* jacobian,
 				       const double t,
 				       topology::SolutionFields* const fields)
{ // integrateJacobian
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BoundaryCondition::verifyConfiguration(mesh);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get boundary mesh.
const pylith::topology::SubMesh&
pylith::bc::Neumann::boundaryMesh(void) const
{ // dataMesh
  assert(0 != _boundaryMesh);

  return *_boundaryMesh;
} // dataMesh

// ----------------------------------------------------------------------
// Get cell field for tractions.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::Neumann::cellField(const char* name,
			       topology::SolutionFields* const fields)
{ // cellField
  assert(0 != _tractions);
  assert(0 != name);

  if (0 == strcasecmp(name, "tractions")) {
    return *_tractions;
  } else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for Neumann BC '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  return *_tractions;
} // cellField


// End of file 
