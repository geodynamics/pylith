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

#include "AbsorbingDampers.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cstring> // USES memcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealSection SubRealSection;
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RestrictVisitor RestrictVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
  _db(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::bc::AbsorbingDampers::deallocate(void)
{ // deallocate
  _db = 0; // :TODO: Use shared pointer
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::AbsorbingDampers::initialize(const topology::Mesh& mesh,
					 const double upDir[3])
{ // initialize
  assert(0 != _boundaryMesh);
  assert(0 != _quadrature);
  assert(0 != _db);

  double_array up(upDir, 3);
  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get damping constants at each quadrature point and rotate to
  // global coordinate frame using orientation information
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const int fiberDim = numQuadPts * spaceDim;
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  delete _parameters;
  _parameters = 
    new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
  assert(0 != _parameters);
  _parameters->add("damping constants", "damping_constants");
  topology::Field<topology::SubMesh>& dampingConsts =
    _parameters->get("damping constants");
  dampingConsts.newSection(cells, fiberDim);
  dampingConsts.allocate();

  logger.stagePop();

  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int jacobianSize = spaceDim * cellDim;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientation(orientationSize);

  // open database with material property information
  _db->open();
  int numValues = 0;
  if (cellDim > 0) {
    numValues = 3;
    const char* valueNames[] = { "density", "vp", "vs" };
    _db->queryVals(valueNames, numValues);
  } else {
    numValues = 2;
    const char* valueNames[] = { "density", "vp" };
    _db->queryVals(valueNames, numValues);
  } // else

  // Container for data returned in query of database
  double_array queryData(numValues);
  double_array quadPtRef(cellDim);
  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Container for damping constants for current cell
  double_array dampingConstsLocal(fiberDim);
  double_array dampingConstsGlobal(fiberDim);

  double_array coordinatesCell(numCorners*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveSubMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double densityScale = _normalizer->densityScale();
  assert(_normalizer->timeScale() > 0);
  const double velocityScale = 
    _normalizer->lengthScale() / _normalizer->timeScale();

  const ALE::Obj<SubRealSection>& dampersSection = dampingConsts.section();
  assert(!dampersSection.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();

  // Compute quadrature information
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_boundaryMesh, cells);
#endif

  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    const double_array& quadPtsNondim = _quadrature->quadPts();
    const double_array& quadPtsRef = _quadrature->quadPtsRef();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), 
				lengthScale);

    dampingConstsGlobal = 0.0;
    for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      // Compute damping constants in normal/tangential coordinates
      const int err = _db->query(&queryData[0], numValues, 
				 &quadPtsGlobal[iQuad*spaceDim], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[iQuad*spaceDim+i];
	msg << ") for absorbing boundary condition " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      // Nondimensionalize damping constants
      const double densityN = 
	_normalizer->nondimensionalize(queryData[0], densityScale);
      const double vpN = 
	_normalizer->nondimensionalize(queryData[1], velocityScale);
      const double vsN = (3 == numValues) ?
	_normalizer->nondimensionalize(queryData[2], velocityScale) :
	0.0;
      
      const double constTangential = densityN * vsN;
      const double constNormal = densityN * vpN;
      const int numTangential = spaceDim-1;
      for (int iDim=0; iDim < numTangential; ++iDim)
	dampingConstsLocal[iDim] = constTangential;
      dampingConstsLocal[spaceDim-1] = constNormal;

      // Compute normal/tangential orientation
      memcpy(&quadPtRef[0], &quadPtsRef[iQuad*cellDim], 
	     cellDim*sizeof(double));
#if defined(PRECOMPUTE_GEOMETRY)
      coordsVisitor.clear();
      sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
#endif
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      for (int iDim=0; iDim < spaceDim; ++iDim) {
	for (int jDim=0; jDim < spaceDim; ++jDim)
	  dampingConstsGlobal[iQuad*spaceDim+iDim] += 
	    dampingConstsLocal[jDim]*orientation[jDim*spaceDim+iDim];
	// Ensure damping constants are positive
	dampingConstsGlobal[iQuad*spaceDim+iDim] = 
	  fabs(dampingConstsGlobal[iQuad*spaceDim+iDim]);
      } // for
    } // for
    dampersSection->updatePoint(*c_iter, &dampingConstsGlobal[0]);
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _parameters);
  assert(0 != fields);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();
  double_array dampersCell(numQuadPts*spaceDim);

  // Get cell information
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<SubRealSection>& dampersSection =
    _parameters->get("damping constants").section();
  assert(!dampersSection.isNull());

  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());
  topology::SubMesh::UpdateAddVisitor residualVisitor(*residualSection,
						      &_cellVector[0]);

  const topology::Field<topology::Mesh>& dispT = 
    fields->get("disp(t)");
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());
  topology::Mesh::RestrictVisitor dispTVisitor(*dispTSection, 
					       numBasis*spaceDim);
  const topology::Field<topology::Mesh>& dispTmdt = 
    fields->get("disp(t-dt)");
  const ALE::Obj<RealSection>& dispTmdtSection = dispTmdt.section();
  assert(!dispTmdtSection.isNull());
  topology::Mesh::RestrictVisitor dispTmdtVisitor(*dispTmdtSection, 
						  numBasis*spaceDim);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveSubMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Get geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    dispTVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, dispTVisitor);
    const double* dispTCell = dispTVisitor.getValues();

    dispTmdtVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, dispTmdtVisitor);
    const double* dispTmdtCell = dispTmdtVisitor.getValues();

    dampersSection->restrictPoint(*c_iter, &dampersCell[0], dampersCell.size());

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] / (2.0 * dt);

      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += 
	      dampersCell[iQuad*spaceDim+iDim] *
	      valIJ * (dispTmdtCell[jBasis*spaceDim+iDim] - \
		       dispTCell[jBasis*spaceDim+iDim]);
        } // for
      } // for

    } // for
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(5*spaceDim))));

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveSubMesh->updateClosure(*c_iter, residualVisitor);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(
				      topology::Jacobian* jacobian,
				      const double t,
				      topology::SolutionFields* const fields)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != jacobian);
  assert(0 != fields);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Get cell information
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<SubRealSection>& dampersSection =
    _parameters->get("damping constants").section();
  assert(!dampersSection.isNull());

  const topology::Field<topology::Mesh>& solution = fields->solution();
  const ALE::Obj<SieveMesh>& sieveMesh = solution.mesh().sieveMesh();
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    solutionSection);
  assert(!globalOrder.isNull());
  topology::Mesh::IndicesVisitor jacobianVisitor(*solutionSection,
						 *globalOrder,
			   (int) pow(sieveMesh->getSieve()->getMaxConeSize(),
				     sieveMesh->depth())*spaceDim);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveSubMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    assert(numQuadPts*spaceDim == dampersSection->getFiberDimension(*c_iter));
    const double* dampingConstsCell = dampersSection->restrictPoint(*c_iter);

    // Compute Jacobian for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] / (2.0 * dt);
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim) {
            const int iBlock = (iBasis*spaceDim + iDim) * (numBasis*spaceDim);
            const int jBlock = (jBasis*spaceDim + iDim);
            _cellMatrix[iBlock+jBlock] += 
	      valIJ * dampingConstsCell[iQuad*spaceDim+iDim];
          } // for
        } // for
      } // for
    } // for
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
    
    // Assemble cell contribution into PETSc Matrix
    jacobianVisitor.clear();
    PetscErrorCode err = updateOperator(jacobianMat, *sieveSubMesh->getSieve(), 
					jacobianVisitor, *c_iter,
					&_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
  } // for

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(
			      const topology::Field<topology::Mesh>& jacobian,
			      const double t,
			      topology::SolutionFields* const fields)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != fields);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Get cell information
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();
  _initCellVector();

  // Get sections
  const ALE::Obj<SubRealSection>& dampersSection =
    _parameters->get("damping constants").section();
  assert(!dampersSection.isNull());

  const topology::Field<topology::Mesh>& solution = fields->solution();
  const ALE::Obj<SieveMesh>& sieveMesh = solution.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());

  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());
  topology::Mesh::UpdateAddVisitor jacobianVisitor(*jacobianSection, 
						   &_cellVector[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveSubMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    assert(numQuadPts*spaceDim == dampersSection->getFiberDimension(*c_iter));
    const double* dampingConstsCell = dampersSection->restrictPoint(*c_iter);

    // Compute Jacobian for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] / (2.0 * dt);
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim) {
            const int iBlock = (iBasis*spaceDim + iDim) * (numBasis*spaceDim);
            const int jBlock = (jBasis*spaceDim + iDim);
            _cellMatrix[iBlock+jBlock] += 
	      valIJ * dampingConstsCell[iQuad*spaceDim+iDim];
          } // for
        } // for
      } // for
    } // for
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
    _lumpCellMatrix();
    
    // Assemble cell contribution into lumped matrix.
    jacobianVisitor.clear();
    sieveSubMesh->updateClosure(*c_iter, jacobianVisitor);
  } // for

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BCIntegratorSubMesh::verifyConfiguration(mesh);
} // verifyConfiguration


// End of file 
