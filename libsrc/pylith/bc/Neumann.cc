// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Neumann.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/FieldsNew.hh" // HOLDSA FieldsNew
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;
typedef pylith::topology::Mesh::RealSection RealSection;

typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::SubMesh>::UpdateAddVisitor UpdateAddVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann::Neumann(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::Neumann::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::Neumann::initialize(const topology::Mesh& mesh,
				const PylithScalar upDir[3])
{ // initialize
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  _queryDatabases();
  _paramsLocalToGlobal(upDir);

  logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::Neumann::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const PylithScalar t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_parameters);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();
  scalar_array tractionsCell(numQuadPts*spaceDim);

  // Get cell information
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  _calculateValue(t);
  const ALE::Obj<SubRealUniformSection>& parametersSection =
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(valueFiberDim == tractionsCell.size());
  assert(valueIndex+valueFiberDim <= parametersFiberDim);

  const ALE::Obj<RealSection>& residualSection = residual.section();
  UpdateAddVisitor residualVisitor(*residualSection, &_cellVector[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    subSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Loop over faces and integrate contribution from each face
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict tractions to cell
    assert(parametersFiberDim == 
	   parametersSection->getFiberDimension(*c_iter));
    const PylithScalar* parametersCell = parametersSection->restrictPoint(*c_iter);
    assert(parametersCell);
    const PylithScalar* tractionsCell = &parametersCell[valueIndex];
    assert(tractionsCell);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for traction bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const PylithScalar valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const PylithScalar valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += 
	      tractionsCell[iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for
    // Assemble cell contribution into field
    residualVisitor.clear();
    subSieveMesh->updateClosure(*c_iter, residualVisitor);

    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  if (1 == mesh.dimension())
    throw std::runtime_error("Neumann boundary conditions are not "
			     "implemented for a 1-D mesh.");

  BCIntegratorSubMesh::verifyConfiguration(mesh);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get cell field for tractions.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::Neumann::cellField(const char* name,
			       topology::SolutionFields* const fields)
{ // cellField
  assert(_parameters);
  assert(name);

  if (0 == strcasecmp(name, "initial-value"))
    return _parameters->get("initial");

  else if (0 == strcasecmp(name, "rate-of-change"))
    return _parameters->get("rate");

  else if (0 == strcasecmp(name, "change-in-value"))
    return _parameters->get("change");

  else if (0 == strcasecmp(name, "rate-start-time"))
    return _parameters->get("rate time");

  else if (0 == strcasecmp(name, "change-start-time"))
    return _parameters->get("change time");

  else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for Neumann BC '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  return _parameters->get("traction"); // Satisfy method definition
} // cellField

// ----------------------------------------------------------------------
// Query databases for parameters.
void
pylith::bc::Neumann::_queryDatabases(void)
{ // _queryDatabases
  assert(_quadrature);
  assert(_boundaryMesh);
  
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar rateScale = pressureScale / timeScale;

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();

  delete _parameters; 
  _parameters = 
    new topology::FieldsNew<topology::SubMesh>(*_boundaryMesh);

  // Create section to hold time dependent values
  _parameters->add("value", "traction", 
		   numQuadPts*spaceDim, topology::FieldBase::MULTI_VECTOR,
		   pressureScale);
  if (_dbInitial) 
    _parameters->add("initial", "initial_traction",
		     numQuadPts*spaceDim, topology::FieldBase::MULTI_VECTOR,
		     pressureScale);
  if (_dbRate) {
    _parameters->add("rate", "traction_rate",
		     numQuadPts*spaceDim, topology::FieldBase::MULTI_VECTOR,
		     rateScale);
    _parameters->add("rate time", "traction_rate__time",
		     numQuadPts, topology::FieldBase::MULTI_SCALAR,
		     timeScale);
  } // if
  if (_dbChange) {
    _parameters->add("change", "change_traction",
		     numQuadPts*spaceDim, topology::FieldBase::MULTI_VECTOR,
		     pressureScale);
    _parameters->add("change time", "change_traction_time",
		     numQuadPts, topology::FieldBase::MULTI_SCALAR,
		     timeScale);
  } // if
  _parameters->allocate(topology::FieldBase::CELLS_FIELD, 1);
  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());

  if (0 != _dbInitial) { // Setup initial values, if provided.
    _dbInitial->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbInitial->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbInitial->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-horiz",
				    "traction-shear-vert",
				    "traction-normal"};
	_dbInitial->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Neumann.");
      } // switch
    _queryDB("initial", _dbInitial, spaceDim, pressureScale);
    _dbInitial->close();
  } // if

  if (0 != _dbRate) { // Setup rate of change of values, if provided.
    _dbRate->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-rate-normal"};
	_dbRate->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-rate-shear", 
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-rate-shear-horiz",
				    "traction-rate-shear-vert",
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Neumann.");
      } // switch
    _queryDB("rate", _dbRate, spaceDim, rateScale);

    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB("rate time", _dbRate, 1, timeScale);
    _dbRate->close();
  } // if

  if (0 != _dbChange) { // Setup change of values, if provided.
    _dbChange->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbChange->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbChange->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-horiz",
				    "traction-shear-vert",
				    "traction-normal"};
	_dbChange->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Neumann.");
      } // switch
    _queryDB("change", _dbChange, spaceDim, pressureScale);

    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB("change time", _dbChange, 1, timeScale);
    _dbChange->close();

    if (0 != _dbTimeHistory)
      _dbTimeHistory->open();
  } // if

} // _queryDatabases

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::Neumann::_queryDB(const char* name,
			      spatialdata::spatialdb::SpatialDB* const db,
			      const int querySize,
			      const PylithScalar scale)
{ // _queryDB
  assert(name);
  assert(db);
  assert(_boundaryMesh);
  assert(_quadrature);
  assert(_parameters);

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  
  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  scalar_array valuesCell(numQuadPts*querySize);
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get sections.
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    subSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex(name);
  const int valueFiberDim = _parameters->sectionFiberDim(name);
  assert(valueIndex+valueFiberDim <= parametersFiberDim);
  scalar_array parametersCell(parametersFiberDim);

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(cs);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  // Compute quadrature information
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_boundaryMesh, cells);
#endif

  // Loop over cells in boundary mesh and perform queries.
  for (SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif
    const scalar_array& quadPtsNondim = _quadrature->quadPts();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    
    valuesCell = 0.0;
    for (int iQuad=0, iSpace=0; 
	 iQuad < numQuadPts;
	 ++iQuad, iSpace+=spaceDim) {
      const int err = db->query(&valuesCell[iQuad*querySize], querySize,
				&quadPtsGlobal[iSpace], spaceDim, cs);
      
      if (err) {
	std::ostringstream msg;
	msg << "Could not find values at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << " " << quadPtsGlobal[i+iSpace];
	msg << ") for traction boundary condition " << _label << "\n"
	    << "using spatial database " << db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if

    } // for
    _normalizer->nondimensionalize(&valuesCell[0], valuesCell.size(),
				   scale);

    // Update section
    assert(parametersFiberDim == parametersSection->getFiberDimension(*c_iter));
    parametersSection->restrictPoint(*c_iter, 
				     &parametersCell[0], parametersCell.size());
    for (int i=0; i < valueFiberDim; ++i)
      parametersCell[valueIndex+i] = valuesCell[i];
    
    parametersSection->updatePoint(*c_iter, &parametersCell[0]);
  } // for
} // _queryDB

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
  pylith::bc::Neumann::_paramsLocalToGlobal(const PylithScalar upDir[3])
{ // _paramsLocalToGlobal
  assert(_boundaryMesh);
  assert(_parameters);
  assert(_quadrature);

  scalar_array up(3);
  for (int i=0; i < 3; ++i)
    up[i] = upDir[i];

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Quadrature related values.
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  scalar_array quadPtRef(cellDim);
  const scalar_array& quadPtsRef = _quadrature->quadPtsRef();
  
  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int jacobianSize = spaceDim * cellDim;
  scalar_array jacobian(jacobianSize);
  PylithScalar jacobianDet = 0;
  scalar_array orientation(orientationSize);

  // Get sections.
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    subSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  const ALE::Obj<SubRealUniformSection>& parametersSection =
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  scalar_array parametersCellLocal(parametersFiberDim);
  scalar_array parametersCellGlobal(parametersFiberDim);

  const int initialIndex = 
    (_dbInitial) ? _parameters->sectionIndex("initial") : -1;
  const int initialFiberDim = 
    (_dbInitial) ? _parameters->sectionFiberDim("initial") : 0;

  const int rateIndex = (_dbRate) ? _parameters->sectionIndex("rate") : -1;
  const int rateFiberDim = 
    (_dbRate) ? _parameters->sectionFiberDim("rate") : 0;

  const int changeIndex = (_dbChange) ? _parameters->sectionIndex("change") : -1;
  const int changeFiberDim = 
    (_dbChange) ? _parameters->sectionFiberDim("change") : 0;

  // Loop over cells in boundary mesh, compute orientations, and then
  // rotate corresponding traction vector from local to global coordinates.
  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif
    // Get parameters
    parametersSection->restrictPoint(*c_iter, 
				     &parametersCellLocal[0], 
				     parametersCellLocal.size());
    // Copy parameter values to global (not all are tranformed)
    parametersCellGlobal = parametersCellLocal;

    for(int iQuad=0, iRef=0, iSpace=0; iQuad < numQuadPts;
	++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      // Compute Jacobian and determinant at quadrature point, then get
      // orientation.
      memcpy(&quadPtRef[0], &quadPtsRef[iRef], cellDim*sizeof(PylithScalar));
#if defined(PRECOMPUTE_GEOMETRY)
      coordsVisitor.clear();
      subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
#endif
      cellGeometry.jacobian(&jacobian, &jacobianDet,
			    coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      if (0 != _dbInitial) {
	// Rotate traction vector from local coordinate system to global
	// coordinate system
	assert(initialIndex >= 0);
	assert(initialFiberDim == numQuadPts*spaceDim);
	PylithScalar* initialGlobal = &parametersCellGlobal[initialIndex+iSpace];
	PylithScalar* initialLocal = &parametersCellLocal[initialIndex+iSpace];
	for(int iDim = 0; iDim < spaceDim; ++iDim) {
	  initialGlobal[iDim] = 0.0;
	  for(int jDim = 0; jDim < spaceDim; ++jDim)
	    initialGlobal[iDim] += 
	      orientation[jDim*spaceDim+iDim] * initialLocal[jDim];
	} // for
      } // if

      if (0 != _dbRate) {
	// Rotate traction vector from local coordinate system to global
	// coordinate system
	assert(rateIndex >= 0);
	assert(rateFiberDim == numQuadPts*spaceDim);
	PylithScalar* rateGlobal = &parametersCellGlobal[rateIndex+iSpace];
	PylithScalar* rateLocal = &parametersCellLocal[rateIndex+iSpace];
	for(int iDim = 0; iDim < spaceDim; ++iDim) {
	  rateGlobal[iDim] = 0.0;
	  for(int jDim = 0; jDim < spaceDim; ++jDim)
	    rateGlobal[iDim] +=
	      orientation[jDim*spaceDim+iDim] * rateLocal[jDim];
	} // for
      } // if

      if (0 != _dbChange) {
	// Rotate traction vector from local coordinate system to global
	// coordinate system
	assert(changeIndex >= 0);
	assert(changeFiberDim == numQuadPts*spaceDim);
	PylithScalar* changeGlobal = &parametersCellGlobal[changeIndex+iSpace];
	PylithScalar* changeLocal = &parametersCellLocal[changeIndex+iSpace];
	for(int iDim = 0; iDim < spaceDim; ++iDim) {
	  changeGlobal[iDim] = 0.0;
	  for(int jDim = 0; jDim < spaceDim; ++jDim)
	    changeGlobal[iDim] +=
	      orientation[jDim*spaceDim+iDim] * changeLocal[jDim];
	} // for
      } // if

    } // for
    
    // Update sections
    assert(parametersFiberDim == parametersSection->getFiberDimension(*c_iter));
    parametersSection->updatePoint(*c_iter, &parametersCellGlobal[0]);
  } // for
} // paramsLocalToGlobal

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of Submesh.
void
pylith::bc::Neumann::_calculateValue(const PylithScalar t)
{ // _calculateValue
  assert(_parameters);
  assert(_boundaryMesh);
  assert(_quadrature);

  const PylithScalar timeScale = _getNormalizer().timeScale();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();

  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  scalar_array parametersCell(parametersFiberDim);
  
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(numQuadPts*spaceDim == valueFiberDim);

  const int initialIndex = 
    (_dbInitial) ? _parameters->sectionIndex("initial") : -1;
  const int initialFiberDim = 
    (_dbInitial) ? _parameters->sectionFiberDim("initial") : 0;

  const int rateIndex = 
    (_dbRate) ? _parameters->sectionIndex("rate") : -1;
  const int rateFiberDim = 
    (_dbRate) ? _parameters->sectionFiberDim("rate") : 0;
  const int rateTimeIndex = 
    (_dbRate) ? _parameters->sectionIndex("rate time") : -1;
  const int rateTimeFiberDim = 
    (_dbRate) ? _parameters->sectionFiberDim("rate time") : 0;

  const int changeIndex = 
    (_dbChange) ? _parameters->sectionIndex("change") : -1;
  const int changeFiberDim = 
    (_dbChange) ? _parameters->sectionFiberDim("change") : 0;
  const int changeTimeIndex = 
    (_dbChange) ? _parameters->sectionIndex("change time") : -1;
  const int changeTimeFiberDim = 
    (_dbChange) ? _parameters->sectionFiberDim("change time") : 0;

  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    
    assert(parametersFiberDim == parametersSection->getFiberDimension(*c_iter));
    parametersSection->restrictPoint(*c_iter, 
				     &parametersCell[0], parametersCell.size());
    for (int i=0; i < valueFiberDim; ++i)
      parametersCell[valueIndex+i] = 0.0;

    // Contribution from initial value
    if (0 != _dbInitial) {
      assert(initialIndex >= 0);
      assert(initialFiberDim == valueFiberDim);
      for (int i=0; i < initialFiberDim; ++i)
	parametersCell[valueIndex+i] += parametersCell[initialIndex+i];
    } // if
    
    // Contribution from rate of change of value
    if (0 != _dbRate) {
      assert(rateIndex >= 0);
      assert(rateFiberDim == valueFiberDim);
      assert(rateTimeIndex >= 0);
      assert(rateTimeFiberDim == numQuadPts);
      
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const PylithScalar tRel = t - parametersCell[rateTimeIndex+iQuad];
	if (tRel > 0.0)  // rate of change integrated over time
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    parametersCell[valueIndex+iQuad*spaceDim+iDim] += 
	      parametersCell[rateIndex+iQuad*spaceDim+iDim] * tRel;
      } // for
    } // if
    
    // Contribution from change of value
    if (0 != _dbChange) {
      assert(changeIndex >= 0);
      assert(changeFiberDim == valueFiberDim);
      assert(changeTimeIndex >= 0);
      assert(changeTimeFiberDim == numQuadPts);

      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const PylithScalar tRel = t - parametersCell[changeTimeIndex+iQuad];
	if (tRel >= 0) { // change in value over time
	  PylithScalar scale = 1.0;
	  if (0 != _dbTimeHistory) {
	    PylithScalar tDim = tRel;
	    _getNormalizer().dimensionalize(&tDim, 1, timeScale);
	    const int err = _dbTimeHistory->query(&scale, tDim);
	    if (0 != err) {
	      std::ostringstream msg;
	      msg << "Error querying for time '" << tDim 
		  << "' in time history database "
		  << _dbTimeHistory->label() << ".";
	      throw std::runtime_error(msg.str());
	    } // if
	  } // if
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    parametersCell[valueIndex+iQuad*spaceDim+iDim] += 
	      parametersCell[changeIndex+iQuad*spaceDim+iDim] * scale;
	} // if
      } // for
    } // if
    
    parametersSection->updatePoint(*c_iter, &parametersCell[0]);
  } // for
}  // _calculateValue


// End of file 
