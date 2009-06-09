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

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann::Neumann(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::Neumann::initialize(const topology::Mesh& mesh,
				const double upDir[3])
{ // initialize
  double_array up(upDir, 3);
  _queryDatabases(up);
  _paramsLocalToGlobal();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Create section for traction vector in global coordinates
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const int fiberDim = spaceDim * numQuadPts;
  
  _parameters =
    new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
  assert(0 != _parameters);
  _parameters->add("traction", "traction");
  topology::Field<topology::SubMesh>& traction = _parameters->get("traction");
  traction.newSection(cells, fiberDim);
  traction.allocate();

  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int jacobianSize = spaceDim * cellDim;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientation(orientationSize);

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
  double_array coordinatesCell(numCorners*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    subSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  const ALE::Obj<SubRealSection>& tractionSection =
    _parameters->get("traction").section();
  assert(!tractionSection.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double pressureScale = _normalizer->pressureScale();

  // Compute quadrature information
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_boundaryMesh, cells);
#endif

  // Loop over cells in boundary mesh, compute orientations, and then
  // compute corresponding traction vector in global coordinates
  // (store values in _tractionGlobal).
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
    const double_array& quadPtsNondim = _quadrature->quadPts();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    
    cellTractionsGlobal = 0.0;
    for(int iQuad=0, iRef=0, iSpace=0; iQuad < numQuadPts;
	++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      // Get traction vector in local coordinate system at quadrature point
      const int err = _db->query(&tractionDataLocal[0], spaceDim,
				 &quadPtsGlobal[iSpace], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find traction values at (";
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
#if defined(PRECOMPUTE_GEOMETRY)
      coordsVisitor.clear();
      subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
#endif
      cellGeometry.jacobian(&jacobian, &jacobianDet,
			    coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
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
    tractionSection->updatePoint(*c_iter, &cellTractionsGlobal[0]);
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::Neumann::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _parameters);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();
  double_array tractionsCell(numQuadPts*spaceDim);

  // Get cell information
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<SubRealSection>& tractionSection =
    _parameters->get("traction").section();
  assert(!tractionSection.isNull());
  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::SubMesh::UpdateAddVisitor residualVisitor(*residualSection,
						      &_cellVector[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
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
    tractionSection->restrictPoint(*c_iter, 
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
    residualVisitor.clear();
    subSieveMesh->updateAdd(*c_iter, residualVisitor);

    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::Neumann::integrateJacobian(topology::Jacobian* jacobian,
 				       const double t,
 				       topology::SolutionFields* const fields)
{ // integrateJacobian
  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BCIntegratorSubMesh::verifyConfiguration(mesh);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get cell field for tractions.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::Neumann::cellField(const char* name,
			       topology::SolutionFields* const fields)
{ // cellField
  assert(0 != _parameters);
  assert(0 != name);

  if (0 == strcasecmp(name, "tractions")) {
    return _parameters->get("traction");
  } else {
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
  const double timeScale = _normalizer().timeScale();
  const double rateScale = valueScale / timeScale;

  const int numBCDOF = _bcDOF.size();
  char** valueNames = (numBCDOF > 0) ? new char*[numBCDOF] : 0;
  char** rateNames = (numBCDOF > 0) ? new char*[numBCDOF] : 0;
  const std::string& valuePrefix = std::string(fieldName) + "-";
  const std::string& ratePrefix = std::string(fieldName) + "-rate-";
  const string_vector& components = _valueComponents();
  for (int i=0; i < numBCDOF; ++i) {
    std::string name = valuePrefix + components[_bcDOF[i]];
    int size = 1 + name.length();
    valueNames[i] = new char[size];
    strcpy(valueNames[i], name.c_str());

    name = ratePrefix + components[_bcDOF[i]];
    size = 1 + name.length();
    rateNames[i] = new char[size];
    strcpy(rateNames[i], name.c_str());
  } // for

  delete _parameters; 
  _parameters = 
    new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);

  // Create section to hold time dependent values
  _parameters->add("value", fieldName);
  topology::Field<topology::SubMesh>& value = _parameters->get("value");
  value.scale(valueScale);
  value.vectorFieldType(topology::FieldBase::OTHER);
  value.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  value.allocate();

  if (0 != _dbInitial) { // Setup initial values, if provided.
    std::string fieldLabel = std::string("initial_") + std::string(fieldName);
    _parameters->add("initial", fieldLabel.c_str());
    topology::Field<topology::SubMesh>& initial = 
      _parameters->get("initial");
    initial.cloneSection(value);
    initial.scale(pressureScale);
    initial.vectorFieldType(topology::FieldBase::OTHER);

    _dbInitial->open();
    _dbInitial->queryVals(valueNames, spaceDim);
    _queryDB(&initial, _dbInitial, spaceDim, pressureScale);
    _dbInitial->close();
  } // if

  if (0 != _dbRate) { // Setup rate of change of values, if provided.
    std::string fieldLabel = std::string("rate_") + std::string(fieldName);
    _parameters->add("rate", fieldLabel.c_str());
    topology::Field<topology::SubMesh>& rate = 
      _parameters->get("rate");
    rate.cloneSection(value);
    rate.scale(rateScale);
    rate.vectorFieldType(topology::FieldBase::OTHER);
    const ALE::Obj<RealSection>& rateSection = rate.section();
    assert(!rateSection.isNull());

    _dbRate->open();
    _dbRate->queryVals(rateNames, numBCDOF);
    _queryDB(&rate, _dbRate, numBCDOF, rateScale);

    std::string timeLabel = 
      std::string("rate_time_") + std::string(fieldName);
    _parameters->add("rate time", timeLabel.c_str());
    topology::Field<topology::SubMesh>& rateTime = 
      _parameters->get("rate time");
    rateTime.newSection(rate, 1);
    rateTime.allocate();
    rateTime.scale(timeScale);
    rateTime.vectorFieldType(topology::FieldBase::SCALAR);

    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB(&rateTime, _dbRate, 1, timeScale);
    _dbRate->close();
  } // if

  if (0 != _dbChange) { // Setup change of values, if provided.
    std::string fieldLabel = std::string("change_") + std::string(fieldName);
    _parameters->add("change", fieldLabel.c_str());
    topology::Field<topology::SubMesh>& change = 
      _parameters->get("change");
    change.cloneSection(value);
    change.scale(valueScale);
    change.vectorFieldType(topology::FieldBase::OTHER);
    const ALE::Obj<RealSection>& changeSection = change.section();
    assert(!changeSection.isNull());

    _dbChange->open();
    _dbChange->queryVals(valueNames, numBCDOF);
    _queryDB(&change, _dbChange, numBCDOF, valueScale);

    std::string timeLabel = 
      std::string("change_time_") + std::string(fieldName);
    _parameters->add("change time", timeLabel.c_str());
    topology::Field<topology::SubMesh>& changeTime = 
      _parameters->get("change time");
    changeTime.newSection(change, 1);
    changeTime.allocate();
    changeTime.scale(timeScale);
    changeTime.vectorFieldType(topology::FieldBase::SCALAR);

    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB(&changeTime, _dbChange, 1, timeScale);
    _dbChange->close();

    if (0 != _dbTimeHistory)
      _dbTimeHistory->open();
  } // if

} // _queryDatabases

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::Neumann::_queryDB(topology::Field<topology::SubMesh>* field,
			      spatialdata::spatialdb::SpatialDB* const db,
			      const int querySize,
			      const double scale)
{ // _queryDB
  assert(0 != field);
  assert(0 != db);
  assert(0 != _boundaryMesh);

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const double lengthScale = _getNormalizer().lengthScale();

  double_array coordsVertex(spaceDim);
  const ALE::Obj<SieveMesh>& sieveMesh = _boundaryMesh->sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  const ALE::Obj<RealSection>& section = field->section();
  assert(!section.isNull());

  double_array valuesCell(querySize);

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Get dimensionalized coordinates of vertex
    coordinates->restrictPoint(_Submesh[iPoint], 
			       &coordsVertex[0], coordsVertex.size());
    _getNormalizer().dimensionalize(&coordsVertex[0], coordsVertex.size(),
				lengthScale);
    int err = db->query(&valuesVertex[0], valuesVertex.size(), 
			&coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Error querying for '" << field->label() << "' at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") using spatial database " << db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    _getNormalizer().nondimensionalize(&valuesVertex[0], valuesVertex.size(),
				   scale);
    section->updatePoint(_Submesh[iPoint], &valuesVertex[0]);
  } // for
} // _queryDB

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of Submesh.
void
pylith::bc::Neumann::_calculateValue(const double t)
{ // _calculateValue
  assert(0 != _parameters);

  const ALE::Obj<RealSection>& valueSection = 
    _parameters->get("value").section();
  assert(!valueSection.isNull());
  valueSection->zero();

  const int numSubmesh = _Submesh.size();
  const int numBCDOF = _bcDOF.size();
  const double timeScale = _getNormalizer().timeScale();

  const ALE::Obj<RealSection>& initialSection = (0 != _dbInitial) ?
    _parameters->get("initial").section() : 0;
  const ALE::Obj<RealSection>& rateSection = ( 0 != _dbRate) ?
    _parameters->get("rate").section() : 0;
  const ALE::Obj<RealSection>& rateTimeSection = (0 != _dbRate) ?
    _parameters->get("rate time").section() : 0;
  const ALE::Obj<RealSection>& changeSection = ( 0 != _dbChange) ?
    _parameters->get("change").section() : 0;
  const ALE::Obj<RealSection>& changeTimeSection = ( 0 != _dbChange) ?
    _parameters->get("change time").section() : 0;

  double_array valuesVertex(numBCDOF);
  double_array bufferVertex(numBCDOF);
  for (int iPoint=0; iPoint < numSubmesh; ++iPoint) {
    const int p_bc = _Submesh[iPoint]; // Get point label
    
    valuesVertex = 0.0;
    
    // Contribution from initial value
    if (0 != _dbInitial) {
      assert(!initialSection.isNull());
      initialSection->restrictPoint(p_bc, 
				    &bufferVertex[0], bufferVertex.size());
      valuesVertex += bufferVertex;
    } // if
    
    // Contribution from rate of change of value
    if (0 != _dbRate) {
      assert(!rateSection.isNull());
      assert(!rateTimeSection.isNull());
      double tRate = 0.0;
      
      rateSection->restrictPoint(p_bc, &bufferVertex[0], bufferVertex.size());
      rateTimeSection->restrictPoint(p_bc, &tRate, 1);
      if (t > tRate) { // rate of change integrated over time
	bufferVertex *= (t - tRate);
	valuesVertex += bufferVertex;
      } // if
    } // if
    
    // Contribution from change of value
    if (0 != _dbChange) {
      assert(!changeSection.isNull());
      assert(!changeTimeSection.isNull());
      double tChange = 0.0;

      changeSection->restrictPoint(p_bc, &bufferVertex[0], bufferVertex.size());
      changeTimeSection->restrictPoint(p_bc, &tChange, 1);
      if (t >= tChange) { // change in value over time
	double scale = 1.0;
	if (0 != _dbTimeHistory) {
	  double tDim = t - tChange;
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
	bufferVertex *= scale;
	valuesVertex += bufferVertex;
      } // if
    } // if

    valueSection->updateAddPoint(p_bc, &valuesVertex[0]);
  } // for
}  // _calculateValue


// End of file 
