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

#include "PointForce.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::PointForce::PointForce(void) :
  _parameters(0),
  _dbInitial(0),
  _dbRate(0),
  _dbRateTime(0),
  _dbChange(0),
  _dbChangeTime(0),
  _dbTimeHistory(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::PointForce::~PointForce(void)
{ // destructor
  if (0 != _dbTimeHistory)
    _dbTimeHistory->close();

  delete _parameters; _parameters = 0;

  _dbInitial = 0; // TODO: Use shared pointers
  _dbRate = 0; // TODO: Use shared pointers
  _dbRateTime = 0; // TODO: Use shared pointers
  _dbChange = 0; // TODO: Use shared pointers
  _dbChangeTime = 0; // TODO: Use shared pointers
  _dbTimeHistory = 0; // TODO: Use shared pointers
} // destructor

// ----------------------------------------------------------------------
// Set indices of vertices with point forces.
void
pylith::bc::PointForce::forceDOF(const int* flags,
				 const int size)
{ // forceDOF
  if (size > 0)
    assert(0 != flags);

  _bcDOF.resize(size);
  for (int i=0; i < size; ++i)
    _bcDOF[i] = flags[i];
} // forceDOF

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::PointForce::initialize(const topology::Mesh& mesh,
				    const double upDir[3])
{ // initialize
  assert(0 != _normalizer);

  if (0 == _bcDOF.size())
    return;

  _getPoints(mesh);

  const double lengthScale = _normalizer->lengthScale();
  const double pressureScale = _normalizer->pressureScale();
  const double forceScale = pressureScale * lengthScale * lengthScale;

  _queryDatabases(mesh, forceScale, "force");
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator that
// do not require assembly over cells, vertices, or processors.
void
pylith::bc::PointForce::integrateResidualAssembled(
			   topology::Field<topology::Mesh>* residual,
			   const double t,
			   topology::SolutionFields* const fields)
{ // integrateResidualAssembled
  assert(0 != residual);
  assert(0 != _parameters);
  assert(0 != _normalizer);

  // Calculate spatial and temporal variation of value for BC.
  _calculateValue(t);

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();

  const topology::Mesh& mesh = residual->mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  double_array residualVertex(spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual->section();
  assert(!residualSection.isNull());

  double_array valuesVertex(numBCDOF);
  const ALE::Obj<RealSection>& valueSection = 
    _parameters->get("value").section();
  assert(!valueSection.isNull());
  
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label.
    residualVertex *= 0.0; // Reset residual contribution to zero.
    
    valueSection->restrictPoint(p_bc, &valuesVertex[0], valuesVertex.size());
    for (int iDOF=0; iDOF < numBCDOF; ++iDOF)
      residualVertex[_bcDOF[iDOF]] += valuesVertex[iDOF];
    residualSection->updateAddPoint(p_bc, &residualVertex[0]);
  } // for

} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of points.
void
pylith::bc::PointForce::_calculateValue(const double t)
{ // _calculateValue
  assert(0 != _parameters);
  assert(0 != _normalizer);

  const ALE::Obj<RealSection>& valueSection = 
    _parameters->get("value").section();
  valueSection->zero();

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  const double timeScale = _normalizer->timeScale();

  const ALE::Obj<RealSection>& changeSection = (0 != _dbChange) ?
    _parameters->get("change").section() : 0;
  const ALE::Obj<RealSection>& changeTimeSection = (0 != _dbChange) ?
    _parameters->get("change time").section() : 0;

  // Contribution from initial value
  if (0 != _dbInitial) {
    double_array values(numBCDOF);
    const ALE::Obj<RealSection>& initialSection = 
      _parameters->get("initial").section();
    assert(!initialSection.isNull());

    for (int iPoint=0; iPoint < numPoints; ++iPoint) {
      const int p_bc = _points[iPoint]; // Get point label

      initialSection->restrictPoint(p_bc, &values[0], values.size());
      valueSection->updateAddPoint(p_bc, &values[0]);
    } // for
  } // if
  
  
  // Contribution from rate of change of value
  if (0 != _dbRate) {
    double_array values(numBCDOF);
    double tRate = 0.0;
    const ALE::Obj<RealSection>& rateSection =
      _parameters->get("rate").section();
    assert(!rateSection.isNull());
    const ALE::Obj<RealSection>& rateTimeSection =
      _parameters->get("rate time").section();
    assert(!rateTimeSection.isNull());

    for (int iPoint=0; iPoint < numPoints; ++iPoint) {
      const int p_bc = _points[iPoint]; // Get point label

      rateSection->restrictPoint(p_bc, &values[0], values.size());
      rateTimeSection->restrictPoint(p_bc, &tRate, 1);
      if (t > tRate) { // rate of change integrated over time
	values *= (t - tRate);
	valueSection->updateAddPoint(p_bc, &values[0]);
      } // if
    } // for
  } // if
      
  // Contribution from change of value
  if (0 != _dbChange) {
    double_array values(numBCDOF);
    double tChange = 0.0;
    const ALE::Obj<RealSection>& changeSection =
      _parameters->get("change").section();
    assert(!changeSection.isNull());
    const ALE::Obj<RealSection>& changeTimeSection =
      _parameters->get("change time").section();
    assert(!changeTimeSection.isNull());

    for (int iPoint=0; iPoint < numPoints; ++iPoint) {
      const int p_bc = _points[iPoint]; // Get point label

      changeSection->restrictPoint(p_bc, &values[0], values.size());
      changeTimeSection->restrictPoint(p_bc, &tChange, 1);
      if (t > tChange) { // change in value over time
	double scale = 1.0;
	if (0 != _dbTimeHistory) {
	  double tDim = t - tChange;
	  _normalizer->dimensionalize(&tDim, 1, timeScale);
	  const int err = _dbTimeHistory->query(&scale, tDim);
	  if (0 != err) {
	    std::ostringstream msg;
	    msg << "ADD SOMETHING HERE";
	    throw std::runtime_error(msg.str());
	  } // if
	} // if
	values *= scale;
	valueSection->updateAddPoint(p_bc, &values[0]);
      } // if
    } // for
  } // if
}  // _calculateValue

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::PointForce::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  if (0 != _dbRate && 0 == _dbRateTime) {
    std::ostringstream msg;
    msg << "Point force boundary condition '" << label() << "',\n has a rate "
	<< "of change spatial database but no rate of change start time "
	<< "spatial database.";
    throw std::runtime_error(msg.str());
  } // if
  if (0 == _dbRate && 0 != _dbRateTime) {
    std::ostringstream msg;
    msg << "Point force boundary condition '" << label() << "',\n has a rate "
	<< "of change start time spatial database but no rate of change "
	<< "spatial database.";
    throw std::runtime_error(msg.str());
  } // if

  if (0 != _dbChange && 0 == _dbChangeTime) {
    std::ostringstream msg;
    msg << "Point force boundary condition '" << label() << "',\n has a "
	<< "change in value spatial database but change in value start time "
	<< "spatial database.";
    throw std::runtime_error(msg.str());
  } // if
  if (0 == _dbChange && 0 != _dbChangeTime) {
    std::ostringstream msg;
    msg << "Point force boundary condition '" << label() << "',\n has a "
	<< "change in value start time spatial database but change in value "
	<< "spatial database.";
    throw std::runtime_error(msg.str());
  } // if
  if (0 == _dbChange && 0 != _dbTimeHistory) {
    std::ostringstream msg;
    msg << "Point force boundary condition '" << label() << "',\n has a "
	<< "time history database but not change in value spatial database.";
    throw std::runtime_error(msg.str());
} // if
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get mesh labels for points associated with point forces.
void
pylith::bc::PointForce::_getPoints(const topology::Mesh& mesh)
{ // _getPoints
  typedef topology::Mesh::IntSection::chart_type chart_type;

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::int_section_type>& groupField = 
    sieveMesh->getIntSection(_label);
  if (groupField.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << _label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  assert(!groupField.isNull());
  const chart_type& chart = groupField->getChart();
  const chart_type::const_iterator& chartEnd = chart.end();
  const int numPoints = groupField->size();
  _points.resize(numPoints);
  int i = 0;
  for(chart_type::const_iterator c_iter = chart.begin();
      c_iter != chartEnd;
      ++c_iter)
    if (groupField->getFiberDimension(*c_iter))
      _points[i++] = *c_iter;
} // _getPoints

// ----------------------------------------------------------------------
// Query databases for parameters.
void
pylith::bc::PointForce::_queryDatabases(const topology::Mesh& mesh,
					const double valueScale,
					const char* fieldName)
{ // _queryDatabases
  assert(0 != _normalizer);
  const double timeScale = _normalizer->timeScale();
  const double rateScale = valueScale / timeScale;

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  char** valueNames = (numBCDOF > 0) ? new char*[numBCDOF] : 0;
  for (int i=0; i < numBCDOF; ++i) {
    std::ostringstream name;
    name << "dof-" << _bcDOF[i];
    const int size = 1 + name.str().length();
    valueNames[i] = new char[size];
    strcpy(valueNames[i], name.str().c_str());
  } // for

  delete _parameters;
  _parameters = new topology::Fields<topology::Field<topology::Mesh> >(mesh);

  if (0 != _dbInitial) { // Setup initial values, if provided.
    std::string fieldLabel = std::string("initial_") + std::string(fieldName);
    _parameters->add("initial", fieldLabel.c_str());
    topology::Field<topology::Mesh>& initial = 
      _parameters->get("initial");
    initial.newSection(_points, numBCDOF);
    initial.allocate();
    initial.scale(valueScale);
    initial.vectorFieldType(topology::FieldBase::OTHER);

    _dbInitial->queryVals(valueNames, numBCDOF);
    _queryDB(&initial, _dbInitial, numBCDOF, valueScale);
  } // if

  if (0 != _dbRate) { // Setup rate of change of values, if provided.
    std::string fieldLabel = std::string("rate_") + std::string(fieldName);
    _parameters->add("rate", fieldLabel.c_str());
    topology::Field<topology::Mesh>& rate = 
      _parameters->get("rate");
    rate.newSection(_points, numBCDOF);
    rate.allocate();
    rate.scale(valueScale);
    rate.vectorFieldType(topology::FieldBase::VECTOR);
    const ALE::Obj<RealSection>& rateSection = rate.section();
    assert(!rateSection.isNull());

    _dbRate->queryVals(valueNames, numBCDOF);
    _queryDB(&rate, _dbRate, numBCDOF, rateScale);

    std::string timeLabel = 
      std::string("rate_time_") + std::string(fieldName);
    _parameters->add("rate time", timeLabel.c_str());
    topology::Field<topology::Mesh>& rateTime = 
      _parameters->get("rate time");
    rateTime.newSection(rateSection->getChart(), 1);
    rateTime.allocate();
    rateTime.scale(timeScale);
    rateTime.vectorFieldType(topology::FieldBase::SCALAR);

    const char* timeNames[1] = { "rate-start-time" };
    _dbRateTime->queryVals(timeNames, 1);
    _queryDB(&rateTime, _dbRateTime, 1, timeScale);
  } // if

  if (0 != _dbChange) { // Setup change of values, if provided.
    // ADD STUFF HERE

    if (0 != _dbTimeHistory)
      _dbTimeHistory->open();
  } // if

  _parameters->add("value", fieldName);
  topology::Field<topology::Mesh>& value = _parameters->get("value");
  value.scale(valueScale);
  value.vectorFieldType(topology::FieldBase::OTHER);
} // _queryDatabases


// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::PointForce::_queryDB(topology::Field<topology::Mesh>* field,
				 spatialdata::spatialdb::SpatialDB* const db,
				 const int querySize,
				 const double scale)
{ // _queryDB
  assert(0 != field);
  assert(0 != db);
  assert(0 != _normalizer);

  const topology::Mesh& mesh = field->mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();

  double_array coordsVertex(spaceDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  const ALE::Obj<RealSection>& section = field->section();
  assert(!section.isNull());

  double_array valuesVertex(querySize);

  db->open();
  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    // Get dimensionalized coordinates of vertex
    coordinates->restrictPoint(_points[iPoint], 
			       &coordsVertex[0], coordsVertex.size());
    _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(),
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
    _normalizer->nondimensionalize(&valuesVertex[0], valuesVertex.size(),
				   scale);
    section->updatePoint(_points[iPoint], &valuesVertex[0]);
  } // for

  db->close();
} // _queryDB


// End of file 
