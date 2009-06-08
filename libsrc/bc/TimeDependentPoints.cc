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

#include "TimeDependentPoints.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
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
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::TimeDependentPoints::TimeDependentPoints(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::TimeDependentPoints::~TimeDependentPoints(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Query databases for parameters.
void
pylith::bc::TimeDependentPoints::_queryDatabases(const topology::Mesh& mesh,
						 const double valueScale,
						 const char* fieldName)
{ // _queryDatabases
  const double timeScale = _getNormalizer().timeScale();
  const double rateScale = valueScale / timeScale;

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  char** valueNames = (numBCDOF > 0) ? new char*[numBCDOF] : 0;
  char** rateNames = (numBCDOF > 0) ? new char*[numBCDOF] : 0;
  const std::string& valuePrefix = std::string(fieldName) + "-";
  const std::string& ratePrefix = std::string(fieldName) + "-rate-";
  for (int i=0; i < numBCDOF; ++i) {
    std::string suffix;
    switch (_bcDOF[i])
      { // switch
      case 0 :
	suffix = "x";
	break;
      case 1 :
	suffix = "y";
	break;
      case 2 :
	suffix = "z";
	break;
      } // switch
    std::string name = valuePrefix + suffix;
    int size = 1 + name.length();
    valueNames[i] = new char[size];
    strcpy(valueNames[i], name.c_str());

    name = ratePrefix + suffix;
    size = 1 + name.length();
    rateNames[i] = new char[size];
    strcpy(rateNames[i], name.c_str());
  } // for

  delete _parameters;
  _parameters = new topology::Fields<topology::Field<topology::Mesh> >(mesh);

  // Create section to hold time dependent values
  _parameters->add("value", fieldName);
  topology::Field<topology::Mesh>& value = _parameters->get("value");
  value.scale(valueScale);
  value.vectorFieldType(topology::FieldBase::OTHER);
  value.newSection(_points, numBCDOF);
  value.allocate();

  if (0 != _dbInitial) { // Setup initial values, if provided.
    std::string fieldLabel = std::string("initial_") + std::string(fieldName);
    _parameters->add("initial", fieldLabel.c_str());
    topology::Field<topology::Mesh>& initial = 
      _parameters->get("initial");
    initial.cloneSection(value);
    initial.scale(valueScale);
    initial.vectorFieldType(topology::FieldBase::OTHER);

    _dbInitial->open();
    _dbInitial->queryVals(valueNames, numBCDOF);
    _queryDB(&initial, _dbInitial, numBCDOF, valueScale);
    _dbInitial->close();
  } // if

  if (0 != _dbRate) { // Setup rate of change of values, if provided.
    std::string fieldLabel = std::string("rate_") + std::string(fieldName);
    _parameters->add("rate", fieldLabel.c_str());
    topology::Field<topology::Mesh>& rate = 
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
    topology::Field<topology::Mesh>& rateTime = 
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
    topology::Field<topology::Mesh>& change = 
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
    topology::Field<topology::Mesh>& changeTime = 
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
pylith::bc::TimeDependentPoints::_queryDB(topology::Field<topology::Mesh>* field,
				 spatialdata::spatialdb::SpatialDB* const db,
				 const int querySize,
				 const double scale)
{ // _queryDB
  assert(0 != field);
  assert(0 != db);

  const topology::Mesh& mesh = field->mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const double lengthScale = _getNormalizer().lengthScale();

  double_array coordsVertex(spaceDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  const ALE::Obj<RealSection>& section = field->section();
  assert(!section.isNull());

  double_array valuesVertex(querySize);

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    // Get dimensionalized coordinates of vertex
    coordinates->restrictPoint(_points[iPoint], 
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
    section->updatePoint(_points[iPoint], &valuesVertex[0]);
  } // for
} // _queryDB

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of points.
void
pylith::bc::TimeDependentPoints::_calculateValue(const double t)
{ // _calculateValue
  assert(0 != _parameters);

  const ALE::Obj<RealSection>& valueSection = 
    _parameters->get("value").section();
  assert(!valueSection.isNull());
  valueSection->zero();

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  const double timeScale = _getNormalizer().timeScale();

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
	  _getNormalizer().dimensionalize(&tDim, 1, timeScale);
	  const int err = _dbTimeHistory->query(&scale, tDim);
	  if (0 != err) {
	    std::ostringstream msg;
	    msg << "Error querying for time '" << tDim << "' in time history database "
		<< _dbTimeHistory->label() << ".";
	    throw std::runtime_error(msg.str());
	  } // if
	} // if
	values *= scale;
	valueSection->updateAddPoint(p_bc, &values[0]);
      } // if
    } // for
  } // if
}  // _calculateValue


// End of file 
