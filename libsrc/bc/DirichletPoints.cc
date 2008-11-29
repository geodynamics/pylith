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

#include "DirichletPoints.hh" // implementation of object methods

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <string.h> // USES strcpy()
#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletPoints::DirichletPoints(void) :
  _tRef(0.0),
  _dbRate(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletPoints::~DirichletPoints(void)
{ // destructor
  _dbRate = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletPoints::initialize(
				   const ALE::Obj<Mesh>& mesh,
				   const spatialdata::geocoords::CoordSys* cs,
				   const double_array& upDir)
{ // initialize
  assert(0 != _db);
  assert(0 != _dbRate);
  assert(!mesh.isNull());
  assert(0 != cs);

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  // Get points associated with boundary condition
  const ALE::Obj<int_section_type>& groupField = mesh->getIntSection(_label);
  if (groupField.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << _label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  assert(!groupField.isNull());
  const int_section_type::chart_type& chart = groupField->getChart();
  const int numPoints = groupField->size();
  _points.resize(numPoints);
  int i = 0;
  for(int_section_type::chart_type::const_iterator c_iter = chart.begin();
      c_iter != chart.end();
      ++c_iter) {
    if (groupField->getFiberDimension(*c_iter)) _points[i++] = *c_iter;
  }

  // Get values for degrees of freedom
  char** valueNames = (numFixedDOF > 0) ? new char*[numFixedDOF] : 0;
  for (int i=0; i < numFixedDOF; ++i) {
    std::ostringstream name;
    name << "dof-" << _fixedDOF[i];
    const int size = 1 + name.str().length();
    valueNames[i] = new char[size];
    strcpy(valueNames[i], name.str().c_str());
  } // for
  _db->open();
  _db->queryVals((const char**) valueNames, numFixedDOF);
  _dbRate->open();
  _dbRate->queryVals((const char**) valueNames, numFixedDOF);
  for (int i=0; i < numFixedDOF; ++i) {
    delete[] valueNames[i]; valueNames[i] = 0;
  } // for
  delete[] valueNames; valueNames = 0;

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const int spaceDim = cs->spaceDim();

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double velocityScale = 
    _normalizer->lengthScale() / _normalizer->timeScale();

  _valuesInitial.resize(numPoints*numFixedDOF);
  _valuesRate.resize(numPoints*numFixedDOF);
  double_array queryValues(numFixedDOF);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    // Get coordinates of vertex
    const real_section_type::value_type* vCoords = 
      coordinates->restrictPoint(_points[iPoint]);
    int err = _db->query(&queryValues[0], numFixedDOF, vCoords, 
				spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find values at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << _db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      _valuesInitial[numFixedDOF*iPoint+iDOF] = 
	_normalizer->nondimensionalize(queryValues[iDOF], lengthScale);

    err = _dbRate->query(&queryValues[0], numFixedDOF, vCoords, 
			 spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find values at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << _dbRate->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      _valuesRate[numFixedDOF*iPoint+iDOF] =
	_normalizer->nondimensionalize(queryValues[iDOF], velocityScale);
  } // for
  _db->close();
  _dbRate->close();
} // initialize

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::DirichletPoints::setConstraintSizes(
				     const ALE::Obj<real_section_type>& field,
				     const ALE::Obj<Mesh>& mesh)
{ // setConstraintSizes
  assert(!field.isNull());
  assert(!mesh.isNull());

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const int numPoints = _points.size();
  _offsetLocal.resize(numPoints);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int fiberDim = field->getFiberDimension(_points[iPoint]);
    const int curNumConstraints = field->getConstraintDimension(_points[iPoint]);
    if (curNumConstraints + numFixedDOF > fiberDim) {
      std::ostringstream msg;
      msg
	<< "Found overly constrained point while setting up constraints for\n"
	<< "DirichletPoints boundary condition '" << _label << "'.\n"
	<< "Number of DOF at point " << _points[iPoint] << " is " << fiberDim
	<< "\nand number of attempted constraints is "
	<< curNumConstraints+numFixedDOF << ".";
      throw std::runtime_error(msg.str());
    } // if
    _offsetLocal[iPoint] = curNumConstraints;
    field->addConstraintDimension(_points[iPoint], numFixedDOF);
  } // for
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::bc::DirichletPoints::setConstraints(
				    const ALE::Obj<real_section_type>& field,
				    const ALE::Obj<Mesh>& mesh)
{ // setConstraints
  assert(!field.isNull());
  assert(!mesh.isNull());

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const Mesh::point_type point = _points[iPoint];

    // Get list of currently constrained DOF
    const int* curFixedDOF = field->getConstraintDof(point);
    const int numTotalConstrained = field->getConstraintDimension(point);

    // Create array holding all constrained DOF
    int_array allFixedDOF(curFixedDOF, numTotalConstrained);

    // Verify other BC has not already constrained DOF
    const int numPrevious = _offsetLocal[iPoint];
    for (int iDOF=0; iDOF < numPrevious; ++iDOF)
      for (int jDOF=0; jDOF < numFixedDOF; ++jDOF)
	if (allFixedDOF[iDOF] == _fixedDOF[jDOF]) {
	  std::ostringstream msg;
	  msg << "Found multiple constraints on degrees of freedom at\n"
	      << "point while setting up constraints for DirichletPoints\n"
	      << "boundary condition '" << _label << "'.\n"
	      << "Degree of freedom " << _fixedDOF[jDOF] 
	      << " is already constrained by another Dirichlet BC.";
	  throw std::runtime_error(msg.str());
	} // if

    // Add in the ones for this DirichletPoints BC
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      assert(_offsetLocal[iPoint]+iDOF < numTotalConstrained);
      allFixedDOF[_offsetLocal[iPoint]+iDOF] = _fixedDOF[iDOF];
    } // for

    // Fill in rest of values not yet set (will be set by
    // another DirichletPoints BC)
    for (int iDOF=_offsetLocal[iPoint]+numFixedDOF; 
	 iDOF < numTotalConstrained; 
	 ++iDOF) {
      assert(iDOF < numTotalConstrained);
      allFixedDOF[iDOF] = 999;
    } // for

    // Sort list of constrained DOF
    //   I need these sorted for my update algorithms to work properly
    std::sort(&allFixedDOF[0], &allFixedDOF[numTotalConstrained]);

    // Update list of constrained DOF
    field->setConstraintDof(point, &allFixedDOF[0]);
  } // for
} // setConstraints

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::DirichletPoints::setField(const double t,
				      const ALE::Obj<real_section_type>& field,
				      const ALE::Obj<Mesh>& mesh)
{ // setField
  assert(!field.isNull());
  assert(!mesh.isNull());

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const int numPoints = _points.size();
  const int fiberDimension = 
    (numPoints > 0) ? field->getFiberDimension(_points[0]) : 0;
  double_array allValues(fiberDimension);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const Mesh::point_type point = _points[iPoint];
    assert(fiberDimension == field->getFiberDimension(point));
    mesh->restrictClosure(field, point, &allValues[0], fiberDimension);
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      allValues[_fixedDOF[iDOF]] = _valuesInitial[iPoint*numFixedDOF+iDOF];
    if (t > _tRef)
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
	allValues[_fixedDOF[iDOF]] += 
	  (t-_tRef) * _valuesRate[iPoint*numFixedDOF+iDOF];
    field->updatePointAll(_points[iPoint], &allValues[0]);
  } // for
} // setField


// End of file 
