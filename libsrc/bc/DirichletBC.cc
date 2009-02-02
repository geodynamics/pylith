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

#include "DirichletBC.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletBC::DirichletBC(void) :
  _tRef(0.0),
  _dbRate(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBC::~DirichletBC(void)
{ // destructor
  _dbRate = 0; // Don't manage memory
} // destructor

// ----------------------------------------------------------------------
// Set indices of fixed degrees of freedom. 
void
pylith::bc::DirichletBC::fixedDOF(const int* flags,
				  const int size)
{ // fixedDOF
  if (size > 0)
    assert(0 != flags);

  _fixedDOF.resize(size);
  for (int i=0; i < size; ++i)
    _fixedDOF[i] = flags[i];
} // fixedDOF

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBC::initialize(const topology::Mesh& mesh,
				    const double upDir[3])
{ // initialize
  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  _getPoints(mesh);
  _setupQueryDatabases();
  _queryDatabases(mesh);
} // initialize

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraintSizes(const topology::Field& field)
{ // setConstraintSizes
  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<MeshRealSection>& section = field.section();
  assert(!section.isNull());

  const int numPoints = _points.size();
  _offsetLocal.resize(numPoints);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int fiberDim = section->getFiberDimension(_points[iPoint]);
    const int curNumConstraints = 
      section->getConstraintDimension(_points[iPoint]);
    if (curNumConstraints + numFixedDOF > fiberDim) {
      std::ostringstream msg;
      msg
	<< "Found overly constrained point while setting up constraints for\n"
	<< "DirichletBC boundary condition '" << _label << "'.\n"
	<< "Number of DOF at point " << _points[iPoint] << " is " << fiberDim
	<< "\nand number of attempted constraints is "
	<< curNumConstraints+numFixedDOF << ".";
      throw std::runtime_error(msg.str());
    } // if
    _offsetLocal[iPoint] = curNumConstraints;
    section->addConstraintDimension(_points[iPoint], numFixedDOF);
  } // for
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraints(const topology::Field& field)
{ // setConstraints
  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<MeshRealSection>& section = field.section();
  assert(!section.isNull());

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];

    // Get list of currently constrained DOF
    const int* curFixedDOF = section->getConstraintDof(point);
    const int numTotalConstrained = section->getConstraintDimension(point);

    // Create array holding all constrained DOF
    int_array allFixedDOF(curFixedDOF, numTotalConstrained);

    // Verify other BC has not already constrained DOF
    const int numPrevious = _offsetLocal[iPoint];
    for (int iDOF=0; iDOF < numPrevious; ++iDOF)
      for (int jDOF=0; jDOF < numFixedDOF; ++jDOF)
	if (allFixedDOF[iDOF] == _fixedDOF[jDOF]) {
	  std::ostringstream msg;
	  msg << "Found multiple constraints on degrees of freedom at\n"
	      << "point while setting up constraints for DirichletBC\n"
	      << "boundary condition '" << _label << "'.\n"
	      << "Degree of freedom " << _fixedDOF[jDOF] 
	      << " is already constrained by another Dirichlet BC.";
	  throw std::runtime_error(msg.str());
	} // if

    // Add in the ones for this DirichletBC BC
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      assert(_offsetLocal[iPoint]+iDOF < numTotalConstrained);
      allFixedDOF[_offsetLocal[iPoint]+iDOF] = _fixedDOF[iDOF];
    } // for

    // Fill in rest of values not yet set (will be set by
    // another DirichletBC BC)
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
    section->setConstraintDof(point, &allFixedDOF[0]);
  } // for
} // setConstraints

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::DirichletBC::setField(const double t,
				  const topology::Field& field)
{ // setField
  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<MeshRealSection>& section = field.section();
  assert(!section.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = field.mesh().sieveMesh();
  assert(!sieveMesh.isNull());

  const int numPoints = _points.size();
  const int fiberDimension = 
    (numPoints > 0) ? section->getFiberDimension(_points[0]) : 0;
  double_array allValues(fiberDimension);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type point = _points[iPoint];
    assert(fiberDimension == section->getFiberDimension(point));
    sieveMesh->restrictClosure(section, point, &allValues[0], fiberDimension);
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      allValues[_fixedDOF[iDOF]] = _valuesInitial[iPoint*numFixedDOF+iDOF];
    if (t > _tRef && 0 != _dbRate)
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
	allValues[_fixedDOF[iDOF]] += 
	  (t-_tRef) * _valuesRate[iPoint*numFixedDOF+iDOF];
    section->updatePointAll(_points[iPoint], &allValues[0]);
  } // for
} // setField

// ----------------------------------------------------------------------
// Get mesh labels for points associated with Dirichlet BC.
void
pylith::bc::DirichletBC::_getPoints(const topology::Mesh& mesh)
{ // _getPoints
  typedef SieveMesh::int_section_type::chart_type chart_type;

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
// Setup initial and rate of change databases for querying.
void
pylith::bc::DirichletBC::_setupQueryDatabases(void)
{ // _setupQueryDatabases
  assert(0 != _db);

  const int numFixedDOF = _fixedDOF.size();
  char** valueNames = (numFixedDOF > 0) ? new char*[numFixedDOF] : 0;
  for (int i=0; i < numFixedDOF; ++i) {
    std::ostringstream name;
    name << "dof-" << _fixedDOF[i];
    const int size = 1 + name.str().length();
    valueNames[i] = new char[size];
    strcpy(valueNames[i], name.str().c_str());
  } // for

  // Setup initial database.
  _db->open();
  _db->queryVals(const_cast<const char**>(valueNames), numFixedDOF);

  // Setup rate database, if provided.
  if (0 != _dbRate) {
    _dbRate->open();
    _dbRate->queryVals((const char**) valueNames, numFixedDOF);
  } // if
  for (int i=0; i < numFixedDOF; ++i) {
    delete[] valueNames[i]; valueNames[i] = 0;
  } // for
  delete[] valueNames; valueNames = 0;
} // _setupQueryDatabases

// ----------------------------------------------------------------------
// Query initial and rate of change databases for values.
void
pylith::bc::DirichletBC::_queryDatabases(const topology::Mesh& mesh)
{ // _queryDatabases
  assert(0 != _db);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<MeshRealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double velocityScale = 
    _normalizer->lengthScale() / _normalizer->timeScale();

  const int numPoints = _points.size();
  const int numFixedDOF = _fixedDOF.size();
  _valuesInitial.resize(numPoints*numFixedDOF);
  if (0 != _dbRate)
    _valuesRate.resize(numPoints*numFixedDOF);

  double_array queryValues(numFixedDOF);
  double_array vCoordsGlobal(spaceDim);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    // Get coordinates of vertex
    coordinates->restrictPoint(_points[iPoint], 
			       &vCoordsGlobal[0], vCoordsGlobal.size());
    _normalizer->dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(),
				lengthScale);
    int err = _db->query(&queryValues[0], numFixedDOF, 
			 &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find initial values at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      _valuesInitial[numFixedDOF*iPoint+iDOF] = 
	_normalizer->nondimensionalize(queryValues[iDOF], lengthScale);

    if (0 != _dbRate) {
      err = _dbRate->query(&queryValues[0], numFixedDOF, 
			   &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find rate values at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << vCoordsGlobal[i];
	msg << ") using spatial database " << _dbRate->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
	_valuesRate[numFixedDOF*iPoint+iDOF] =
	  _normalizer->nondimensionalize(queryValues[iDOF], velocityScale);
    } // if
  } // for
  _db->close();
  if (0 != _dbRate)
    _dbRate->close();
} // _queryDatabases


// End of file 
