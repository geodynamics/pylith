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
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
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
pylith::bc::PointForce::PointForce(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::PointForce::~PointForce(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set indices of vertices with point forces.
void
pylith::bc::PointForce::forceDOF(const int* flags,
				 const int size)
{ // forceDOF
  if (size > 0)
    assert(0 != flags);

  _forceDOF.resize(size);
  for (int i=0; i < size; ++i)
    _forceDOF[i] = flags[i];
} // forceDOF

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::PointForce::initialize(const topology::Mesh& mesh,
				    const double upDir[3])
{ // initialize
  const int numForceDOF = _forceDOF.size();
  if (0 == numForceDOF)
    return;

  _getPoints(mesh);
  _setupQueryDatabases();
  _queryDatabases(mesh);
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator that
// do not require assembly over cells, vertices, or processors.
void
integrateResidualAssembled(const topology::Field<topology::Mesh>& residual,
			   const double t,
			   topology::SolutionFields* const fields)
{ // integrateResidualAssembled
  const int numPoints = _points.size();
  const int numForceDOF = _forceDOF.size();

  const ALE::Obj<RealSection>& amplitudeSection = 
    _parameters->get("amplitude").section();
  assert(!amplitudeSection.isNull());

  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  const ALE::Obj<RealSection>& startTimeSection =  (0 == _dbStartTime) ?
    0 : _parameters->get("start time").section();

  double_array forcesVertex(spaceDim);
  double scale = 1.0;
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_force = _points[iPoint];
    assert(numForceDOF == amplitudeSection->getFiberDimension(p_force));
    const double* amplitudeVertex = amplitudeSection->restrictPoint(p_force);
    assert(0 != amplitudeVertex);

    if (0 != _dbStartTime) {
      assert(!startTimeSection.isNull()); // Expect section with start times
      assert(0 != _dbTimeAmp); // Expect database with temporal evolution

      assert(1 == startTimeSection->getFiberDimension(p_force));
      const double tRef = *startTimeSection->restrictPoint(p_force);
      const double tRel = t - tRef;
      if (tRel > 0.0) {
	_normalizer->dimensionalize(&tRel, 1, timeScale);
	scale = _dbTimeAmp->query(tRel);
      } // if
    } // if
    
    for (int iDOF=0; iDOF < numForceDOF; ++iDOF)
      forcesVertex[_forceDOF[iDOF]] = amplitudeVertex[iDOF]*scale;
  } // for
} // interateResidualAssembled

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
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
// Setup initial and rate of change databases for querying.
void
pylith::bc::PointForce::_setupQueryDatabases(void)
{ // _setupQueryDatabases
  assert(0 != _db);

  const int numForceDOF = _forceDOF.size();
  char** valueNames = (numForceDOF > 0) ? new char*[numForceDOF] : 0;
  for (int i=0; i < numForceDOF; ++i) {
    std::ostringstream name;
    name << "dof-" << _forceDOF[i];
    const int size = 1 + name.str().length();
    valueNames[i] = new char[size];
    strcpy(valueNames[i], name.str().c_str());
  } // for

  // Setup initial database.
  _db->open();
  _db->queryVals(const_cast<const char**>(valueNames), numForceDOF);

  // Setup rate database, if provided.
  if (0 != _dbRate) {
    _dbRate->open();
    _dbRate->queryVals((const char**) valueNames, numForceDOF);
  } // if
  for (int i=0; i < numForceDOF; ++i) {
    delete[] valueNames[i]; valueNames[i] = 0;
  } // for
  delete[] valueNames; valueNames = 0;
} // _setupQueryDatabases

// ----------------------------------------------------------------------
// Query initial and rate of change databases for values.
void
pylith::bc::PointForce::_queryDatabases(const topology::Mesh& mesh)
{ // _queryDatabases
  assert(0 != _db);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<RealSection>& coordinates = 
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
  const int numForceDOF = _forceDOF.size();
  _valuesInitial.resize(numPoints*numForceDOF);
  if (0 != _dbRate)
    _valuesRate.resize(numPoints*numForceDOF);

  double_array queryValues(numForceDOF);
  double_array vCoordsGlobal(spaceDim);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    // Get coordinates of vertex
    coordinates->restrictPoint(_points[iPoint], 
			       &vCoordsGlobal[0], vCoordsGlobal.size());
    _normalizer->dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(),
				lengthScale);
    int err = _db->query(&queryValues[0], numForceDOF, 
			 &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find initial values at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    for (int iDOF=0; iDOF < numForceDOF; ++iDOF)
      _valuesInitial[numForceDOF*iPoint+iDOF] = 
	_normalizer->nondimensionalize(queryValues[iDOF], lengthScale);
  } // for
  _db->close();
} // _queryDatabases


// End of file 
