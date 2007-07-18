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

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

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
// Initialize boundary condition. Determine orienation and compute stress
// vector at integration points.
void
pylith::bc::Neumann::initialize(const ALE::Obj<ALE::Mesh>& mesh,
				const spatialdata::geocoords::CoordSys* cs,
				const double_array& upDir,
				const double_array& normalDir);
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _db);
  assert(!mesh.isNull());
  assert(0 != cs);

  if (3 != upDir.size())
    throw std::runtime_error("Up direction for surface orientation must be "
			     "a vector with 3 components.");
  if (3 != normalDir.size())
    throw std::runtime_error("Normal direction for surface orientation must be "
			     "a vector with 3 components.");

  // Extract submesh associated with surface
  const ALE::Obj<ALE::Mesh> submesh =
    ALE::Selection<ALE::Mesh>::submesh(mesh, mesh->getIntSection(_label));
  if (submesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << _label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  assert(!submesh.isNull());

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<sieve_type>& sieve = submesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = submesh->heightStratum(1);
  const int numCells = cells->size();

  // Create section for stress vector in global coordinates
  const int spaceDim = cs->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int stressVecSize = spaceDim * numQuadPts;
  // Not sure this part is correct
  _stressVecGlobal = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_stressVecGlobal.isNull());
  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const int numCorners = sieve->nCone(*c_iter, mesh->depth())->size();

  for(Mesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    _stressVecGlobal->setFiberDimension(*c_iter, stressVecSize);
  }
  mesh->allocate(_stressVecGlobal);

  // Set up orientation information
  const int orientationSize = spaceDim * spaceDim;
  double orientation[orientationSize];
  const int surfaceDim = mesh->getDimension()-1;
  orient_fn_type orientFn;
  const ALE::Obj<real_section_type>& coordinates =
    mesh->getRealSection("coordinates");
  switch (surfaceDim)
    { // switch
    case 0 :
      orientFn = _orient1D;
      break;
    case 1 :
      orientFn = _orient2D;
      break;
    case 2 :
      orientFn = _orient3D;
      break;
    default :
      assert(0);
    } // switch

  // Loop over cell faces, compute orientations, and then compute
  // corresponding stress vector in global coordinates.
  // Store values in _stressVecGlobal.
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    _quadrature->computeGeometry(mesh, coordinates, *c_iter);
    for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
    // FINISH FIXING FROM HERE

  const int_section_type::chart_type& chart = groupField->getChart();
  const int numPoints = chart.size();
  _points.resize(numPoints);
  int i = 0;
  for(int_section_type::chart_type::iterator c_iter = chart.begin();
      c_iter != chart.end();
      ++c_iter) {
    _points[i++] = *c_iter;
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
  for (int i=0; i < numFixedDOF; ++i) {
    delete[] valueNames[i]; valueNames[i] = 0;
  } // for
  delete[] valueNames; valueNames = 0;

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const int spaceDim = cs->spaceDim();

  _values.resize(numPoints*numFixedDOF);
  double_array queryValues(numFixedDOF);
  for (int iPoint=0, i=0; iPoint < numPoints; ++iPoint) {
    // Get coordinates of vertex
    const real_section_type::value_type* vCoords = 
      coordinates->restrictPoint(_points[iPoint]);
    int err = _db->query(&queryValues[0], numFixedDOF, vCoords, spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find values at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << _db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      _values[i++] = queryValues[iDOF];
  } // for
  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::Dirichlet::setConstraintSizes(const ALE::Obj<real_section_type>& field,
					  const ALE::Obj<ALE::Mesh>& mesh)
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
      msg << "Found overly constrained point while setting up constraints for Dirichlet "
	  << "boundary condition '" << _label << "'.\n" << "Number of DOF at point "
	  << _points[iPoint] << " is " << fiberDim << " and number of attempted constraints is "
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
pylith::bc::Dirichlet::setConstraints(const ALE::Obj<real_section_type>& field,
				      const ALE::Obj<ALE::Mesh>& mesh)
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

    // Add in the ones for this Dirichlet BC
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      allFixedDOF[_offsetLocal[iPoint]+iDOF] = _fixedDOF[iDOF];

    // Fill in rest of values not yet set (will be set by another Dirichlet BC)
    for (int iDOF=_offsetLocal[iPoint]+numFixedDOF; 
	 iDOF < numTotalConstrained; 
	 ++iDOF)
      allFixedDOF[_offsetLocal[iPoint]+iDOF] = 999;

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
pylith::bc::Dirichlet::setField(const double t,
				const ALE::Obj<real_section_type>& field,
				const ALE::Obj<ALE::Mesh>& mesh)
{ // setField
  assert(!field.isNull());
  assert(!mesh.isNull());

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const int numPoints = _points.size();
  for (int iPoint=0, i=0; iPoint < numPoints; ++iPoint) {
    const Mesh::point_type point = _points[iPoint];
    const int fiberDimension = field->getFiberDimension(point);
    double_array allValues(fiberDimension);
    mesh->restrict(field, point, &allValues[0], fiberDimension);
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      allValues[_fixedDOF[iDOF]] = _values[i++];
    field->updatePointAll(_points[iPoint], &allValues[0]);
  } // for
} // setField


// End of file 
