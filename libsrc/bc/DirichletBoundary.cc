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

#include "DirichletBoundary.hh" // implementation of object methods

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <Selection.hh> // USES submesh algorithms

#include <string.h> // USES strcpy()
#include <strings.h> // USES strcasecmp()
#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletBoundary::DirichletBoundary(void) :
  _tRef(0.0),
  _dbRate(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBoundary::~DirichletBoundary(void)
{ // destructor
  _dbRate = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBoundary::initialize(
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

  // Extract submesh associated with boundary
  _boundaryMesh = 
    ALE::Selection<Mesh>::submeshV(mesh, mesh->getIntSection(_label));
  if (_boundaryMesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for Dirichlet boundary "
	<< "condition '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  _boundaryMesh->setRealSection("coordinates", 
				mesh->getRealSection("coordinates"));
  // Create the parallel overlap
  Obj<Mesh::send_overlap_type> sendParallelMeshOverlap = _boundaryMesh->getSendOverlap();
  Obj<Mesh::recv_overlap_type> recvParallelMeshOverlap = _boundaryMesh->getRecvOverlap();
  Mesh::renumbering_type&      renumbering             = mesh->getRenumbering();
  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<Mesh::point_type,Mesh::point_type> > globalPoints(renumbering);

  ALE::OverlapBuilder<>::constructOverlap(globalPoints, renumbering, sendParallelMeshOverlap, recvParallelMeshOverlap);
  _boundaryMesh->setCalculatedOverlap(true);

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

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _boundaryMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const int spaceDim = cs->spaceDim();

  _values = new real_section_type(_boundaryMesh->comm(), 
				  _boundaryMesh->debug());
  _values->addSpace(); // initial values
  _values->addSpace(); // rate of change of values
  _values->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(), vertices->end()), *std::max_element(vertices->begin(), vertices->end())+1));
  _values->setFiberDimension(vertices, 2*numFixedDOF);
  _values->setFiberDimension(vertices, numFixedDOF, 0); // initial values
  _values->setFiberDimension(vertices, numFixedDOF, 1); // rate of change
  _boundaryMesh->allocate(_values);

  double_array queryValues(2*numFixedDOF);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double velocityScale = 
    _normalizer->lengthScale() / _normalizer->timeScale();

  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    // Get coordinates of vertex
    const real_section_type::value_type* vCoords = 
      coordinates->restrictPoint(*v_iter);
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
    for (int i=0; i < numFixedDOF; ++i)
      _normalizer->nondimensionalize(queryValues[i], lengthScale);

    err = _dbRate->query(&queryValues[numFixedDOF], numFixedDOF, vCoords, 
			 spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find values at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << _dbRate->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    for (int i=0; i < numFixedDOF; ++i)
      _normalizer->nondimensionalize(queryValues[numFixedDOF+i], velocityScale);

    _values->updatePoint(*v_iter, &queryValues[0]);
  } // for
  _db->close();
  _dbRate->close();
} // initialize

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::DirichletBoundary::setConstraintSizes(
				     const ALE::Obj<real_section_type>& field,
				     const ALE::Obj<Mesh>& mesh)
{ // setConstraintSizes
  assert(!field.isNull());
  assert(!mesh.isNull());

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _boundaryMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  _offsetLocal = new int_section_type(_boundaryMesh->comm(), 
				      _boundaryMesh->debug());
  _offsetLocal->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(), vertices->end()), *std::max_element(vertices->begin(), vertices->end())+1));
  _offsetLocal->setFiberDimension(vertices, 1);
  _boundaryMesh->allocate(_offsetLocal);

  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    const int fiberDim = field->getFiberDimension(*v_iter);
    const int curNumConstraints = field->getConstraintDimension(*v_iter);
    if (curNumConstraints + numFixedDOF > fiberDim) {
      std::ostringstream msg;
      msg << "Found overly constrained point while setting up constraints "
	  << "for DirichletBoundary boundary condition '" << _label << "'.\n"
	  << "Number of DOF at point " << *v_iter << " is " << fiberDim
	  << " and number of attempted constraints is "
	  << curNumConstraints+numFixedDOF << ".";
      throw std::runtime_error(msg.str());
    } // if
    _offsetLocal->updatePoint(*v_iter, &curNumConstraints);
    field->addConstraintDimension(*v_iter, numFixedDOF);
  } // for
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::bc::DirichletBoundary::setConstraints(
				    const ALE::Obj<real_section_type>& field,
				    const ALE::Obj<Mesh>& mesh)
{ // setConstraints
  assert(!field.isNull());
  assert(!mesh.isNull());

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _boundaryMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    // Get list of currently constrained DOF
    const int* curFixedDOF = field->getConstraintDof(*v_iter);
    const int numTotalConstrained = field->getConstraintDimension(*v_iter);

    // Create array holding all constrained DOF
    int_array allFixedDOF(curFixedDOF, numTotalConstrained);

    const int_section_type::value_type* offset = 
      _offsetLocal->restrictPoint(*v_iter);

    // Verify other BC has not already constrained DOF
    const int numPrevious = offset[0];
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

    // Add in the ones for this DirichletBoundary BC
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      assert(offset[0]+iDOF < numTotalConstrained);
      allFixedDOF[offset[0]+iDOF] = _fixedDOF[iDOF];
    } // for

    // Fill in rest of values not yet set
    // (will be set by another DirichletBoundary BC)
    for (int iDOF=offset[0]+numFixedDOF; iDOF < numTotalConstrained; ++iDOF) {
      assert(iDOF < numTotalConstrained);
      allFixedDOF[iDOF] = 999;
    } // for

    // Sort list of constrained DOF
    // I need these sorted for my update algorithms to work properly
    std::sort(&allFixedDOF[0], &allFixedDOF[numTotalConstrained]);

    // Update list of constrained DOF
    field->setConstraintDof(*v_iter, &allFixedDOF[0]);
  } // for
} // setConstraints

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::DirichletBoundary::setField(const double t,
				      const ALE::Obj<real_section_type>& field,
				      const ALE::Obj<Mesh>& mesh)
{ // setField
  assert(!field.isNull());
  assert(!mesh.isNull());

  const int numFixedDOF = _fixedDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _boundaryMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  const int fiberDimension = 
    (vertices->size() > 0) ? field->getFiberDimension(*vertices->begin()) : 0;

  double_array fieldValues(fiberDimension);

  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    assert(fiberDimension == field->getFiberDimension(*v_iter));
    mesh->restrictClosure(field, *v_iter, &fieldValues[0], fiberDimension);

    const real_section_type::value_type* values = 
      _values->restrictPoint(*v_iter);

    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      fieldValues[_fixedDOF[iDOF]] = values[iDOF];
    if (t > _tRef)
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
	fieldValues[_fixedDOF[iDOF]] += (t-_tRef) * values[numFixedDOF+iDOF];
    field->updatePointAll(*v_iter, &fieldValues[0]);
  } // for
} // setField

// ----------------------------------------------------------------------
// Get vertex field of BC initial or rate of change of values.
const ALE::Obj<pylith::real_section_type>&
pylith::bc::DirichletBoundary::vertexField(VectorFieldEnum* fieldType,
					   const char* name,
					   const ALE::Obj<Mesh>& mesh,
					   topology::FieldsManager* const fields)
{ // getVertexField
  assert(0 != fieldType);
  assert(0 != name);
  assert(!_boundaryMesh.isNull());
  assert(!_values.isNull());  
  assert(0 != _normalizer);

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _boundaryMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  ALE::Obj<real_section_type> field = 0;
  int fiberDim = 0;
  double scale = 0.0;
  if (0 == strcasecmp(name, "initial")) {
    *fieldType = VECTOR_FIELD;
    field = _values->getFibration(0);
    fiberDim = 
      (vertices->size() > 0) ? field->getFiberDimension(*vertices->begin()) : 0;
    scale = _normalizer->lengthScale();
  } else if (0 == strcasecmp(name, "rate-of-change")) {
    *fieldType = VECTOR_FIELD;
    field = _values->getFibration(0);
    fiberDim = 
      (vertices->size() > 0) ? field->getFiberDimension(*vertices->begin()) : 0;
    scale = _normalizer->lengthScale() / _normalizer->timeScale();
  } else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for Dirichlet BC '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Allocate buffer if necessary
  if (_buffer.isNull()) {
    _buffer = new real_section_type(_boundaryMesh->comm(), 
				    _boundaryMesh->debug());
    _buffer->setChart(real_section_type::chart_type(
				    *std::min_element(vertices->begin(),
						      vertices->end()),
				    *std::max_element(vertices->begin(),
						      vertices->end())+1));
    _buffer->setFiberDimension(vertices, fiberDim);
    _boundaryMesh->allocate(_buffer);
  } // if

  // dimensionalize values
  double_array pointValues(fiberDim);
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    assert(fiberDim == field->getFiberDimension(*v_iter));
    assert(fiberDim == _buffer->getFiberDimension(*v_iter));
    const real_section_type::value_type* values = 
      field->restrictPoint(*v_iter);
    for (int i=0; i < fiberDim; ++i)
      pointValues[i] = _normalizer->dimensionalize(values[i], scale);
    _buffer->updatePointAll(*v_iter, &pointValues[0]);
  } // for

  return _buffer;
} // getVertexField


// End of file 
