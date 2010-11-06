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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletBC.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletBC::DirichletBC(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBC::~DirichletBC(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletBC::deallocate(void)
{ // deallocate
  TimeDependentPoints::deallocate();
  feassemble::Constraint::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBC::initialize(const topology::Mesh& mesh,
				    const double upDir[3])
{ // initialize
  if (0 == _bcDOF.size())
    return;

  _getPoints(mesh);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  _queryDatabases(mesh, lengthScale, "displacement");
} // initialize

// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraintSizes(const topology::Field<topology::Mesh>& field)
{ // setConstraintSizes
  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<RealSection>& section = field.section();
  assert(!section.isNull());

  // Set constraints in field
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

  // We only worry about the conventional DOF in a split field associated with
  // fibrations.
  const int numFibrations = section->getNumSpaces(); // > 1 for split field
  if (numFibrations > 0)
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      const int fibration = _bcDOF[iDOF];
      for (int iPoint=0; iPoint < numPoints; ++iPoint)
        section->addConstraintDimension(_points[iPoint], 1, fibration);
    } // for
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::bc::DirichletBC::setConstraints(const topology::Field<topology::Mesh>& field)
{ // setConstraints
  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  const ALE::Obj<RealSection>& section = field.section();
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
	if (allFixedDOF[iDOF] == _bcDOF[jDOF]) {
	  std::ostringstream msg;
	  msg << "Found multiple constraints on degrees of freedom at\n"
	      << "point while setting up constraints for DirichletBC\n"
	      << "boundary condition '" << _label << "'.\n"
	      << "Degree of freedom " << _bcDOF[jDOF] 
	      << " is already constrained by another Dirichlet BC.";
	  throw std::runtime_error(msg.str());
	} // if

    // Add in the ones for this DirichletBC BC
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      assert(_offsetLocal[iPoint]+iDOF < numTotalConstrained);
      allFixedDOF[_offsetLocal[iPoint]+iDOF] = _bcDOF[iDOF];
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

  // We only worry about the conventional DOF in a split field associated with
  // fibrations.
  const int numFibrations = section->getNumSpaces(); // > 1 for split field
  if (numFibrations > 0) {
    int zero = 0;
    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
      const int fibration = _bcDOF[iDOF];
      for (int iPoint=0; iPoint < numPoints; ++iPoint)
        section->setConstraintDof(_points[iPoint], &zero, fibration);
    } // for
  } // if
} // setConstraints

// ----------------------------------------------------------------------
// Set values in field.
void
pylith::bc::DirichletBC::setField(const double t,
				  const topology::Field<topology::Mesh>& field)
{ // setField
  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  // Calculate spatial and temporal variation of value for BC.
  _calculateValue(t);

  const int numPoints = _points.size();

  assert(_parameters);
  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(valueFiberDim == numFixedDOF);
  assert(valueIndex+valueFiberDim <= parametersFiberDim);

  const ALE::Obj<RealSection>& fieldSection = field.section();
  assert(!fieldSection.isNull());
  const int fiberDimension = 
    (numPoints > 0) ? fieldSection->getFiberDimension(_points[0]) : 0;
  double_array fieldVertex(fiberDimension);
  
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type p_bc = _points[iPoint];

    fieldSection->restrictPoint(p_bc, &fieldVertex[0], fieldVertex.size());

    assert(parametersFiberDim == parametersSection->getFiberDimension(p_bc));
    const double* parametersVertex = parametersSection->restrictPoint(p_bc);

    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      fieldVertex[_bcDOF[iDOF]] = parametersVertex[valueIndex+iDOF];

    fieldSection->updatePointAll(p_bc, &fieldVertex[0]);
  } // for
} // setField

// ----------------------------------------------------------------------
// Set increment in values from t0 to t1 in field.
void
pylith::bc::DirichletBC::setFieldIncr(const double t0,
				      const double t1,
				      const topology::Field<topology::Mesh>& field)
{ // setFieldIncr
  assert(_useSolnIncr);

  const int numFixedDOF = _bcDOF.size();
  if (0 == numFixedDOF)
    return;

  // Calculate spatial and temporal variation of value for BC.
  _calculateValueIncr(t0, t1);

  const int numPoints = _points.size();

  assert(_parameters);
  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(valueFiberDim == numFixedDOF);
  assert(valueIndex+valueFiberDim <= parametersFiberDim);

  const ALE::Obj<RealSection>& fieldSection = field.section();
  assert(!fieldSection.isNull());
  const int fiberDimension = 
    (numPoints > 0) ? fieldSection->getFiberDimension(_points[0]) : 0;
  double_array fieldVertex(fiberDimension);
  
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const SieveMesh::point_type p_bc = _points[iPoint];

    fieldSection->restrictPoint(p_bc, &fieldVertex[0], fieldVertex.size());

    assert(parametersFiberDim == parametersSection->getFiberDimension(p_bc));
    const double* parametersVertex = parametersSection->restrictPoint(p_bc);

    for (int iDOF=0; iDOF < numFixedDOF; ++iDOF)
      fieldVertex[_bcDOF[iDOF]] = parametersVertex[valueIndex+iDOF];

    fieldSection->updatePointAll(p_bc, &fieldVertex[0]);
  } // for

} // setFieldIncr

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletBC::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BoundaryCondition::verifyConfiguration(mesh);
  TimeDependent::verifyConfiguration(mesh);
} // verifyConfiguration


// End of file 
