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

#include "PointForce.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
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
pylith::bc::PointForce::PointForce(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::PointForce::~PointForce(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::PointForce::deallocate(void)
{ // deallocate
  TimeDependentPoints::deallocate();
  feassemble::Integrator<feassemble::Quadrature<topology::Mesh> >::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::PointForce::initialize(const topology::Mesh& mesh,
				    const PylithScalar upDir[3])
{ // initialize
  if (0 == _bcDOF.size())
    return;

  _getPoints(mesh);

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;

  _queryDatabases(mesh, forceScale, "force");
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::PointForce::integrateResidual(
			   const topology::Field<topology::Mesh>& residual,
			   const PylithScalar t,
			   topology::SolutionFields* const fields)
{ // integrateResidualAssembled
  assert(_parameters);
  assert(_normalizer);

  // Calculate spatial and temporal variation of value for BC.
  _calculateValue(t);

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();

  const topology::Mesh& mesh = residual.mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  scalar_array residualVertex(spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  // Get global order
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
        residualSection);
  assert(!globalOrder.isNull());

  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(valueIndex+valueFiberDim <= parametersFiberDim);
  assert(valueFiberDim == numBCDOF);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label.

    // Contribute to residual only if point is local.
    if (!globalOrder->isLocal(p_bc))
      continue;

    residualVertex *= 0.0; // Reset residual contribution to zero.
    
    assert(parametersFiberDim == parametersSection->getFiberDimension(p_bc));
    const PylithScalar* parametersVertex = parametersSection->restrictPoint(p_bc);
    assert(parametersVertex);

    for (int iDOF=0; iDOF < numBCDOF; ++iDOF)
      residualVertex[_bcDOF[iDOF]] += parametersVertex[valueIndex+iDOF];

    assert(residualVertex.size() == residualSection->getFiberDimension(p_bc));
    residualSection->updateAddPoint(p_bc, &residualVertex[0]);
  } // for
} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::PointForce::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BoundaryCondition::verifyConfiguration(mesh);
  TimeDependent::verifyConfiguration(mesh);
} // verifyConfiguration


// End of file 
