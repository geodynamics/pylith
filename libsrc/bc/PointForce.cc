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
#include "pylith/topology/Fields.hh" // USES Fields
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
// Initialize boundary condition.
void
pylith::bc::PointForce::initialize(const topology::Mesh& mesh,
				    const double upDir[3])
{ // initialize
  if (0 == _bcDOF.size())
    return;

  _getPoints(mesh);

  assert(0 != _normalizer);
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
// Verify configuration is acceptable.
void
pylith::bc::PointForce::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BoundaryCondition::verifyConfiguration(mesh);
  TimeDependent::verifyConfiguration(mesh);
} // verifyConfiguration


// End of file 
