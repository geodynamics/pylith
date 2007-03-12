// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "ExplicitElasticity.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial

#include "petscmat.h" // USES PetscMat
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ExplicitElasticity::ExplicitElasticity(void) :
  _material(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ExplicitElasticity::~ExplicitElasticity(void)
{ // destructor
  delete _material; _material = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::ExplicitElasticity::ExplicitElasticity(const ExplicitElasticity& i) :
  IntegratorExplicit(i),
  _material(0)
{ // copy constructor
  if (0 != i._material)
    _material = i._material->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set material.
void
pylith::feassemble::ExplicitElasticity::material(
				       const materials::ElasticMaterial* m)
{ // material
  delete _material; _material = (0 != m) ? m->clone() : 0;
} // material

// ----------------------------------------------------------------------
// Integrate residual term (b) for dynamic elasticity term for 3-D
// finite elements.
void
pylith::feassemble::ExplicitElasticity::integrateResidual(
			      const ALE::Obj<real_section_type>& residual,
			      const ALE::Obj<real_section_type>& dispT,
			      const ALE::Obj<real_section_type>& dispTmdt,
			      const ALE::Obj<real_section_type>& coordinates)
{ // integrateResidual
  assert(0 != _quadrature);

  // Get information about section
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = dispT->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

  // Get parameters used in integration.
  const double dt = _dt;

  // Allocate vector for cell values (if necessary)
  _initCellVector();

  for (topology_type::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(coordinates, *cellIter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const real_section_type::value_type* dispTCell = 
      dispT->restrict(patch, *cellIter);
    const real_section_type::value_type* dispTmdtCell = 
      dispTmdt->restrict(patch, *cellIter);

    // Get cell geometry information
    const int numQuadPts = _quadrature->numQuadPts();
    const double* basis = _quadrature->basis();
    const double* quadPts = _quadrature->quadPts();
    const double* quadWts = _quadrature->quadWts();
    const double* jacobianDet = _quadrature->jacobianDet();
    const int numBasis = _quadrature->numCorners();
    const int spaceDim = _quadrature->spaceDim();

    // Get material physical properties at quadrature points for this cell
    _material->calcProperties(*cellIter, patch, numQuadPts);
    const double* density = _material->density();
    const double* elasticConsts = _material->elasticConsts();

    // Compute action for cell

    // Compute action for inertial terms
    const double dt2 = dt*dt;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad] / dt2;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const int iBlock = iBasis * spaceDim;
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const int jBlock = jBasis * spaceDim;
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBlock+iDim] += 
	      valIJ * 2.0 * (dispTCell[jBlock+iDim] - 
			     dispTmdtCell[jBlock+iDim]);
        } // for
      } // for
    } // for

    // Compute action for elastic terms
    // ADD STUFF HERE

    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+4*spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into field
    residual->updateAdd(patch, *cellIter, _cellVector);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ExplicitElasticity::integrateJacobian(
			     PetscMat* mat,
			     const ALE::Obj<real_section_type>& dispT,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrateJacobian
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute lumped matrix associated with operator.
void
pylith::feassemble::ExplicitElasticity::integrateJacobian(
			     const ALE::Obj<real_section_type>& fieldOut,
			     const ALE::Obj<real_section_type>& dispT,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrateJacobian
} // integrateJacobian

// ----------------------------------------------------------------------
// Setup material property parameters by querying database.
void
pylith::feassemble::ExplicitElasticity::initialize(ALE::Obj<ALE::Mesh>& mesh,
						   spatialdata::geocoords::CoordSys* cs)
{ // initialize
  assert(0 != _material);
  _material->initialize(mesh, cs, _quadrature);
} // initialize


// End of file 
