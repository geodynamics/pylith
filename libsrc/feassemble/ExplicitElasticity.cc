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
// Integrate constant term (b) for dynamic elasticity term for 3-D
// finite elements.
void
pylith::feassemble::ExplicitElasticity::integrateConstant(
			      const ALE::Obj<real_section_type>& b,
			      const ALE::Obj<real_section_type>& dispT,
			      const ALE::Obj<real_section_type>& dispTmdt,
			      const ALE::Obj<Mesh>& mesh)
{ // integrateConstant
  assert(0 != _quadrature);
  assert(!b.isNull());
  assert(!dispT.isNull());
  assert(!dispTmdt.isNull());
  assert(!mesh.isNull());

  // Get information about section
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double* quadWts = _quadrature->quadWts();
  const int numBasis = _quadrature->numCorners();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  // Allocate vector for cell values (if necessary)
  _initCellVector();
  _material->initCellData(numQuadPts);

  for (Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *cellIter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const real_section_type::value_type* dispTCell = 
      mesh->restrict(dispT, *cellIter);
    const real_section_type::value_type* dispTmdtCell = 
      mesh->restrict(dispTmdt, *cellIter);

    // Get cell geometry information that depends on cell
    const double* basis = _quadrature->basis();
    const double* basisDeriv = _quadrature->basisDeriv();
    const double* jacobianDet = _quadrature->jacobianDet();

    // Compute action for cell

    // Compute action for inertial terms
    const double* density = _material->calcDensity(*cellIter, numQuadPts);
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
	      valIJ * (2.0 * dispTCell[jBlock+iDim] - 
		       dispTmdtCell[jBlock+iDim]);
        } // for
      } // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(5*spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    /** :TODO:
     *
     * If cellDim and spaceDim are different, we need to transform
     * displacements into cellDim, compute action, and transform
     * result back into spaceDim. Can we get this from the inverse of
     * the Jacobian?
     */
    if (cellDim != spaceDim)
      throw std::logic_error("Not implemented yet.");

    // Compute action for elastic terms
    if (1 == cellDim) {
      // Compute total strains
      const int stressSize = _material->stressSize();
      assert(1 == stressSize);
      const int strainSize = stressSize * numQuadPts;
      double* totalStrain = (strainSize > 0) ? new double[strainSize] : 0;
      memset(totalStrain, 0, strainSize*sizeof(double));
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis)
	  totalStrain[iQuad] += basisDeriv[iQ+iBasis] * dispTCell[iBasis];
      } // for
      const double* stress = 
	_material->calcStress(*cellIter, totalStrain, numQuadPts, 
			      spaceDim);

      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const int iStress = iQuad*stressSize;
	const double s11 = stress[iStress  ];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double N1 = wt*basisDeriv[iQ+iBasis*cellDim  ];
	  _cellVector[iBlock  ] -= N1*s11;
	} // for
      } // for
      delete[] totalStrain; totalStrain = 0;
      PetscErrorCode err = PetscLogFlops(numQuadPts*(1+numBasis*5));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } else if (2 == cellDim) {
      // Compute total strains
      const int stressSize = _material->stressSize();
      assert(3 == stressSize);
      const int strainSize = stressSize * numQuadPts;
      double* totalStrain = (strainSize > 0) ? new double[strainSize] : 0;
      memset(totalStrain, 0, strainSize*sizeof(double));
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  totalStrain[iQuad  ] += 
	    basisDeriv[iQ+iBasis  ] * dispTCell[iBasis  ];
	  totalStrain[iQuad+1] += 
	    basisDeriv[iQ+iBasis+1] * dispTCell[iBasis+1];
	  totalStrain[iQuad+2] += 
	    0.5 * (basisDeriv[iQ+iBasis+1] * dispTCell[iBasis  ] +
		   basisDeriv[iQ+iBasis  ] * dispTCell[iBasis+1]);
	} // for
      } // for
      const double* stress = 
	_material->calcStress(*cellIter, totalStrain, numQuadPts, 
			      spaceDim);
      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const int iStress = iQuad*stressSize;
	const double s11 = stress[iStress  ];
	const double s22 = stress[iStress+1];
	const double s12 = stress[iStress+2];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double N1 = wt*basisDeriv[iQ+iBasis*cellDim  ];
	  const double N2 = wt*basisDeriv[iQ+iBasis*cellDim+1];
	  _cellVector[iBlock  ] -= N1*s11 + N2*s12;
	  _cellVector[iBlock+1] -= N1*s12 + N2*s22;
	} // for
      } // for
      err = PetscLogFlops(numQuadPts*(1+numBasis*(8+2+9)));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } else if (3 == cellDim) {
      // Compute total strains
      const int stressSize = _material->stressSize();
      assert(6 == stressSize);
      const int strainSize = stressSize*numQuadPts;
      double* totalStrain = (strainSize > 0) ? new double[strainSize] : 0;
      memset(totalStrain, 0, strainSize*sizeof(double));
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  totalStrain[iQuad  ] += 
	    basisDeriv[iQ+iBasis  ] * dispTCell[iBasis  ];
	  totalStrain[iQuad+1] += 
	    basisDeriv[iQ+iBasis+1] * dispTCell[iBasis+1];
	  totalStrain[iQuad+2] += 
	    basisDeriv[iQ+iBasis+2] * dispTCell[iBasis+2];
	  totalStrain[iQuad+3] += 
	    0.5 * (basisDeriv[iQ+iBasis+1] * dispTCell[iBasis  ] +
		   basisDeriv[iQ+iBasis  ] * dispTCell[iBasis+1]);
	  totalStrain[iQuad+4] += 
	    0.5 * (basisDeriv[iQ+iBasis+2] * dispTCell[iBasis+1] +
		   basisDeriv[iQ+iBasis+1] * dispTCell[iBasis+2]);
	  totalStrain[iQuad+5] += 
	    0.5 * (basisDeriv[iQ+iBasis+2] * dispTCell[iBasis  ] +
		   basisDeriv[iQ+iBasis  ] * dispTCell[iBasis+2]);
	} // for
      } // for
      const double* stress = 
	_material->calcStress(*cellIter, totalStrain, numQuadPts, spaceDim);
      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const int iStress = iQuad*stressSize;
	const double s11 = stress[iStress  ];
	const double s22 = stress[iStress+1];
	const double s33 = stress[iStress+2];
	const double s12 = stress[iStress+3];
	const double s23 = stress[iStress+4];
	const double s13 = stress[iStress+5];

	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double N1 = wt*basisDeriv[iQ+iBasis*cellDim+0];
	  const double N2 = wt*basisDeriv[iQ+iBasis*cellDim+1];
	  const double N3 = wt*basisDeriv[iQ+iBasis*cellDim+2];
	  _cellVector[iBlock  ] -= N1*s11 + N2*s12 + N3*s13;
	  _cellVector[iBlock+1] -= N1*s12 + N2*s22 + N3*s23;
	  _cellVector[iBlock+2] -= N1*s13 + N2*s23 + N3*s33;
	} // for
      } // for
      err = PetscLogFlops(numQuadPts*(1+numBasis*(3+12)));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } // if/else

    // Assemble cell contribution into field
    mesh->updateAdd(b, *cellIter, _cellVector);
  } // for
} // integrateConstant

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ExplicitElasticity::integrateJacobian(
			     PetscMat* mat,
			     const ALE::Obj<real_section_type>& dispT,
			     const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != mat);
  assert(!dispT.isNull());
  assert(!mesh.isNull());

  // Get information about section
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double* quadWts = _quadrature->quadWts();
  const int numBasis = _quadrature->numCorners();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vector for cell values (if necessary)
  _initCellMatrix();
  _material->initCellData(numQuadPts);

  for (Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *cellIter);

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double* basis = _quadrature->basis();
    const double* jacobianDet = _quadrature->jacobianDet();

    // Get material physical properties at quadrature points for this cell
    const double* density = _material->calcDensity(*cellIter, numQuadPts);

    // Compute Jacobian for cell

    // Compute Jacobian for inertial terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad] / dt2;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim) {
            const int iBlock = (iBasis * spaceDim + iDim) * (spaceDim * numBasis);
            const int jBlock = (jBasis * spaceDim + iDim);
            _cellMatrix[iBlock+jBlock] += valIJ;
          } // for
        } // for
      } // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into field
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", dispT->getAtlas());

    err = updateOperator(*mat, mesh, dispT, globalOrder,
			 *cellIter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for
} // integrateJacobian


// End of file 
