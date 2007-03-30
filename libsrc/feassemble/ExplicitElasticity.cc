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
                  const ALE::Obj<Mesh>& m,
			      const ALE::Obj<real_section_type>& b,
			      const ALE::Obj<real_section_type>& dispT,
			      const ALE::Obj<real_section_type>& dispTmdt,
			      const ALE::Obj<real_section_type>& coordinates)
{ // integrateConstant
  assert(0 != _quadrature);

  // Get information about section
  const ALE::Obj<Mesh::label_sequence>& cells    = m->heightStratum(0);
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;

  // Allocate vector for cell values (if necessary)
  _initCellVector();

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double* quadWts = _quadrature->quadWts();
  const int numBasis = _quadrature->numCorners();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  for (Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(m, coordinates, *cellIter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const real_section_type::value_type* dispTCell = 
      m->restrict(dispT, *cellIter);
    const real_section_type::value_type* dispTmdtCell = 
      m->restrict(dispTmdt, *cellIter);

    // Get cell geometry information that depends on cell
    const double* basis = _quadrature->basis();
    const double* basisDeriv = _quadrature->basisDeriv();
    const double* jacobianDet = _quadrature->jacobianDet();

    // Get material physical properties at quadrature points for this cell
    _material->calcProperties(*cellIter, numQuadPts);
    const double* density = _material->density();
    const double* elasticConsts = _material->elasticConsts();
    const int numElasticConsts = _material->numElasticConsts();

    // Compute action for cell

    // Compute action for inertial terms
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
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+4*spaceDim))));
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
      assert(1 == numElasticConsts);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double C1111 = elasticConsts[iQuad];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double valI = wt*basisDeriv[iQ+iBasis]*C1111;
	  for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	    const int jBlock = jBasis * spaceDim;
	    const double valIJ = valI * basisDeriv[iQ+jBasis];
	    _cellVector[iBlock] -= valIJ * dispTCell[jBlock];
	  } // for
	} // for
      } // for      
      PetscErrorCode err = 
	PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*3)));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } else if (2 == cellDim) {
      assert(6 == numElasticConsts);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const int iConst = iQuad*numElasticConsts;
	const double C1111 = elasticConsts[iConst+0];
	const double C1122 = elasticConsts[iConst+1];
	const double C1112 = elasticConsts[iConst+2];
	const double C2222 = elasticConsts[iConst+3];
	const double C2212 = elasticConsts[iConst+4];
	const double C1212 = elasticConsts[iConst+5];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double Nip = wt*basisDeriv[iQ+iBasis*cellDim+0];
	  const double Niq = wt*basisDeriv[iQ+iBasis*cellDim+1];
	  for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	    const int jBlock = jBasis * spaceDim;
	    const double Njp = basisDeriv[iQ+jBasis*cellDim+0];
	    const double Njq = basisDeriv[iQ+jBasis*cellDim+1];
	    const double ki0j0 = 
	      C1111 * Nip * Njp + C1112 * Niq * Njp +
	      C1112 * Nip * Njq + C1212 * Niq * Njq;
	    const double ki0j1 =
	      C1122 * Nip * Njq + C2212 * Niq * Njq +
	      C1112 * Nip * Njp + C1212 * Niq * Njp;
	    const double ki1j1 =
	      C2222 * Niq * Njq + C2212 * Nip * Njq +
	      C2212 * Niq * Njp + C1212 * Nip * Njp;
	    _cellVector[iBlock  ] -=
	      ki0j0 * dispTCell[jBlock  ] + ki0j1 * dispTCell[jBlock+1];
	    _cellVector[iBlock+1] -=
	      ki0j1 * dispTCell[jBlock  ] + ki1j1 * dispTCell[jBlock+1];
	  } // for
	} // for
      } // for
      PetscErrorCode err = 
	PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*(2+3*11+2*4))));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } else if (3 == cellDim) {
      assert(21 == numElasticConsts);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const int iConst = iQuad*numElasticConsts;
	const double C1111 = elasticConsts[iConst+ 0];
	const double C1122 = elasticConsts[iConst+ 1];
	const double C1133 = elasticConsts[iConst+ 2];
	const double C1112 = elasticConsts[iConst+ 3];
	const double C1123 = elasticConsts[iConst+ 4];
	const double C1113 = elasticConsts[iConst+ 5];
	const double C2222 = elasticConsts[iConst+ 6];
	const double C2233 = elasticConsts[iConst+ 7];
	const double C2212 = elasticConsts[iConst+ 8];
	const double C2223 = elasticConsts[iConst+ 9];
	const double C2213 = elasticConsts[iConst+10];
	const double C3333 = elasticConsts[iConst+11];
	const double C3312 = elasticConsts[iConst+12];
	const double C3323 = elasticConsts[iConst+13];
	const double C3313 = elasticConsts[iConst+14];
	const double C1212 = elasticConsts[iConst+15];
	const double C1223 = elasticConsts[iConst+16];
	const double C1213 = elasticConsts[iConst+17];
	const double C2323 = elasticConsts[iConst+18];
	const double C2313 = elasticConsts[iConst+19];
	const double C1313 = elasticConsts[iConst+20];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double Nip = wt*basisDeriv[iQ+iBasis*cellDim+0];
	  const double Niq = wt*basisDeriv[iQ+iBasis*cellDim+1];
	  const double Nir = wt*basisDeriv[iQ+iBasis*cellDim+2];
	  for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	    const int jBlock = jBasis * spaceDim;
	    const double Njp = basisDeriv[iQ+jBasis*cellDim+0];
	    const double Njq = basisDeriv[iQ+jBasis*cellDim+1];
	    const double Njr = basisDeriv[iQ+jBasis*cellDim+2];
	    const double ki0j0 = 
	      C1111 * Nip * Njp + C1112 * Niq * Njp + C1113 * Nir * Njp +
	      C1112 * Nip * Njq + C1212 * Niq * Njq + C1213 * Nir * Njq +
	      C1113 * Nip * Njr + C1213 * Niq * Njr + C1313 * Nir * Njr;
	    const double ki0j1 =
	      C1122 * Nip * Njq + C2212 * Niq * Njq + C2213 * Nir * Njq +
	      C1112 * Nip * Njp + C1212 * Niq * Njp + C1213 * Nir * Njp +
	      C1123 * Nip * Njr + C1223 * Niq * Njr + C2313 * Nir * Njr;
	    const double ki0j2 =
	      C1133 * Nip * Njr + C3312 * Niq * Njr + C3313 * Nir * Njr +
	      C1123 * Nip * Njq + C1223 * Niq * Njq + C2313 * Nir * Njq +
	      C1113 * Nip * Njp + C1213 * Niq * Njp + C1313 * Nir * Njp;
	    const double ki1j1 =
	      C2222 * Niq * Njq + C2212 * Nip * Njq + C2223 * Nir * Njq +
	      C2212 * Niq * Njp + C1212 * Nip * Njp + C1223 * Nir * Njp +
	      C2223 * Niq * Njr + C1223 * Nip * Njr + C2323 * Nir * Njr;
	    const double ki1j2 =
	      C2233 * Niq * Njr + C3312 * Nip * Njr + C3323 * Nir * Njr +
	      C2223 * Niq * Njq + C1223 * Nip * Njq + C2323 * Nir * Njq +
	      C2213 * Niq * Njp + C1213 * Nip * Njp + C2313 * Nir * Njp;
	    const double ki2j2 =
	      C3333 * Nir * Njr + C3323 * Niq * Njr + C3313 * Nip * Njr +
	      C3323 * Nir * Njq + C2323 * Niq * Njq + C2313 * Nip * Njq +
	      C3313 * Nir * Njp + C2313 * Niq * Njp + C1313 * Nip * Njp;

	    _cellVector[iBlock  ] -=
	      ki0j0 * dispTCell[jBlock  ] + 
	      ki0j1 * dispTCell[jBlock+1] +
	      ki0j2 * dispTCell[jBlock+2];
	    _cellVector[iBlock+1] -=
	      ki0j1 * dispTCell[jBlock  ] + 
	      ki1j1 * dispTCell[jBlock+1] +
	      ki1j2 * dispTCell[jBlock+2];
	    _cellVector[iBlock+2] -=
	      ki0j2 * dispTCell[jBlock  ] + 
	      ki1j2 * dispTCell[jBlock+1] +
	      ki2j2 * dispTCell[jBlock+2];
	  } // for
	} // for
      } // for
      PetscErrorCode err = 
	PetscLogFlops(numQuadPts*(1+numBasis*(3+numBasis*(3+6*26+3*6))));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } // if/else

    // Assemble cell contribution into field
    m->updateAdd(b, *cellIter, _cellVector);
  } // for
} // integrateConstant

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ExplicitElasticity::integrateJacobian(
			     PetscMat* mat,
                 const ALE::Obj<Mesh>& m,
			     const ALE::Obj<real_section_type>& dispT,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrateJacobian
  assert(0 != _quadrature);

  // Get information about section
  const ALE::Obj<Mesh::label_sequence>& cells    = m->heightStratum(0);
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;

  // Allocate vector for cell values (if necessary)
  _initCellMatrix();

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double* quadWts = _quadrature->quadWts();
  const int numBasis = _quadrature->numCorners();
  const int spaceDim = _quadrature->spaceDim();
  for (Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(m, coordinates, *cellIter);

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double* basis = _quadrature->basis();
    const double* jacobianDet = _quadrature->jacobianDet();

    // Get material physical properties at quadrature points for this cell
    _material->calcProperties(*cellIter, numQuadPts);
    const double* density = _material->density();

    // Compute Jacobian for cell

    // Compute Jacobian for inertial terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad] / dt2;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const int iBlock = (iBasis * spaceDim) * (spaceDim * numBasis);
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const int jBlock = (jBasis * spaceDim);
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellMatrix[iBlock+jBlock+iDim*spaceDim+iDim] += valIJ;
        } // for
      } // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into field
    const ALE::Obj<Mesh::order_type>& globalOrder = m->getFactory()->getGlobalOrder(m, "default", dispT->getAtlas());

    err = updateOperator(*mat, m, dispT, globalOrder, *cellIter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for
} // integrateJacobian


// End of file 
