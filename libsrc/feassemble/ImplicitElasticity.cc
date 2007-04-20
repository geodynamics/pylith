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

#include "ImplicitElasticity.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/utils/array.hh" // USES double_array

#include "petscmat.h" // USES PetscMat
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ImplicitElasticity::ImplicitElasticity(void) :
  _material(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ImplicitElasticity::~ImplicitElasticity(void)
{ // destructor
  delete _material; _material = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::ImplicitElasticity::ImplicitElasticity(const ImplicitElasticity& i) :
  IntegratorImplicit(i),
  _material(0)
{ // copy constructor
  if (0 != i._material)
    _material = i._material->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set material.
void
pylith::feassemble::ImplicitElasticity::material(
				       const materials::ElasticMaterial* m)
{ // material
  delete _material; _material = (0 != m) ? m->clone() : 0;
} // material

// ----------------------------------------------------------------------
// Integrate portion of RHS that is assumed to remain constant over the
// entire time step.  At present, this only includes body forces.
// Note:  The vector produced by this computation remains unchanged for
// small strain computations, so it might be reasonable to keep the vector
// separate from the other contributions to the RHS.

void
pylith::feassemble::ImplicitElasticity::integrateConstant(
			      const ALE::Obj<real_section_type>& b,
			      const ALE::Obj<real_section_type>& dispT,
			      const ALE::Obj<Mesh>& mesh)
{ // integrateConstant
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(!b.isNull());
  assert(!dispT.isNull());
  assert(!mesh.isNull());

  /* Comment out the guts of this code until we figure out where grav
     section is coming from.  Code should be essentially correct, though.
     It may also be necessary to put in a switch to turn gravity off.

  // Get information about section
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  // Allocate vector for cell values (if necessary)
  _initCellVector();

  // Loop over cells.
  for (Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *cellIter);

    // Set cell data in material
    _material->initCellData(*cellIter, numQuadPts);

    // Reset element vector to zero
    _resetCellVector();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Get density at quadrature points for this cell
    const std::vector<double_array>& density = _material->calcDensity();

    // Compute action for cell

    // Compute action for element body forces
    if (!grav.isNull()) {
      const real_section_type::value_type* gravCell =
	mesh->restrict(grav, cell);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double valI = wt*basis[iQ+iBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim) {
	    _cellVector[iBlock+iDim] += valI*gravCell[iDim];
	  } // for
	} // for
      } // for
      PetscErrorCode err =
	PetscLogFlops(numQuadPts*(2+numBasis*(2+2*spaceDim)));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } // if
	      

    // Assemble cell contribution into field
    mesh->updateAdd(b, *cellIter, _cellVector);
  } // for
  */
} // integrateConstant

// ----------------------------------------------------------------------
// Compute residual for quasi-static finite elements.
void
pylith::feassemble::ImplicitElasticity::integrateResidual(
			      const ALE::Obj<real_section_type>& b,
			      const ALE::Obj<real_section_type>& dispT,
			      const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(!b.isNull());
  assert(!dispT.isNull());
  assert(!mesh.isNull());

  // Get information about section
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  // Allocate vector for cell values (if necessary)
  _initCellVector();

  // Allocate vector for total strain
  int tensorSize = 0;
  if (1 == cellDim)
    tensorSize = 1;
  else if (2 == cellDim)
    tensorSize = 3;
  else if (3 == cellDim)
    tensorSize = 6;
  else
    throw std::logic_error("Tensor size not implemented for given cellDim.");
  std::vector<double_array> totalStrain(numQuadPts);
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    totalStrain[iQuad].resize(tensorSize);
    totalStrain[iQuad] = 0.0;
  } // for

  // Loop over cells
  for (Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *cellIter);

    // Set cell data in material
    _material->initCellData(*cellIter, numQuadPts);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const real_section_type::value_type* dispTCell = 
      mesh->restrict(dispT, *cellIter);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    if (cellDim != spaceDim)
      throw std::logic_error("Not implemented yet.");

    // Compute B(transpose) * sigma, first computing strains
    if (1 == cellDim) {
      // Compute total strains
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis)
	  totalStrain[iQuad] += basisDeriv[iQ+iBasis] * dispTCell[iBasis];
      const std::vector<double_array>& stress = 
	_material->calcStress(totalStrain);

      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double s11 = stress[iQuad][0];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double N1 = wt*basisDeriv[iQ+iBasis*cellDim  ];
	  _cellVector[iBlock  ] -= N1*s11;
	} // for
      } // for
      PetscErrorCode err = PetscLogFlops(numQuadPts*(1+numBasis*5));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");

    } else if (2 == cellDim) {
      // Compute total strains
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  totalStrain[iQuad][0] += 
	    basisDeriv[iQ+iBasis  ] * dispTCell[iBasis  ];
	  totalStrain[iQuad][1] += 
	    basisDeriv[iQ+iBasis+1] * dispTCell[iBasis+1];
	  totalStrain[iQuad][2] += 
	    0.5 * (basisDeriv[iQ+iBasis+1] * dispTCell[iBasis  ] +
		   basisDeriv[iQ+iBasis  ] * dispTCell[iBasis+1]);
	} // for
      } // for
      const std::vector<double_array>& stress = 
	_material->calcStress(totalStrain);

      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double s11 = stress[iQuad][0];
	const double s22 = stress[iQuad][1];
	const double s12 = stress[iQuad][2];
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
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  totalStrain[iQuad][0] += 
	    basisDeriv[iQ+iBasis  ] * dispTCell[iBasis  ];
	  totalStrain[iQuad][1] += 
	    basisDeriv[iQ+iBasis+1] * dispTCell[iBasis+1];
	  totalStrain[iQuad][2] += 
	    basisDeriv[iQ+iBasis+2] * dispTCell[iBasis+2];
	  totalStrain[iQuad][3] += 
	    0.5 * (basisDeriv[iQ+iBasis+1] * dispTCell[iBasis  ] +
		   basisDeriv[iQ+iBasis  ] * dispTCell[iBasis+1]);
	  totalStrain[iQuad][4] += 
	    0.5 * (basisDeriv[iQ+iBasis+2] * dispTCell[iBasis+1] +
		   basisDeriv[iQ+iBasis+1] * dispTCell[iBasis+2]);
	  totalStrain[iQuad][5] += 
	    0.5 * (basisDeriv[iQ+iBasis+2] * dispTCell[iBasis  ] +
		   basisDeriv[iQ+iBasis  ] * dispTCell[iBasis+2]);
	} // for
      } // for
      const std::vector<double_array>& stress = 
	_material->calcStress(totalStrain);

      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double s11 = stress[iQuad][0];
	const double s22 = stress[iQuad][1];
	const double s33 = stress[iQuad][2];
	const double s12 = stress[iQuad][3];
	const double s23 = stress[iQuad][4];
	const double s13 = stress[iQuad][5];

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
} // integrateResidual


// ----------------------------------------------------------------------
// Compute stiffness matrix.
void
pylith::feassemble::ImplicitElasticity::integrateJacobian(
			     PetscMat* mat,
			     const ALE::Obj<real_section_type>& dispT,
			     const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != _mat);
  assert(!dispT.isNull());
  assert(!mesh.isNull());

  // Get information about section
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  const int numBasis = _quadrature->numCorners();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vector for cell values (if necessary)
  _initCellMatrix();

  // Loop over cells
  for (Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *cellIter);

    // Set cell data in material
    _material->initCellData(*cellIter, numQuadPts);

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute Jacobian for cell, specific for each geometry type
    if (cellDim != spaceDim)
      throw std::logic_error("Not implemented yet.")

	// Need to finish fixing from here***************************
    // 1D Case
    if (1 == cellDim) {
      assert(1 == numElasticConsts);
      // Compute strains
      /* for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * numBasis;
	  for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  const int jBlock = jBasis; */

      // Get "elasticity" matrix at quadrature points for this cell
      _material->calcDerivElastic(*cellIter, patch, numQuadPts);
      const double* elasticConsts = _material->elasticConsts();
      const int numElasticConsts = _material->numElasticConsts();

      // Compute Jacobian for consistent tangent matrix
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double C1111 = elasticConsts[iQuad];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * spaceDim;
	  const double valI = wt*basisDeriv[iQ+iBasis]*C1111;
	  for (int jBasis=0; jBasis < numBasis; ++jBasis] {
	    const int jBlock = jBasis * spaceDim;
	    const double valIJ = valI * basisDeriv[iQ+jBasis];
	    _cellMatrix[iBlock][jBlock] += valIJ;
	  } // for
	} // for
      } // for
      PetscErrorCode err = 
        PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*3)));
      if (err)
        throw std::runtime_error("Logging PETSc flops failed.");

    // 2D Case
    } else if (2 == cellDim) {
      assert(6 == numElasticConsts);
      // Compute strains
      /* for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * numBasis;
	  for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  const int jBlock = jBasis; */

      // Get "elasticity" matrix at quadrature points for this cell
      _material->calcDerivElastic(*cellIter, patch, numQuadPts);
      const double* elasticConsts = _material->elasticConsts();
      const int numElasticConsts = _material->numElasticConsts();

      // Compute Jacobian for consistent tangent matrix
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
            _cellMatrix[iBlock  ][jBlock  ] += ki0j0;
            _cellMatrix[iBlock  ][jBlock+1] += ki0j1;
            _cellMatrix[iBlock+1][jBlock  ] += ki0j1;
            _cellMatrix[iBlock+1][jBlock+1] += ki1j1;
          } // for
        } // for
      } // for
      PetscErrorCode err = 
        PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*(3*11+4))));
      if (err)
        throw std::runtime_error("Logging PETSc flops failed.");

    // 3D Case
    } else if (3 == cellDim) {
      assert(21 == numElasticConsts);
      // Compute strains
      /* for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const int iBlock = iBasis * numBasis;
	  for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  const int jBlock = jBasis; */

      // Get "elasticity" matrix at quadrature points for this cell
      _material->calcDerivElastic(*cellIter, patch, numQuadPts);
      const double* elasticConsts = _material->elasticConsts();
      const int numElasticConsts = _material->numElasticConsts();

      // Compute Jacobian for consistent tangent matrix
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

	    _cellMatrix[iblock  ][jBlock  ] += ki0j0;
	    _cellMatrix[iblock+1][jBlock  ] += ki0j1;
	    _cellMatrix[iblock+2][jBlock  ] += ki0j2;
	    _cellMatrix[iblock  ][jBlock+1] += ki0j1;
	    _cellMatrix[iblock+1][jBlock+1] += ki1j1;
	    _cellMatrix[iblock+2][jBlock+1] += ki1j2;
	    _cellMatrix[iblock  ][jBlock+2] += ki0j2;
	    _cellMatrix[iblock+1][jBlock+2] += ki1j2;
	    _cellMatrix[iblock+2][jBlock+2] += ki2j2;
          } // for
        } // for
      } // for
      PetscErrorCode err = 
        PetscLogFlops(numQuadPts*(1+numBasis*(3+numBasis*(6*26+9))));
      if (err)
        throw std::runtime_error("Logging PETSc flops failed.");
    } // if/else

    // Assemble cell contribution into field
    err = assembleMatrix(*mat, *cellIter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for
} // integrateJacobian


// End of file 
