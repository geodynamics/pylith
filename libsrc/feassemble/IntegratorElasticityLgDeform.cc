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

#include "IntegratorElasticityLgDeform.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/array.hh" // USES double_array

#include <cstring> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorElasticityLgDeform::IntegratorElasticityLgDeform(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorElasticityLgDeform::~IntegratorElasticityLgDeform(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Determine whether we need to recompute the Jacobian.
bool
pylith::feassemble::IntegratorElasticityLgDeform::needNewJacobian(void)
{ // needNewJacobian
  _needNewJacobian = IntegratorElasticity::needNewJacobian();
  return _needNewJacobian;
} // needNewJacobian

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorElasticityLgDeform::updateStateVars(
				      const double t,
				      topology::SolutionFields* const fields)
{ // updateStateVars
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != fields);

  // No need to update state vars if material doesn't have any.
  if (!_material->hasStateVars())
    return;

  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D;
  } else {
      std::cerr << "Bad cell dimension '" << cellDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad cell dimension in IntegratorElasticityLgDeform.");
  } // else

  // Allocate arrays for cell data.
  double_array dispCell(numBasis*spaceDim);
  double_array strainCell(numQuadPts*tensorSize);
  double_array deformCell(numQuadPts*spaceDim*spaceDim);
  deformCell = 0.0;
  strainCell = 0.0;

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get fields
  const topology::Field<topology::Mesh>& disp = fields->get("disp(t)");
  const ALE::Obj<RealSection>& dispSection = disp.section();
  assert(!dispSection.isNull());
  topology::Mesh::RestrictVisitor dispVisitor(*dispSection, 
					      dispCell.size(), &dispCell[0]);

  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Retrieve geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Restrict input fields to cell
    dispVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispVisitor);

    // Get cell geometry information that depends on cell
    const double_array& basisDeriv = _quadrature->basisDeriv();
  
    // Compute deformation tensor.
    _calcDeformation(&deformCell, basisDeriv, coordinatesCell, dispCell,
		     numBasis, numQuadPts, spaceDim);

    // Compute strains
    calcTotalStrainFn(&strainCell, deformCell, numQuadPts);

    // Update material state
    _material->updateStateVars(strainCell, *c_iter);
  } // for
} // updateStateVars

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticityLgDeform::_calcStrainStressField(
				 topology::Field<topology::Mesh>* field,
				 const char* name,
				 topology::SolutionFields* const fields)
{ // _calcStrainStressField
  assert(0 != field);
  assert(0 != _quadrature);
  assert(0 != _material);

  const bool calcStress = (0 == strcasecmp(name, "stress")) ? true : false;
    
  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D;
  } else {
      std::cerr << "Bad cell dimension '" << cellDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad cell dimension in IntegratorElasticityLgDeform.");
  } // else
  
  // Allocate arrays for cell data.
  double_array dispCell(numBasis*spaceDim);
  double_array deformCell(numQuadPts*spaceDim*spaceDim);
  double_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  double_array stressCell(numQuadPts*tensorSize);
  stressCell = 0.0;

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = field->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get field
  const topology::Field<topology::Mesh>& disp = fields->get("disp(t)");
  const ALE::Obj<RealSection>& dispSection = disp.section();
  assert(!dispSection.isNull());
  topology::Mesh::RestrictVisitor dispVisitor(*dispSection, 
					      dispCell.size(), &dispCell[0]);
    
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  const ALE::Obj<RealSection>& fieldSection = field->section();
  assert(!fieldSection.isNull());

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Retrieve geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Restrict input fields to cell
    dispVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispVisitor);

    // Get cell geometry information that depends on cell
    const double_array& basisDeriv = _quadrature->basisDeriv();
    
    // Compute deformation tensor.
    _calcDeformation(&deformCell, basisDeriv, coordinatesCell, dispCell,
		     numBasis, numQuadPts, spaceDim);

    // Compute strains
    calcTotalStrainFn(&strainCell, deformCell, numQuadPts);

    if (!calcStress) 
      fieldSection->updatePoint(*c_iter, &strainCell[0]);
    else {
      _material->retrievePropsAndVars(*c_iter);
      stressCell = _material->calcStress(strainCell);
      fieldSection->updatePoint(*c_iter, &stressCell[0]);
    } // else
  } // for
} // _calcStrainStressField

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticityLgDeform::_calcStressFromStrain(
				   topology::Field<topology::Mesh>* field)
{ // _calcStressFromStrain
  assert(0 != field);
  assert(0 != _quadrature);
  assert(0 != _material);

  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  
  // Allocate arrays for cell data.
  double_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  double_array stressCell(numQuadPts*tensorSize);
  stressCell = 0.0;

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = field->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get field
  const ALE::Obj<RealSection>& fieldSection = field->section();
  assert(!fieldSection.isNull());

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    fieldSection->restrictPoint(*c_iter, &strainCell[0], strainCell.size());
    _material->retrievePropsAndVars(*c_iter);
    stressCell = _material->calcStress(strainCell);
    fieldSection->updatePoint(*c_iter, &stressCell[0]);
  } // for
} // _calcStressFromStrain

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 1-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityResidual1D(
				     const double_array& stress)
{ // _elasticityResidual1D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();

  assert(1 == cellDim);
  assert(quadWts.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad];
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const double N1 = wt*basisDeriv[iQuad*numBasis+iBasis  ];
      _cellVector[iBasis*spaceDim  ] -= N1*s11;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*5));
} // _elasticityResidual1D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 2-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityResidual2D(
				     const double_array& stress)
{ // _elasticityResidual2D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(2 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int stressSize = 3;

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad*stressSize+0];
    const double s22 = stress[iQuad*stressSize+1];
    const double s12 = stress[iQuad*stressSize+2];
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const double N1 = wt*basisDeriv[iQ+iBasis*spaceDim  ];
      const double N2 = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      _cellVector[iBasis*spaceDim  ] -= N1*s11 + N2*s12;
      _cellVector[iBasis*spaceDim+1] -= N1*s12 + N2*s22;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(8+2+9)));
} // _elasticityResidual2D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 3-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityResidual3D(
				     const double_array& stress)
{ // _elasticityResidual3D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(3 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int stressSize = 6;
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad*stressSize+0];
    const double s22 = stress[iQuad*stressSize+1];
    const double s33 = stress[iQuad*stressSize+2];
    const double s12 = stress[iQuad*stressSize+3];
    const double s23 = stress[iQuad*stressSize+4];
    const double s13 = stress[iQuad*stressSize+5];
    
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const int iBlock = iBasis*spaceDim;
      const double N1 = wt*basisDeriv[iQ+iBasis*spaceDim+0];
      const double N2 = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const double N3 = wt*basisDeriv[iQ+iBasis*spaceDim+2];

      _cellVector[iBlock  ] -= N1*s11 + N2*s12 + N3*s13;
      _cellVector[iBlock+1] -= N1*s12 + N2*s22 + N3*s23;
      _cellVector[iBlock+2] -= N1*s13 + N2*s23 + N3*s33;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(3+12)));
} // _elasticityResidual3D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 1-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityJacobian1D(
			       const double_array& elasticConsts)
{ // _elasticityJacobian1D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(1 == cellDim);
  assert(quadWts.size() == numQuadPts);
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double C1111 = elasticConsts[iQuad];
    for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
      const double valI = wt*basisDeriv[iQ+iBasis]*C1111;
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const double valIJ = valI * basisDeriv[iQ+jBasis];
	const int iBlock = iBasis*spaceDim * (numBasis*spaceDim);
	const int jBlock = jBasis*spaceDim;
	_cellMatrix[iBlock+jBlock] += valIJ;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*3)));
} // _elasticityJacobian1D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 2-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityJacobian2D(
			       const double_array& elasticConsts)
{ // _elasticityJacobian2D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(2 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int numConsts = 6;

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const double C1111 = elasticConsts[iQuad*numConsts+0];
    const double C1122 = elasticConsts[iQuad*numConsts+1];
    const double C1112 = elasticConsts[iQuad*numConsts+2]/2.0;
    const double C2222 = elasticConsts[iQuad*numConsts+3];
    const double C2212 = elasticConsts[iQuad*numConsts+4]/2.0;
    const double C1212 = elasticConsts[iQuad*numConsts+5]/2.0;
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const double Nip = wt*basisDeriv[iQ+iBasis*spaceDim  ];
      const double Niq = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const int iBlock = (iBasis*spaceDim  ) * (numBasis*spaceDim);
      const int iBlock1 = (iBasis*spaceDim+1) * (numBasis*spaceDim);
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const double Njp = basisDeriv[iQ+jBasis*spaceDim  ];
	const double Njq = basisDeriv[iQ+jBasis*spaceDim+1];
	const double ki0j0 = 
	  C1111 * Nip * Njp + C1112 * Niq * Njp +
	  C1112 * Nip * Njq + C1212 * Niq * Njq;
	const double ki0j1 =
	  C1122 * Nip * Njq + C2212 * Niq * Njq +
	  C1112 * Nip * Njp + C1212 * Niq * Njp;
	const double ki1j0 =
	  C1122 * Niq * Njp + C2212 * Niq * Njq +
	  C1112 * Nip * Njp + C1212 * Nip * Njq;
	const double ki1j1 =
	  C2222 * Niq * Njq + C2212 * Nip * Njq +
	  C2212 * Niq * Njp + C1212 * Nip * Njp;
	const int jBlock = (jBasis*spaceDim  );
	const int jBlock1 = (jBasis*spaceDim+1);
	_cellMatrix[iBlock +jBlock ] += ki0j0;
	_cellMatrix[iBlock +jBlock1] += ki0j1;
	_cellMatrix[iBlock1+jBlock ] += ki1j0;
	_cellMatrix[iBlock1+jBlock1] += ki1j1;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*(3*11+4))));
} // _elasticityJacobian2D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 3-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityJacobian3D(
			       const double_array& elasticConsts)
{ // _elasticityJacobian3D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(3 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int numConsts = 21;

  // Compute Jacobian for consistent tangent matrix
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const double C1111 = elasticConsts[iQuad*numConsts+ 0];
    const double C1122 = elasticConsts[iQuad*numConsts+ 1];
    const double C1133 = elasticConsts[iQuad*numConsts+ 2];
    const double C1112 = elasticConsts[iQuad*numConsts+ 3] / 2.0;
    const double C1123 = elasticConsts[iQuad*numConsts+ 4] / 2.0;
    const double C1113 = elasticConsts[iQuad*numConsts+ 5] / 2.0;
    const double C2222 = elasticConsts[iQuad*numConsts+ 6];
    const double C2233 = elasticConsts[iQuad*numConsts+ 7];
    const double C2212 = elasticConsts[iQuad*numConsts+ 8] / 2.0;
    const double C2223 = elasticConsts[iQuad*numConsts+ 9] / 2.0;
    const double C2213 = elasticConsts[iQuad*numConsts+10] / 2.0;
    const double C3333 = elasticConsts[iQuad*numConsts+11];
    const double C3312 = elasticConsts[iQuad*numConsts+12] / 2.0;
    const double C3323 = elasticConsts[iQuad*numConsts+13] / 2.0;
    const double C3313 = elasticConsts[iQuad*numConsts+14] / 2.0;
    const double C1212 = elasticConsts[iQuad*numConsts+15] / 2.0;
    const double C1223 = elasticConsts[iQuad*numConsts+16] / 2.0;
    const double C1213 = elasticConsts[iQuad*numConsts+17] / 2.0;
    const double C2323 = elasticConsts[iQuad*numConsts+18] / 2.0;
    const double C2313 = elasticConsts[iQuad*numConsts+19] / 2.0;
    const double C1313 = elasticConsts[iQuad*numConsts+20] / 2.0;
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const double Nip = wt*basisDeriv[iQ+iBasis*spaceDim+0];
      const double Niq = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const double Nir = wt*basisDeriv[iQ+iBasis*spaceDim+2];
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const double Njp = basisDeriv[iQ+jBasis*spaceDim+0];
	const double Njq = basisDeriv[iQ+jBasis*spaceDim+1];
	const double Njr = basisDeriv[iQ+jBasis*spaceDim+2];
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
	const double ki1j0 =
	  C1122 * Niq * Njp + C1112 * Nip * Njp + C1123 * Nir * Njp +
	  C2212 * Niq * Njq + C1212 * Nip * Njq + C1223 * Nir * Njq +
	  C2213 * Niq * Njr + C1213 * Nip * Njr + C2313 * Nir * Njr;
	const double ki1j1 =
	  C2222 * Niq * Njq + C2212 * Nip * Njq + C2223 * Nir * Njq +
	  C2212 * Niq * Njp + C1212 * Nip * Njp + C1223 * Nir * Njp +
	  C2223 * Niq * Njr + C1223 * Nip * Njr + C2323 * Nir * Njr;
	const double ki1j2 =
	  C2233 * Niq * Njr + C3312 * Nip * Njr + C3323 * Nir * Njr +
	  C2223 * Niq * Njq + C1223 * Nip * Njq + C2323 * Nir * Njq +
	  C2213 * Niq * Njp + C1213 * Nip * Njp + C2313 * Nir * Njp;
	const double ki2j0 =
	  C1133 * Nir * Njp + C1123 * Niq * Njp + C1113 * Nip * Njp +
	  C3312 * Nir * Njq + C1223 * Niq * Njq + C1213 * Nip * Njq +
	  C3313 * Nir * Njr + C2313 * Niq * Njr + C1313 * Nip * Njr; 
	const double ki2j1 =
	  C2233 * Nir * Njq + C2223 * Niq * Njq + C2213 * Nip * Njq +
	  C3312 * Nir * Njp + C1223 * Niq * Njp + C1213 * Nip * Njp +
	  C3323 * Nir * Njr + C2323 * Niq * Njr + C2313 * Nip * Njr; 
	const double ki2j2 =
	  C3333 * Nir * Njr + C3323 * Niq * Njr + C3313 * Nip * Njr +
	  C3323 * Nir * Njq + C2323 * Niq * Njq + C2313 * Nip * Njq +
	  C3313 * Nir * Njp + C2313 * Niq * Njp + C1313 * Nip * Njp;
	const int iBlock = iBasis*spaceDim * (numBasis*spaceDim);
	const int iBlock1 = (iBasis*spaceDim+1) * (numBasis*spaceDim);
	const int iBlock2 = (iBasis*spaceDim+2) * (numBasis*spaceDim);
	const int jBlock = jBasis*spaceDim;
	const int jBlock1 = jBasis*spaceDim+1;
	const int jBlock2 = jBasis*spaceDim+2;
	_cellMatrix[iBlock +jBlock ] += ki0j0;
	_cellMatrix[iBlock +jBlock1] += ki0j1;
	_cellMatrix[iBlock +jBlock2] += ki0j2;
	_cellMatrix[iBlock1+jBlock ] += ki1j0;
	_cellMatrix[iBlock1+jBlock1] += ki1j1;
	_cellMatrix[iBlock1+jBlock2] += ki1j2;
	_cellMatrix[iBlock2+jBlock ] += ki2j0;
	_cellMatrix[iBlock2+jBlock1] += ki2j1;
	_cellMatrix[iBlock2+jBlock2] += ki2j2;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(3+numBasis*(6*26+9))));
} // _elasticityJacobian3D

// ----------------------------------------------------------------------
// Calculate Green-Lagrange strain tensor at quadrature points of a 1-D cell.
void 
pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain1D(
					      double_array* strain,
					      const double_array& deform,
					      const int numQuadPts)
{ // _calcTotalStrain1D
  // Green-Lagrange strain tensor = 1/2 ( X^T X - I )
  // X: deformation tensor
  // I: identity matrix

  assert(0 != strain);

  const int dim = 1;
  const int strainSize = 1;
  assert(deform.size() == numQuadPts*dim*dim);
  assert(strain->size() == numQuadPts*strainSize);


  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    (*strain)[iQuad] = 0.5*(deform[iQuad]*deform[iQuad] - 1.0);
} // _calcTotalStrain1D
  
// ----------------------------------------------------------------------
// Calculate Green-Lagrange strain tensor at quadrature points of a 2-D cell.
void 
pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D(
					      double_array* strain,
					      const double_array& deform,
					      const int numQuadPts)
{ // _calcTotalStrain2D
  // Green-Lagrange strain tensor = 1/2 ( X^T X - I )
  // X: deformation tensor
  // I: identity matrix

  assert(0 != strain);

  const int dim = 2;
  const int deformSize = dim*dim;
  const int strainSize = 3;
  assert(deform.size() == numQuadPts*deformSize);
  assert(strain->size() == numQuadPts*strainSize);

  for (int iQuad=0, iDeform=0, iStrain=0;
       iQuad < numQuadPts;
       ++iQuad, iDeform+=deformSize, iStrain+=strainSize) {
    (*strain)[iStrain  ] =
      0.5 * (deform[iDeform  ]*deform[iDeform  ] + 
	     deform[iDeform+2]*deform[iDeform+2] - 1.0);
    (*strain)[iStrain+1] =
      0.5 * (deform[iDeform+1]*deform[iDeform+1] + 
	     deform[iDeform+3]*deform[iDeform+3] - 1.0);
    (*strain)[iStrain+2] =
      0.5 * (deform[iDeform  ]*deform[iDeform+1] + 
	     deform[iDeform+2]*deform[iDeform+3]);
  } // for
} // _calcTotalStrain2D
  
// ----------------------------------------------------------------------
// Calculate Green-Lagrange strain tensor at quadrature points of a 3-D cell.
void 
pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D(
					      double_array* strain,
					      const double_array& deform,
					      const int numQuadPts)
{ // _calcTotalStrain3D
  // Green-Lagrange strain tensor = 1/2 ( X^T X - I )
  // X: deformation tensor
  // I: identity matrix

  assert(0 != strain);

  const int dim = 3;
  const int deformSize = dim*dim;
  const int strainSize = 6;
  assert(deform.size() == numQuadPts*dim*dim);
  assert(strain->size() == numQuadPts*strainSize);

  for (int iQuad=0, iDeform=0, iStrain=0;
       iQuad < numQuadPts;
       ++iQuad, iDeform+=deformSize, iStrain+=strainSize) {
    (*strain)[iStrain  ] =
      0.5 * (deform[iDeform  ]*deform[iDeform  ] +
	     deform[iDeform+3]*deform[iDeform+3] +
	     deform[iDeform+6]*deform[iDeform+6] - 1.0);
    (*strain)[iStrain+1] =
      0.5 * (deform[iDeform+1]*deform[iDeform+1] +
	     deform[iDeform+4]*deform[iDeform+4] +
	     deform[iDeform+7]*deform[iDeform+7] - 1.0);
    (*strain)[iStrain+2] =
      0.5 * (deform[iDeform+2]*deform[iDeform+2] +
	     deform[iDeform+5]*deform[iDeform+5] +
	     deform[iDeform+8]*deform[iDeform+8] - 1.0);
    (*strain)[iStrain+3] =
      0.5 * (deform[iDeform  ]*deform[iDeform+1] +
	     deform[iDeform+3]*deform[iDeform+4] +
	     deform[iDeform+6]*deform[iDeform+7]);
    (*strain)[iStrain+4] =
      0.5 * (deform[iDeform+1]*deform[iDeform+2] +
	     deform[iDeform+4]*deform[iDeform+5] +
	     deform[iDeform+7]*deform[iDeform+8]);
    (*strain)[iStrain+5] =
      0.5 * (deform[iDeform+0]*deform[iDeform+2] +
	     deform[iDeform+3]*deform[iDeform+5] +
	     deform[iDeform+6]*deform[iDeform+8]);
  } // for
} // _calcTotalStrain3D
  
// ----------------------------------------------------------------------
// Calculate deformation tensor.
void 
pylith::feassemble::IntegratorElasticityLgDeform::_calcDeformation(
					      double_array* deform,
					      const double_array& basisDeriv,
					      const double_array& vertices,
					      const double_array& disp,
					      const int numBasis,
					      const int numQuadPts,
					      const int dim)
{ // _calcDeformation
  assert(0 != deform);

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  const int deformSize = dim*dim;

  (*deform) = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis)
      for (int iDim=0, indexD=0; iDim < dim; ++iDim)
	for (int jDim=0; jDim < dim; ++jDim, ++indexD)
	  (*deform)[iQuad*deformSize+indexD] += 
	    basisDeriv[iQ+iBasis*dim+jDim] *
	    (vertices[iBasis*dim+iDim] + disp[iBasis*dim+iDim]);
} // _calcDeformation
  

// End of file 
