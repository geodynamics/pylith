// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

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
  PYLITH_METHOD_BEGIN;

  _needNewJacobian = IntegratorElasticity::needNewJacobian();
  PYLITH_METHOD_RETURN(_needNewJacobian);
} // needNewJacobian

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorElasticityLgDeform::updateStateVars(const PylithScalar t,
								  topology::SolutionFields* const fields)
{ // updateStateVars
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_material);
  assert(fields);

  // No need to update state vars if material doesn't have any.
  if (!_material->hasStateVars())
    PYLITH_METHOD_END;

  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (2 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D;
  } else {
    std::ostringstream msg;
    msg << "Bad cell dimension '" << cellDim << " in IntegratorElasticityLgDeform::updateStateVars()'." << std::endl;
    throw std::logic_error(msg.str());
  } // else

  // Allocate arrays for cell data.
  scalar_array dispCellTmp(numBasis*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);
  scalar_array deformCell(numQuadPts*spaceDim*spaceDim);
  deformCell = 0.0;
  strainCell = 0.0;

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  scalar_array dispCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"), "displacement");
  dispVisitor.optimizeClosure();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  _material->createPropsAndVarsVisitors();

  // Loop over cells
  for (PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Retrieve geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);
    const scalar_array& basisDeriv = _quadrature->basisDeriv();

    // Get physical properties and state variables for cell.
    _material->retrievePropsAndVars(cell);

    dispVisitor.getClosure(&dispCell, cell);
  
    // Compute deformation tensor.
    _calcDeformation(&deformCell, basisDeriv, &dispCell[0], numBasis, numQuadPts, spaceDim);

    // Compute strains
    calcTotalStrainFn(&strainCell, deformCell, numQuadPts);

    // Update material state
    _material->updateStateVars(strainCell, cell);
  } // for
  _material->destroyPropsAndVarsVisitors();

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticityLgDeform::_calcStrainStressField(topology::Field* field,
									 const char* name,
									 topology::SolutionFields* const fields)
{ // _calcStrainStressField
  PYLITH_METHOD_BEGIN;

  assert(field);
  assert(_quadrature);
  assert(_material);

  const bool calcStress = (0 == strcasecmp(name, "stress") || 0 == strcasecmp(name, "cauchy_stress")) ? true : false;
  const bool calcCauchyStress = (0 == strcasecmp(name, "cauchy_stress")) ? true : false;

    
  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (2 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D;
  } else {
    std::ostringstream msg;
    msg << "Bad cell dimension '" << cellDim << " in IntegratorElasticityLgDeform::_calcStrainStressField()'." << std::endl;
    throw std::logic_error(msg.str());
  } // else
  
  // Allocate arrays for cell data.
  scalar_array deformCell(numQuadPts*spaceDim*spaceDim);
  const int tensorCellSize = numQuadPts*tensorSize;
  scalar_array strainCell(tensorCellSize);
  scalar_array stressCell(tensorCellSize);
  scalar_array stressCauchyCell(tensorCellSize);
  deformCell = 0.0;
  strainCell = 0.0;
  stressCell = 0.0;
  stressCauchyCell = 0.0;

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Setup field visitors.
  scalar_array dispCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"), "displacement");
  dispVisitor.optimizeClosure();

  topology::VecVisitorMesh fieldVisitor(*field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  _material->createPropsAndVarsVisitors();

  // Loop over cells
  for (PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Retrieve geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);
    const scalar_array& basisDeriv = _quadrature->basisDeriv();

    // Restrict input fields to cell
    dispVisitor.getClosure(&dispCell, cell);

    // Compute deformation tensor.
    _calcDeformation(&deformCell, basisDeriv, &dispCell[0], numBasis, numQuadPts, spaceDim);

    // Compute strains
    calcTotalStrainFn(&strainCell, deformCell, numQuadPts);

    const PetscInt off = fieldVisitor.sectionOffset(cell);
    assert(tensorCellSize == fieldVisitor.sectionDof(cell));
    if (!calcStress) {
      for (int i=0; i < tensorCellSize; ++i) {
	fieldArray[off+i] = strainCell[i];
      } // for
    } else {
      _material->retrievePropsAndVars(cell);
      stressCell = _material->calcStress(strainCell);

      if (calcCauchyStress) {
	if (2 == spaceDim) {
	  _calcCauchyStress2D(&stressCauchyCell, stressCell, deformCell, numQuadPts);
	} else if (3 == spaceDim) {
	  _calcCauchyStress3D(&stressCauchyCell, stressCell, deformCell, numQuadPts);
	} // if/else
	for (int i=0; i < tensorCellSize; ++i) {
	  fieldArray[off+i] = stressCauchyCell[i];
	} // for
      } else {
	for (int i=0; i < tensorCellSize; ++i) {
	  fieldArray[off+i] = stressCell[i];
	} // for
      } // if/else
    } // else
  } // for
  _material->destroyPropsAndVarsVisitors();

  PYLITH_METHOD_END;
} // _calcStrainStressField

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 2-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityResidual2D(const scalar_array& stress,
									const scalar_array& disp)
{ // _elasticityResidual2D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(2 == cellDim);
  assert(quadWts.size() == size_t(numQuadPts));
  const int stressSize = 3;

#if 1
  // Compute deformation gradient tensor.
  scalar_array deform(numQuadPts*spaceDim*spaceDim);
  _calcDeformation(&deform, basisDeriv, &disp[0], numBasis, numQuadPts, spaceDim);
#endif

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const int iQ = iQuad*numBasis*spaceDim;
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    const PylithScalar s11 = stress[iQuad*stressSize+0];
    const PylithScalar s22 = stress[iQuad*stressSize+1];
    const PylithScalar s12 = stress[iQuad*stressSize+2];


#if 1
    // Convert second Piola-Kirchoff stress, S, to first Piola-Kirchoff stress, P
    // 
    // P = S X^t, Pij = Sik Xjk
    const int iD = iQuad*spaceDim*spaceDim;
    const PylithScalar p11 = s11*deform[iD  ] + s12*deform[iD+1];
    const PylithScalar p12 = s11*deform[iD+2] + s12*deform[iD+3];
    const PylithScalar p21 = s12*deform[iD  ] + s22*deform[iD+1];
    const PylithScalar p22 = s12*deform[iD+2] + s22*deform[iD+3];

    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      const PylithScalar Nip = wt*basisDeriv[iQ+iB  ];
      const PylithScalar Niq = wt*basisDeriv[iQ+iB+1];
    
      _cellVector[iB  ] -= Nip*p11 + Niq*p21;
      _cellVector[iB+1] -= Nip*p12 + Niq*p22;
#else


    PylithScalar l11 = 0.0;
    PylithScalar l12 = 0.0;
    PylithScalar l21 = 0.0;
    PylithScalar l22 = 0.0;
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      l11 += basisDeriv[iQ+iB  ] * disp[iB  ];
      l12 += basisDeriv[iQ+iB+1] * disp[iB  ];
      l21 += basisDeriv[iQ+iB  ] * disp[iB+1];
      l22 += basisDeriv[iQ+iB+1] * disp[iB+1];
    } // for

    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      const PylithScalar Nip = wt*basisDeriv[iQ+iB  ];
      const PylithScalar Niq = wt*basisDeriv[iQ+iB+1];

      // Generated using Maxima (see jacobian2d_lgdeform.wxm)
      _cellVector[iB  ] -= (1+l11)*Nip*s11 + l12*Niq*s22 + ((1+l11)*Niq + l12*Nip)*s12;
      _cellVector[iB+1] -= l21*Nip*s11 + (1+l22)*Niq*s22 + (l21*Niq + (1+l22)*Nip)*s12;
#endif
    } // for
  } // for

  PetscLogFlops(numQuadPts*(1+numBasis*(4*2+2+2*13)));
} // _elasticityResidual2D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 3-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityResidual3D(const scalar_array& stress,
									const scalar_array& disp)
{ // _elasticityResidual3D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(3 == cellDim);
  assert(quadWts.size() == size_t(numQuadPts));
  const int stressSize = 6;
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const int iQ = iQuad*numBasis*spaceDim;
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    const PylithScalar s11 = stress[iQuad*stressSize+0];
    const PylithScalar s22 = stress[iQuad*stressSize+1];
    const PylithScalar s33 = stress[iQuad*stressSize+2];
    const PylithScalar s12 = stress[iQuad*stressSize+3];
    const PylithScalar s23 = stress[iQuad*stressSize+4];
    const PylithScalar s13 = stress[iQuad*stressSize+5];
    
    PylithScalar l11 = 0.0;
    PylithScalar l12 = 0.0;
    PylithScalar l13 = 0.0;
    PylithScalar l21 = 0.0;
    PylithScalar l22 = 0.0;
    PylithScalar l23 = 0.0;
    PylithScalar l31 = 0.0;
    PylithScalar l32 = 0.0;
    PylithScalar l33 = 0.0;
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      l11 += basisDeriv[iQ+iB  ] * disp[iB  ];
      l12 += basisDeriv[iQ+iB+1] * disp[iB  ];
      l13 += basisDeriv[iQ+iB+2] * disp[iB  ];
      l21 += basisDeriv[iQ+iB  ] * disp[iB+1];
      l22 += basisDeriv[iQ+iB+1] * disp[iB+1];
      l23 += basisDeriv[iQ+iB+2] * disp[iB+1];
      l31 += basisDeriv[iQ+iB  ] * disp[iB+2];
      l32 += basisDeriv[iQ+iB+1] * disp[iB+2];
      l33 += basisDeriv[iQ+iB+2] * disp[iB+2];
    } // for
    
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      const PylithScalar Nip = wt*basisDeriv[iQ+iB  ];
      const PylithScalar Niq = wt*basisDeriv[iQ+iB+1];
      const PylithScalar Nir = wt*basisDeriv[iQ+iB+2];

      // Generated using Maxima (see jacobian3d_lgdeform.wxm)
      _cellVector[iB  ] -= l13*Nir*s33 + (l12*Nir+l13*Niq)*s23 + l12*Niq*s22 + ((l11+1)*Nir + l13*Nip)*s13 + ((l11+1)*Niq+l12*Nip)*s12 + (l11+1)*Nip*s11;
      _cellVector[iB+1] -= l23*Nir*s33 + ((l22+1)*Nir+l23*Niq)*s23 + (l22+1)*Niq*s22 + (l21*Nir+l23*Nip)*s13 + (l21*Niq+(l22+1)*Nip)*s12 + l21*Nip*s11;
      _cellVector[iB+2] -= (l33+1)*Nir*s33 + (l32*Nir+(l33+1)*Niq)*s23 + l32*Niq*s22 + (l31*Nir+(l33+1)*Nip)*s13 + (l31*Niq+l32*Nip)*s12 + l31*Nip*s11;
    } // for
  } // for

  PetscLogFlops(numQuadPts*(1+numBasis*(2*9+3+3*27)));
} // _elasticityResidual3D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 2-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityJacobian2D(const scalar_array& elasticConsts,
									const scalar_array& stress,
									const scalar_array& disp)
{ // _elasticityJacobian2D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  const int tensorSize = _material->tensorSize();
  
  assert(2 == cellDim);
  assert(quadWts.size() == size_t(numQuadPts));
  const int numConsts = 9;

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const int iQ = iQuad*numBasis*spaceDim;
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const int iC = iQuad*numConsts;
    const PylithScalar C1111 = elasticConsts[iC+0];
    const PylithScalar C1122 = elasticConsts[iC+1];
    const PylithScalar C1112 = elasticConsts[iC+2] / 2.0; // 2*mu -> mu
    const PylithScalar C2211 = elasticConsts[iC+3];
    const PylithScalar C2222 = elasticConsts[iC+4];
    const PylithScalar C2212 = elasticConsts[iC+5] / 2.0;
    const PylithScalar C1211 = elasticConsts[iC+6];
    const PylithScalar C1222 = elasticConsts[iC+7];
    const PylithScalar C1212 = elasticConsts[iC+8] / 2.0;

    const int iS = iQuad*tensorSize;
    const PylithScalar s11 = stress[iS+0];
    const PylithScalar s22 = stress[iS+1];
    const PylithScalar s12 = stress[iS+2];

    PylithScalar l11 = 0.0;
    PylithScalar l12 = 0.0;
    PylithScalar l21 = 0.0;
    PylithScalar l22 = 0.0;
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      l11 += basisDeriv[iQ+iB  ] * disp[iB  ];
      l12 += basisDeriv[iQ+iB+1] * disp[iB  ];
      l21 += basisDeriv[iQ+iB  ] * disp[iB+1];
      l22 += basisDeriv[iQ+iB+1] * disp[iB+1];
    } // for

    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      const PylithScalar Nip = wt*basisDeriv[iQ+iB  ];
      const PylithScalar Niq = wt*basisDeriv[iQ+iB+1];

      const int iBlock  = (iB  ) * (numBasis*spaceDim);
      const int iBlock1 = (iB+1) * (numBasis*spaceDim);

      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const int jB = jBasis*spaceDim;
	const PylithScalar Njp = basisDeriv[iQ+jB  ];
	const PylithScalar Njq = basisDeriv[iQ+jB+1];

	const PylithScalar Ki0j0 = 
	  l12*Niq*(l12*Njq*C2222 + 
		   ((l11+1)*Njq+l12*Njp)*C2212 + 
		   (l11+1)*Njp*C2211) + 
	  ((l11+1)*Niq+l12*Nip)*(l12*Njq*C1222 + 
				 ((l11+1)*Njq+l12*Njp)*C1212 + 
				 (l11+1)*Njp*C1211) + 
	  (l11+1)*Nip*(l12*Njq*C1122 + 
		       ((l11+1)*Njq+l12*Njp)*C1112 + 
		       (l11+1)*Njp*C1111);

	const PylithScalar Ki0j1 =
	  l12*Niq*((l22+1.0)*Njq*C2222 + 
		   (l21*Njq+(l22+1.0)*Njp)*C2212 + 
		   l21*Njp*C2211) + 
	  ((l11+1.0)*Niq+l12*Nip)*((l22+1.0)*Njq*C1222 + 
				 (l21*Njq+(l22+1.0)*Njp)*C1212 + 
				 l21*Njp*C1211) + 
	  (l11+1.0)*Nip*((l22+1.0)*Njq*C1122 + 
		       (l21*Njq+(l22+1.0)*Njp)*C1112 + 
		       l21*Njp*C1111);

	const PylithScalar Ki1j0 =
	  (l22+1.0)*Niq*(l12*Njq*C2222 + 
		       ((l11+1.0)*Njq+l12*Njp)*C2212 + 
		       (l11+1.0)*Njp*C2211) + 
	  (l21*Niq+(l22+1.0)*Nip)*(l12*Njq*C1222 + 
				 ((l11+1.0)*Njq+l12*Njp)*C1212 + 
				 (l11+1.0)*Njp*C1211) + 
	  l21*Nip*(l12*Njq*C1122 + 
		   ((l11+1.0)*Njq+l12*Njp)*C1112 + 
		   (l11+1.0)*Njp*C1111);

	const PylithScalar Ki1j1 =
	  (l22+1.0)*Niq*((l22+1.0)*Njq*C2222 + 
		       (l21*Njq+(l22+1.0)*Njp)*C2212 + 
		       l21*Njp*C2211) + 
	  (l21*Niq+(l22+1.0)*Nip)*((l22+1.0)*Njq*C1222 + 
				 (l21*Njq+(l22+1.0)*Njp)*C1212 + 
				 l21*Njp*C1211) + 
	  l21*Nip*((l22+1.0)*Njq*C1122 + 
		   (l21*Njq+(l22+1.0)*Njp)*C1112 + 
		   l21*Njp*C1111);

#if 0
	const PylithScalar Knl = 0.0;
#else
	const PylithScalar Knl = (Njp*s11 + Njq*s12)*Nip + (Njp*s12 + Njq*s22)*Niq;
#endif

	const int jBlock = (jB);
	const int jBlock1 = (jB+1);
	_cellMatrix[iBlock +jBlock ] += Ki0j0 + Knl;
	_cellMatrix[iBlock +jBlock1] += Ki0j1;
	_cellMatrix[iBlock1+jBlock ] += Ki1j0;
	_cellMatrix[iBlock1+jBlock1] += Ki1j1 + Knl;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*(3*11+4))));
} // _elasticityJacobian2D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 3-D cells.
void
pylith::feassemble::IntegratorElasticityLgDeform::_elasticityJacobian3D(const scalar_array& elasticConsts,
									const scalar_array& stress,
									const scalar_array& disp)
{ // _elasticityJacobian3D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  const int tensorSize = _material->tensorSize();

  assert(3 == cellDim);
  assert(quadWts.size() == size_t(numQuadPts));
  assert(6 == tensorSize);
  const int numConsts = 36;

  // Compute Jacobian for consistent tangent matrix
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const int iQ = iQuad*numBasis*spaceDim;
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const int iC = iQuad*numConsts;
    const PylithScalar C1111 = elasticConsts[iC+ 0];
    const PylithScalar C1122 = elasticConsts[iC+ 1];
    const PylithScalar C1133 = elasticConsts[iC+ 2];
    const PylithScalar C1112 = elasticConsts[iC+ 3] / 2.0;
    const PylithScalar C1123 = elasticConsts[iC+ 4] / 2.0;
    const PylithScalar C1113 = elasticConsts[iC+ 5] / 2.0;
    const PylithScalar C2211 = elasticConsts[iC+ 6];
    const PylithScalar C2222 = elasticConsts[iC+ 7];
    const PylithScalar C2233 = elasticConsts[iC+ 8];
    const PylithScalar C2212 = elasticConsts[iC+ 9] / 2.0;
    const PylithScalar C2223 = elasticConsts[iC+10] / 2.0;
    const PylithScalar C2213 = elasticConsts[iC+11] / 2.0;
    const PylithScalar C3311 = elasticConsts[iC+12];
    const PylithScalar C3322 = elasticConsts[iC+13];
    const PylithScalar C3333 = elasticConsts[iC+14];
    const PylithScalar C3312 = elasticConsts[iC+15] / 2.0;
    const PylithScalar C3323 = elasticConsts[iC+16] / 2.0;
    const PylithScalar C3313 = elasticConsts[iC+17] / 2.0;
    const PylithScalar C1211 = elasticConsts[iC+18];
    const PylithScalar C1222 = elasticConsts[iC+19];
    const PylithScalar C1233 = elasticConsts[iC+20];
    const PylithScalar C1212 = elasticConsts[iC+21] / 2.0;
    const PylithScalar C1223 = elasticConsts[iC+22] / 2.0;
    const PylithScalar C1213 = elasticConsts[iC+23] / 2.0;
    const PylithScalar C2311 = elasticConsts[iC+24];
    const PylithScalar C2322 = elasticConsts[iC+25];
    const PylithScalar C2333 = elasticConsts[iC+26];
    const PylithScalar C2312 = elasticConsts[iC+27] / 2.0;
    const PylithScalar C2323 = elasticConsts[iC+28] / 2.0;
    const PylithScalar C2313 = elasticConsts[iC+29] / 2.0;
    const PylithScalar C1311 = elasticConsts[iC+30];
    const PylithScalar C1322 = elasticConsts[iC+31];
    const PylithScalar C1333 = elasticConsts[iC+32];
    const PylithScalar C1312 = elasticConsts[iC+33] / 2.0;
    const PylithScalar C1323 = elasticConsts[iC+34] / 2.0;
    const PylithScalar C1313 = elasticConsts[iC+35] / 2.0;

    const int iS = iQuad*tensorSize;
    const PylithScalar s11 = stress[iS+0];
    const PylithScalar s22 = stress[iS+1];
    const PylithScalar s33 = stress[iS+2];
    const PylithScalar s12 = stress[iS+3];
    const PylithScalar s23 = stress[iS+4];
    const PylithScalar s13 = stress[iS+5];

    PylithScalar l11 = 0.0;
    PylithScalar l12 = 0.0;
    PylithScalar l13 = 0.0;
    PylithScalar l21 = 0.0;
    PylithScalar l22 = 0.0;
    PylithScalar l23 = 0.0;
    PylithScalar l31 = 0.0;
    PylithScalar l32 = 0.0;
    PylithScalar l33 = 0.0;
    for (int kBasis=0; kBasis < numBasis; ++kBasis) {
      const int kB = kBasis*spaceDim;
      l11 += basisDeriv[iQ+kB  ] * disp[kB  ];
      l12 += basisDeriv[iQ+kB+1] * disp[kB  ];
      l13 += basisDeriv[iQ+kB+2] * disp[kB  ];
      l21 += basisDeriv[iQ+kB  ] * disp[kB+1];
      l22 += basisDeriv[iQ+kB+1] * disp[kB+1];
      l23 += basisDeriv[iQ+kB+2] * disp[kB+1];
      l31 += basisDeriv[iQ+kB  ] * disp[kB+2];
      l32 += basisDeriv[iQ+kB+1] * disp[kB+2];
      l33 += basisDeriv[iQ+kB+2] * disp[kB+2];
    } // for
    
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const int iB = iBasis*spaceDim;
      const PylithScalar Nip = wt*basisDeriv[iQ+iB+0];
      const PylithScalar Niq = wt*basisDeriv[iQ+iB+1];
      const PylithScalar Nir = wt*basisDeriv[iQ+iB+2];
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const int jB = jBasis*spaceDim;
	const PylithScalar Njp = basisDeriv[iQ+jB+0];
	const PylithScalar Njq = basisDeriv[iQ+jB+1];
	const PylithScalar Njr = basisDeriv[iQ+jB+2];



	// Generated using Maxima (see jacobian3d_lgdeform.wxm)
	const PylithScalar Ki0j0 = 
	  l13*Nir*(l13*Njr*C3333 + 
		   (l12*Njr+l13*Njq)*C3323 + 
		   ((l11+1)*Njr+l13*Njp)*C3313 + 
		   ((l11+1)*Njq+l12*Njp)*C3312 + 
		   l12*Njq*C3322 + 
		   (l11+1)*Njp*C3311) + 
	  (l12*Nir+l13*Niq)*(l13*Njr*C2333 + 
			     (l12*Njr+l13*Njq)*C2323 + 
			     ((l11+1)*Njr+l13*Njp)*C2313 + 
			     l12*Njq*C2322 + 
			     ((l11+1)*Njq+l12*Njp)*C2312 + 
			     (l11+1)*Njp*C2311) + 
	  ((l11+1)*Nir+l13*Nip)*(l13*Njr*C1333 + 
				 (l12*Njr+l13*Njq)*C1323 + 
				 l12*Njq*C1322 + 
				 ((l11+1)*Njr+l13*Njp)*C1313 + 
				 ((l11+1)*Njq+l12*Njp)*C1312 + 
				 (l11+1)*Njp*C1311) + 
	  ((l11+1)*Niq+l12*Nip)*(l13*Njr*C1233 + 
				 l12*Njq*C1222 + 
				 (l12*Njr+l13*Njq)*C1223 + 
				 ((l11+1)*Njr+l13*Njp)*C1213 + 
				 ((l11+1)*Njq+l12*Njp)*C1212 + 
				 (l11+1)*Njp*C1211) + 
	  l12*Niq*(l13*Njr*C2233 + 
		   (l12*Njr+l13*Njq)*C2223 + 
		   l12*Njq*C2222 + 
		   ((l11+1)*Njr+l13*Njp)*C2213 + 
		   ((l11+1)*Njq+l12*Njp)*C2212 + 
		   (l11+1)*Njp*C2211) + 
	  (l11+1)*Nip*(l13*Njr*C1133 + 
		       (l12*Njr+l13*Njq)*C1123 + 
		       l12*Njq*C1122 + 
		       ((l11+1)*Njr+l13*Njp)*C1113 + 
		       ((l11+1)*Njq+l12*Njp)*C1112 + 
		       (l11+1)*Njp*C1111);

	const PylithScalar Ki0j1 =
	  l13*Nir*(l23*Njr*C3333 + 
		   ((l22+1)*Njr+l23*Njq)*C3323 + 
		   (l21*Njr+l23*Njp)*C3313 + 
		   (l21*Njq+(l22+1)*Njp)*C3312 + 
		   (l22+1)*Njq*C3322 + 
		   l21*Njp*C3311) + 
	  (l12*Nir+l13*Niq)*(l23*Njr*C2333 + 
			     ((l22+1)*Njr+l23*Njq)*C2323 + 
			     (l21*Njr+l23*Njp)*C2313 + 
			     (l22+1)*Njq*C2322 + 
			     (l21*Njq+(l22+1)*Njp)*C2312 + 
			     l21*Njp*C2311) + 
	  ((l11+1)*Nir+l13*Nip)*(l23*Njr*C1333 + 
				 ((l22+1)*Njr+l23*Njq)*C1323 + 
				 (l22+1)*Njq*C1322 + 
				 (l21*Njr+l23*Njp)*C1313 + 
				 (l21*Njq+(l22+1)*Njp)*C1312 + 
				 l21*Njp*C1311) + 
	  ((l11+1)*Niq+l12*Nip)*(l23*Njr*C1233 + 
				 (l22+1)*Njq*C1222 + 
				 ((l22+1)*Njr+l23*Njq)*C1223 + 
				 (l21*Njr+l23*Njp)*C1213 + 
				 (l21*Njq+(l22+1)*Njp)*C1212 + 
				 l21*Njp*C1211) + 
	  l12*Niq*(l23*Njr*C2233 + 
		   ((l22+1)*Njr+l23*Njq)*C2223 + 
		   (l22+1)*Njq*C2222 + 
		   (l21*Njr+l23*Njp)*C2213 + 
		   (l21*Njq+(l22+1)*Njp)*C2212 + 
		   l21*Njp*C2211) + 
	  (l11+1)*Nip*(l23*Njr*C1133 + 
		       ((l22+1)*Njr+l23*Njq)*C1123 + 
		       (l22+1)*Njq*C1122 + 
		       (l21*Njr+l23*Njp)*C1113 + 
		       (l21*Njq+(l22+1)*Njp)*C1112 + 
		       l21*Njp*C1111);

	const PylithScalar Ki0j2 =
	  l13*Nir*((l33+1)*Njr*C3333 + 
		   (l32*Njr+(l33+1)*Njq)*C3323 + 
		   (l31*Njr+(l33+1)*Njp)*C3313 + 
		   (l31*Njq+l32*Njp)*C3312 + 
		   l32*Njq*C3322 + 
		   l31*Njp*C3311) + 
	  (l12*Nir+l13*Niq)*((l33+1)*Njr*C2333 + 
			     (l32*Njr+(l33+1)*Njq)*C2323 + 
			     (l31*Njr+(l33+1)*Njp)*C2313 + 
			     l32*Njq*C2322 + 
			     (l31*Njq+l32*Njp)*C2312 + 
			     l31*Njp*C2311) + 
	  ((l11+1)*Nir+l13*Nip)*((l33+1)*Njr*C1333 + 
				 (l32*Njr+(l33+1)*Njq)*C1323 + 
				 l32*Njq*C1322 + 
				 (l31*Njr+(l33+1)*Njp)*C1313 + 
				 (l31*Njq+l32*Njp)*C1312 + 
				 l31*Njp*C1311) + 
	  ((l11+1)*Niq+l12*Nip)*((l33+1)*Njr*C1233 + 
				 l32*Njq*C1222 + 
				 (l32*Njr+(l33+1)*Njq)*C1223 + 
				 (l31*Njr+(l33+1)*Njp)*C1213 + 
				 (l31*Njq+l32*Njp)*C1212 + 
				 l31*Njp*C1211) + 
	  l12*Niq*((l33+1)*Njr*C2233 + 
		   (l32*Njr+(l33+1)*Njq)*C2223 + 
		   l32*Njq*C2222 + 
		   (l31*Njr+(l33+1)*Njp)*C2213 + 
		   (l31*Njq+l32*Njp)*C2212 + 
		   l31*Njp*C2211) + 
	  (l11+1)*Nip*((l33+1)*Njr*C1133 + 
		       (l32*Njr+(l33+1)*Njq)*C1123 + 
		       l32*Njq*C1122 + 
		       (l31*Njr+(l33+1)*Njp)*C1113 + 
		       (l31*Njq+l32*Njp)*C1112 + 
		       l31*Njp*C1111);

	const PylithScalar Ki1j0 =
	  l23*Nir*(l13*Njr*C3333 + 
		   (l12*Njr+l13*Njq)*C3323 + 
		   ((l11+1)*Njr+l13*Njp)*C3313 + 
		   ((l11+1)*Njq+l12*Njp)*C3312 + 
		   l12*Njq*C3322 + 
		   (l11+1)*Njp*C3311) + 
	  ((l22+1)*Nir+l23*Niq)*(l13*Njr*C2333 + 
				 (l12*Njr+l13*Njq)*C2323 + 
				 ((l11+1)*Njr+l13*Njp)*C2313 + 
				 l12*Njq*C2322 + 
				 ((l11+1)*Njq+l12*Njp)*C2312 + 
				 (l11+1)*Njp*C2311) + 
	  (l21*Nir+l23*Nip)*(l13*Njr*C1333 + 
			     (l12*Njr+l13*Njq)*C1323 + 
			     l12*Njq*C1322 + 
			     ((l11+1)*Njr+l13*Njp)*C1313 + 
			     ((l11+1)*Njq+l12*Njp)*C1312 + 
			     (l11+1)*Njp*C1311) + 
	  (l21*Niq+(l22+1)*Nip)*(l13*Njr*C1233 + 
				 l12*Njq*C1222 + 
				 (l12*Njr+l13*Njq)*C1223 + 
				 ((l11+1)*Njr+l13*Njp)*C1213 + 
				 ((l11+1)*Njq+l12*Njp)*C1212 + 
				 (l11+1)*Njp*C1211) + 
	  (l22+1)*Niq*(l13*Njr*C2233 + 
		       (l12*Njr+l13*Njq)*C2223 + 
		       l12*Njq*C2222 + 
		       ((l11+1)*Njr+l13*Njp)*C2213 + 
		       ((l11+1)*Njq+l12*Njp)*C2212 + 
		       (l11+1)*Njp*C2211) + 
	  l21*Nip*(l13*Njr*C1133 + 
		   (l12*Njr+l13*Njq)*C1123 + 
		   l12*Njq*C1122 + 
		   ((l11+1)*Njr+l13*Njp)*C1113 + 
		   ((l11+1)*Njq+l12*Njp)*C1112 + 
		   (l11+1)*Njp*C1111);

	const PylithScalar Ki1j1 =
	  l23*Nir*(l23*Njr*C3333 + 
		   ((l22+1)*Njr+l23*Njq)*C3323 + 
		   (l21*Njr+l23*Njp)*C3313 + 
		   (l21*Njq+(l22+1)*Njp)*C3312 + 
		   (l22+1)*Njq*C3322 + 
		   l21*Njp*C3311) + 
	  ((l22+1)*Nir+l23*Niq)*(l23*Njr*C2333 + 
				 ((l22+1)*Njr+l23*Njq)*C2323 + 
				 (l21*Njr+l23*Njp)*C2313 + 
				 (l22+1)*Njq*C2322 + 
				 (l21*Njq+(l22+1)*Njp)*C2312 + 
				 l21*Njp*C2311) + 
	  (l21*Nir+l23*Nip)*(l23*Njr*C1333 + 
			     ((l22+1)*Njr+l23*Njq)*C1323 + 
			     (l22+1)*Njq*C1322 + 
			     (l21*Njr+l23*Njp)*C1313 + 
			     (l21*Njq+(l22+1)*Njp)*C1312 + 
			     l21*Njp*C1311) + 
	  (l21*Niq+(l22+1)*Nip)*(l23*Njr*C1233 + 
				 (l22+1)*Njq*C1222 + 
				 ((l22+1)*Njr+l23*Njq)*C1223 + 
				 (l21*Njr+l23*Njp)*C1213 + 
				 (l21*Njq+(l22+1)*Njp)*C1212 + 
				 l21*Njp*C1211) + 
	  (l22+1)*Niq*(l23*Njr*C2233 + 
		       ((l22+1)*Njr+l23*Njq)*C2223 + 
		       (l22+1)*Njq*C2222 + 
		       (l21*Njr+l23*Njp)*C2213 + 
		       (l21*Njq+(l22+1)*Njp)*C2212 + 
		       l21*Njp*C2211) + 
	  l21*Nip*(l23*Njr*C1133 + 
		   ((l22+1)*Njr+l23*Njq)*C1123 + 
		   (l22+1)*Njq*C1122 + 
		   (l21*Njr+l23*Njp)*C1113 + 
		   (l21*Njq+(l22+1)*Njp)*C1112 + 
		   l21*Njp*C1111);

	const PylithScalar Ki1j2 =
	  l23*Nir*((l33+1)*Njr*C3333 + 
		   (l32*Njr+(l33+1)*Njq)*C3323 + 
		   (l31*Njr+(l33+1)*Njp)*C3313 + 
		   (l31*Njq+l32*Njp)*C3312 + 
		   l32*Njq*C3322 + 
		   l31*Njp*C3311) + 
	  ((l22+1)*Nir+l23*Niq)*((l33+1)*Njr*C2333 + 
				 (l32*Njr+(l33+1)*Njq)*C2323 + 
				 (l31*Njr+(l33+1)*Njp)*C2313 + 
				 l32*Njq*C2322 + 
				 (l31*Njq+l32*Njp)*C2312 + 
				 l31*Njp*C2311) + 
	  (l21*Nir+l23*Nip)*((l33+1)*Njr*C1333 +
			     (l32*Njr+(l33+1)*Njq)*C1323 +
			     l32*Njq*C1322 +
			     (l31*Njr+(l33+1)*Njp)*C1313 +
			     (l31*Njq+l32*Njp)*C1312 + 
			     l31*Njp*C1311) + 
	  (l21*Niq+(l22+1)*Nip)*((l33+1)*Njr*C1233 + 
				 l32*Njq*C1222 + 
				 (l32*Njr+(l33+1)*Njq)*C1223 + 
				 (l31*Njr+(l33+1)*Njp)*C1213 + 
				 (l31*Njq+l32*Njp)*C1212 + 
				 l31*Njp*C1211) +
	  (l22+1)*Niq*((l33+1)*Njr*C2233 + 
		       (l32*Njr+(l33+1)*Njq)*C2223 + 
		       l32*Njq*C2222 + 
		       (l31*Njr+(l33+1)*Njp)*C2213 + 
		       (l31*Njq+l32*Njp)*C2212 + 
		       l31*Njp*C2211) + 
	  l21*Nip*((l33+1)*Njr*C1133 + 
		   (l32*Njr+(l33+1)*Njq)*C1123 + 
		   l32*Njq*C1122 + 
		   (l31*Njr+(l33+1)*Njp)*C1113 + 
		   (l31*Njq+l32*Njp)*C1112 + 
		   l31*Njp*C1111);

	const PylithScalar Ki2j0 =
	  (l33+1)*Nir*(l13*Njr*C3333 + 
		       (l12*Njr+l13*Njq)*C3323 + 
		       ((l11+1)*Njr+l13*Njp)*C3313 + 
		       ((l11+1)*Njq+l12*Njp)*C3312 + 
		       l12*Njq*C3322 + 
		       (l11+1)*Njp*C3311) + 
	  (l32*Nir+(l33+1)*Niq)*(l13*Njr*C2333 + 
				 (l12*Njr+l13*Njq)*C2323 + 
				 ((l11+1)*Njr+l13*Njp)*C2313 + 
				 l12*Njq*C2322 + 
				 ((l11+1)*Njq+l12*Njp)*C2312 + 
				 (l11+1)*Njp*C2311) + 
	  (l31*Nir+(l33+1)*Nip)*(l13*Njr*C1333 + 
				 (l12*Njr+l13*Njq)*C1323 + 
				 l12*Njq*C1322 + 
				 ((l11+1)*Njr+l13*Njp)*C1313 + 
				 ((l11+1)*Njq+l12*Njp)*C1312 + 
				 (l11+1)*Njp*C1311) +
	  (l31*Niq+l32*Nip)*(l13*Njr*C1233 + 
			     l12*Njq*C1222 + 
			     (l12*Njr+l13*Njq)*C1223 + 
			     ((l11+1)*Njr+l13*Njp)*C1213 + 
			     ((l11+1)*Njq+l12*Njp)*C1212 + 
			     (l11+1)*Njp*C1211) + 
	  l32*Niq*(l13*Njr*C2233 +
		   (l12*Njr+l13*Njq)*C2223 + 
		   l12*Njq*C2222 + 
		   ((l11+1)*Njr+l13*Njp)*C2213 + 
		   ((l11+1)*Njq+l12*Njp)*C2212 + 
		   (l11+1)*Njp*C2211) + 
	  l31*Nip*(l13*Njr*C1133 + 
		   (l12*Njr+l13*Njq)*C1123 + 
		   l12*Njq*C1122 + 
		   ((l11+1)*Njr+l13*Njp)*C1113 + 
		   ((l11+1)*Njq+l12*Njp)*C1112 + 
		   (l11+1)*Njp*C1111);

	const PylithScalar Ki2j1 =
	  (l33+1)*Nir*(l23*Njr*C3333 + 
		       ((l22+1)*Njr+l23*Njq)*C3323 + 
		       (l21*Njr+l23*Njp)*C3313 + 
		       (l21*Njq+(l22+1)*Njp)*C3312 + 
		       (l22+1)*Njq*C3322 + 
		       l21*Njp*C3311) + 
	  (l32*Nir+(l33+1)*Niq)*(l23*Njr*C2333 + 
				 ((l22+1)*Njr+l23*Njq)*C2323 + 
				 (l21*Njr+l23*Njp)*C2313 + 
				 (l22+1)*Njq*C2322 + 
				 (l21*Njq+(l22+1)*Njp)*C2312 + 
				 l21*Njp*C2311) + 
	  (l31*Nir+(l33+1)*Nip)*(l23*Njr*C1333 + 
				 ((l22+1)*Njr+l23*Njq)*C1323 + 
				 (l22+1)*Njq*C1322 + 
				 (l21*Njr+l23*Njp)*C1313 + 
				 (l21*Njq+(l22+1)*Njp)*C1312 + 
				 l21*Njp*C1311) + 
	  (l31*Niq+l32*Nip)*(l23*Njr*C1233 + 
			     (l22+1)*Njq*C1222 + 
			     ((l22+1)*Njr+l23*Njq)*C1223 + 
			     (l21*Njr+l23*Njp)*C1213 + 
			     (l21*Njq+(l22+1)*Njp)*C1212 + 
			     l21*Njp*C1211) + 
	  l32*Niq*(l23*Njr*C2233 + 
		   ((l22+1)*Njr+l23*Njq)*C2223 + 
		   (l22+1)*Njq*C2222 + 
		   (l21*Njr+l23*Njp)*C2213 + 
		   (l21*Njq+(l22+1)*Njp)*C2212 + 
		   l21*Njp*C2211) + 
	  l31*Nip*(l23*Njr*C1133 + 
		   ((l22+1)*Njr+l23*Njq)*C1123 + 
		   (l22+1)*Njq*C1122 + 
		   (l21*Njr+l23*Njp)*C1113 + 
		   (l21*Njq+(l22+1)*Njp)*C1112 + 
		   l21*Njp*C1111);

	const PylithScalar Ki2j2 =
	  (l33+1)*Nir*((l33+1)*Njr*C3333 + 
		       (l32*Njr+(l33+1)*Njq)*C3323 + 
		       (l31*Njr+(l33+1)*Njp)*C3313 + 
		       (l31*Njq+l32*Njp)*C3312 + 
		       l32*Njq*C3322 + 
		       l31*Njp*C3311) + 
	  (l32*Nir+(l33+1)*Niq)*((l33+1)*Njr*C2333 + 
				 (l32*Njr+(l33+1)*Njq)*C2323 + 
				 (l31*Njr+(l33+1)*Njp)*C2313 + 
				 l32*Njq*C2322 + 
				 (l31*Njq+l32*Njp)*C2312 + 
				 l31*Njp*C2311) + 
	  (l31*Nir+(l33+1)*Nip)*((l33+1)*Njr*C1333 + 
				 (l32*Njr+(l33+1)*Njq)*C1323 + 
				 l32*Njq*C1322 + 
				 (l31*Njr+(l33+1)*Njp)*C1313 + 
				 (l31*Njq+l32*Njp)*C1312 + 
				 l31*Njp*C1311) + 
	  (l31*Niq+l32*Nip)*((l33+1)*Njr*C1233 + 
			     l32*Njq*C1222 + 
			     (l32*Njr+(l33+1)*Njq)*C1223 + 
			     (l31*Njr+(l33+1)*Njp)*C1213 + 
			     (l31*Njq+l32*Njp)*C1212 + 
			     l31*Njp*C1211) + 
	  l32*Niq*((l33+1)*Njr*C2233 + 
		   (l32*Njr+(l33+1)*Njq)*C2223 + 
		   l32*Njq*C2222 + 
		   (l31*Njr+(l33+1)*Njp)*C2213 + 
		   (l31*Njq+l32*Njp)*C2212 + 
		   l31*Njp*C2211) + 
	  l31*Nip*((l33+1)*Njr*C1133 + 
		   (l32*Njr+(l33+1)*Njq)*C1123 + 
		   l32*Njq*C1122 + 
		   (l31*Njr+(l33+1)*Njp)*C1113 + 
		   (l31*Njq+l32*Njp)*C1112 + 
		   l31*Njp*C1111);

	const PylithScalar Knl = 
	  Nir*(Njr*s33+Njq*s23+Njp*s13) + 
	  Niq*(Njr*s23+Njq*s22+Njp*s12) + 
	  Nip*(Njr*s13+Njq*s12+Njp*s11);

	const int iBlock = iB * (numBasis*spaceDim);
	const int iBlock1 = (iB+1) * (numBasis*spaceDim);
	const int iBlock2 = (iB+2) * (numBasis*spaceDim);
	const int jBlock = jB;
	const int jBlock1 = jB+1;
	const int jBlock2 = jB+2;
	_cellMatrix[iBlock +jBlock ] += Ki0j0 + Knl;
	_cellMatrix[iBlock +jBlock1] += Ki0j1;
	_cellMatrix[iBlock +jBlock2] += Ki0j2;
	_cellMatrix[iBlock1+jBlock ] += Ki1j0;
	_cellMatrix[iBlock1+jBlock1] += Ki1j1 + Knl;
	_cellMatrix[iBlock1+jBlock2] += Ki1j2;
	_cellMatrix[iBlock2+jBlock ] += Ki2j0;
	_cellMatrix[iBlock2+jBlock1] += Ki2j1;
	_cellMatrix[iBlock2+jBlock2] += Ki2j2 + Knl;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(3+numBasis*(6*26+9))));
} // _elasticityJacobian3D

// ----------------------------------------------------------------------
// Calculate Green-Lagrange strain tensor at quadrature points of a 2-D cell.
void 
pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D(scalar_array* strain,
								     const scalar_array& deform,
								     const int numQuadPts)
{ // _calcTotalStrain2D
  // Green-Lagrange strain tensor = 1/2 ( X^T X - I )
  // X: deformation gradient tensor
  // I: identity matrix

  assert(strain);

  const int dim = 2;
  const int deformSize = dim*dim;
  const int strainSize = 3;
  assert(deform.size() == size_t(numQuadPts*deformSize));
  assert(strain->size() == size_t(numQuadPts*strainSize));

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
pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D(scalar_array* strain,
								     const scalar_array& deform,
								     const int numQuadPts)
{ // _calcTotalStrain3D
  // Green-Lagrange strain tensor = 1/2 ( X^T X - I )
  // X: deformation gradient tensor
  // I: identity matrix

  assert(strain);

  const int dim = 3;
  const int deformSize = dim*dim;
  const int strainSize = 6;
  assert(deform.size() == size_t(numQuadPts*dim*dim));
  assert(strain->size() == size_t(numQuadPts*strainSize));

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
pylith::feassemble::IntegratorElasticityLgDeform::_calcDeformation(scalar_array* deform,
								   const scalar_array& basisDeriv,
								   const PylithScalar* disp,
								   const int numBasis,
								   const int numQuadPts,
								   const int dim)
{ // _calcDeformation
  assert(deform);

  assert(basisDeriv.size() == size_t(numQuadPts*numBasis*dim));

  // Xij = d(u_i)/d(x_j) + delta_ij
  //
  // u_i = Nb * u^b_i
  //
  /// Xij = dNb/dxj * u^b_i + delta_ij

  const int deformSize = dim*dim;

  (*deform) = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis) {
      for (int iDim=0, indexD=0; iDim < dim; ++iDim) {
	for (int jDim=0; jDim < dim; ++jDim, ++indexD) {
	  const int iB = iBasis*dim;
	  (*deform)[iQuad*deformSize+indexD] += basisDeriv[iQ+iB+jDim] * disp[iB+iDim];
	} // for
      } // for
    } // for
    for (int iDim=0; iDim < dim; ++iDim) {
      (*deform)[iQuad*deformSize+iDim*dim+iDim] += 1.0;
    } // for
  } // for

} // _calcDeformation


// ----------------------------------------------------------------------
// Calculate 2-D Cauchy stress from 2nd Piola-Kirchoff stress.
void
pylith::feassemble::IntegratorElasticityLgDeform::_calcCauchyStress2D(scalar_array* cauchyStress,
								      const scalar_array& stress,
								      const scalar_array& deform,
								      const int numQuadPts)
{ // _calcCauchyStress2D
  PYLITH_METHOD_BEGIN;

  const int spaceDim = 2;
  const int tensorSize = 3;

  assert(cauchyStress);
  assert(stress.size() == size_t(numQuadPts*tensorSize));
  assert(cauchyStress->size() == stress.size());
  assert(deform.size() == size_t(numQuadPts*spaceDim*spaceDim));

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {

    const PylithScalar* X = &deform[iQuad*spaceDim*spaceDim];
    const PylithScalar detX = X[0]*X[3] - X[1]*X[2];
	    
    // Stress
    // i=0, j=0, index=0
    // i=1, j=1, index=1
    // i=0, j=1, index=2
    
    // Deformation gradient
    // i=0, j=0, index=0
    // i=0, j=1, index=1
    // i=1, j=0, index=2
    // i=1, j=1, index=3
    
    // S00
    (*cauchyStress)[iQuad*tensorSize+0] =
      (X[0]*X[0]*stress[0] +
       X[0]*X[1]*stress[2] +
       X[1]*X[0]*stress[2] +
       X[1]*X[1]*stress[1]) / detX;
    
    // S01
    (*cauchyStress)[iQuad*tensorSize+2] =
      (X[0]*X[2]*stress[0] +
       X[0]*X[3]*stress[2] +
       X[1]*X[2]*stress[2] +
       X[1]*X[3]*stress[1]) / detX;
    
    // S11
    (*cauchyStress)[iQuad*tensorSize+1] =
      (X[2]*X[2]*stress[0] +
       X[2]*X[3]*stress[2] +
       X[3]*X[2]*stress[2] +
       X[3]*X[3]*stress[1]) / detX;

  } // for

  PYLITH_METHOD_END;
} // _calcCauchySterss2D


// ----------------------------------------------------------------------
// Calculate 3-D Cauchy stress from 2nd Piola-Kirchoff stress.
void
pylith::feassemble::IntegratorElasticityLgDeform::_calcCauchyStress3D(scalar_array* cauchyStress,
								      const scalar_array& stress,
								      const scalar_array& deform,
								      const int numQuadPts)
{ // _calcCauchyStress3D
  PYLITH_METHOD_BEGIN;

  const int spaceDim = 3;
  const int tensorSize = 5;

  assert(cauchyStress);
  assert(stress.size() == size_t(numQuadPts*tensorSize));
  assert(cauchyStress->size() == stress.size());
  assert(deform.size() == size_t(numQuadPts*spaceDim*spaceDim));

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const PylithScalar* X = &deform[iQuad*spaceDim*spaceDim];
    const PylithScalar detX = 
        X[0] * (X[4]*X[8] - X[5]*X[7])
      - X[1] * (X[5]*X[6] - X[3]*X[8])
      + X[2] * (X[3]*X[7] - X[4]*X[6]);
    
    // Stress
    // i=0, j=0, index=0
    // i=1, j=1, index=1
    // i=2, j=2, index=2
    // i=0, j=1, index=3
    // i=1, j=2, index=4
    // i=0, j=2, index=5
    
    // Deformation gradient
    // i=0, j=0, index=0
    // i=0, j=1, index=1
    // i=0, j=2, index=2
    // i=1, j=0, index=3
    // i=1, j=1, index=4
    // i=1, j=2, index=5
    // i=2, j=0, index=6
    // i=2, j=1, index=7
    // i=2, j=2, index=8
    
    // S00
    (*cauchyStress)[iQuad*tensorSize+0] = 
      (X[0]*X[0]*stress[0] +
       X[0]*X[1]*stress[3] +
       X[0]*X[2]*stress[5] +
       X[1]*X[0]*stress[3] +
       X[1]*X[1]*stress[1] +
       X[1]*X[2]*stress[4] +
       X[2]*X[0]*stress[5] +
       X[2]*X[1]*stress[4] +
       X[2]*X[2]*stress[2]) / detX;

    // S01
    (*cauchyStress)[iQuad*tensorSize+3] = 
      (X[0]*X[3]*stress[0] +
       X[0]*X[4]*stress[3] +
       X[0]*X[5]*stress[5] +
       X[1]*X[3]*stress[3] +
       X[1]*X[4]*stress[1] +
       X[1]*X[5]*stress[4] +
       X[2]*X[3]*stress[5] +
       X[2]*X[4]*stress[4] +
       X[2]*X[5]*stress[2]) / detX;

    // S02
    (*cauchyStress)[iQuad*tensorSize+5] = 
      (X[0]*X[6]*stress[0] +
       X[0]*X[7]*stress[3] +
       X[0]*X[8]*stress[5] +
       X[1]*X[6]*stress[3] +
       X[1]*X[7]*stress[1] +
       X[1]*X[8]*stress[4] +
       X[2]*X[6]*stress[5] +
       X[2]*X[7]*stress[4] +
       X[2]*X[8]*stress[2]) / detX;

    // S11
    (*cauchyStress)[iQuad*tensorSize+1] = 
      (X[3]*X[3]*stress[0] +
       X[3]*X[4]*stress[3] +
       X[3]*X[5]*stress[5] +
       X[4]*X[3]*stress[3] +
       X[4]*X[4]*stress[1] +
       X[4]*X[5]*stress[4] +
       X[5]*X[3]*stress[5] +
       X[5]*X[4]*stress[4] +
       X[5]*X[5]*stress[2]) / detX;

    // S12
    (*cauchyStress)[iQuad*tensorSize+4] = 
      (X[3]*X[6]*stress[0] +
       X[3]*X[7]*stress[3] +
       X[3]*X[8]*stress[5] +
       X[4]*X[6]*stress[3] +
       X[4]*X[7]*stress[1] +
       X[4]*X[8]*stress[4] +
       X[5]*X[6]*stress[5] +
       X[5]*X[7]*stress[4] +
       X[5]*X[8]*stress[2]) / detX;

    // S22
    (*cauchyStress)[iQuad*tensorSize+2] = 
      (X[6]*X[6]*stress[0] +
       X[6]*X[7]*stress[3] +
       X[6]*X[8]*stress[5] +
       X[7]*X[6]*stress[3] +
       X[7]*X[7]*stress[1] +
       X[7]*X[8]*stress[4] +
       X[8]*X[6]*stress[5] +
       X[8]*X[7]*stress[4] +
       X[8]*X[8]*stress[2]) / detX;
    
  } // for

  PYLITH_METHOD_END;
} // _calcCauchyStress3D


// End of file 
