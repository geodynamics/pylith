/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/Viscoelasticity.hh
 *
 * Generic viscoelastic functions and kernels.
 *
 **Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

#if !defined(pylith_fekernels_viscoelasticity_hh)
#define pylith_fekernels_viscoelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::Viscoelasticity {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Viscous strain coefficient function for Maxwell viscoelastic materials.
     *
     * @param[in] dt Time step size.
     * @param[in] maxwellTime Relaxation time for material.
     *
     * @returns Viscous strain coefficient.
     */
    static
    PylithScalar maxwellViscousStrainCoeff(const PylithScalar dt,
                                           const PylithScalar maxwellTime);

	/** 2D scalar inner product of two tensors represented as vectors.
	 *
	 *  @param[in] tensor1 First tensor (tens1_xx, tens1_yy, tens1_xy)
	 *  @param[in] tensor2 Second tensor (tens2_xx, tens2_yy, tens2_xy)
	 *
	 * @returns scalar inner product of the two tensors.
	 */
	static
	PylithScalar scalarProduct2D(const PylithScalar* tensor1,
								 const PylithScalar* tensor2);

	/** 2DPS scalar inner product of two tensors represented as vectors.
	 * In this case all 3 normal components are present, but only one shear component.
	 *
	 *  @param[in] tensor1 First tensor (tens1_xx, tens1_yy, tens1_zz, tens1_xy)
	 *  @param[in] tensor2 Second tensor (tens2_xx, tens2_yy, tens2_zz, tens2_xy)
	 *
	 * @returns scalar inner product of the two tensors.
	 */
	static
	PylithScalar scalarProduct2DPS(const PylithScalar* tensor1,
								   const PylithScalar* tensor2);

	/** 3D scalar inner product of two tensors represented as vectors.
	 *
	 *  @param[in] tensor1 First tensor (tens1_xx, tens1_yy, tens1_zz, tens1_xy, tens1_yz, tens1_xz)
	 *  @param[in] tensor2 Second tensor (tens2_xx, tens2_yy, tens2_zz, tens2_xy, tens2_yz, tens2_xz)
	 *
	 * @returns scalar inner product of the two tensors.
	 */
	static
	PylithScalar scalarProduct3D(const PylithScalar* tensor1,
								 const PylithScalar* tensor2);

}; // Viscoelasticity

#endif // pylith_fekernels_viscoelasticity_hh

// End of file
