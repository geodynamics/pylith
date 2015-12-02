/* ====================================================================== 
 * Kernels for stress.
 *
 * \int_V \nabla \vec{\phi}_u : \tensor{\sigma} \, dV

g1_IsotropicLinearElasticityPlaneStrain

Jg3_dispdisp_IsotropicLinearElasticityPlaneStrain

 * ====================================================================== 
 */

/** f1 entry function for isotropic linear elasticity in 3-D.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [2].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param x Coordinates of point evaluation.
 * @param f1 Result [dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1)]
 */
void
pylith_fekernels_f1_IsotropicLinearElasticity3D(const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar f1[]);
					      

/** g3_uu entry function for isotropic linear elasticity in 3-D.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [2].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param utshift Coefficient for dF/ds_t term in Jacobian.
 * @param x Coordinates of point evaluation.
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_g3_uu_IsotropicLinearElasticity3D(const PylithInt dim,
						   const PylithInt numS,
						   const PylithInt numA,
						   const PylithInt sOff[],
						   const PylithInt sOff_x[],
						   const PylithScalar s[],
						   const PylithScalar s_t[],
						   const PylithScalar s_x[],
						   const PylithInt aOff[],
						   const PylithInt aOff_x[],
						   const PylithScalar a[],
						   const PylithScalar a_t[],
						   const PylithScalar a_x[],
						   const PylithReal t,
						   const PylithReal utshift,
						   const PylithScalar x[],
						   PylithScalar g3[]);
					      

/** Calculate volumetric stress for 3-D isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [1].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param x Coordinates of point evaluation.
 * @param f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(const PylithInt dim,
							      const PylithInt numS,
							      const PylithInt numA,
							      const PylithInt sOff[],
							      const PylithInt sOff_x[],
							      const PylithScalar s[],
							      const PylithScalar s_t[],
							      const PylithScalar s_x[],
							      const PylithInt aOff[],
							      const PylithInt aOff_x[],
							      const PylithScalar a[],
							      const PylithScalar a_t[],
							      const PylithScalar a_x[],
							      const PylithReal t,
							      const PylithScalar x[],
							      PylithScalar stress[]);


/** Calculate deviatoric stress for 3-D isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [1].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param x Coordinates of point evaluation.
 * @param f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(const PylithInt dim,
							      const PylithInt numS,
							      const PylithInt numA,
							      const PylithInt sOff[],
							      const PylithInt sOff_x[],
							      const PylithScalar s[],
							      const PylithScalar s_t[],
							      const PylithScalar s_x[],
							      const PylithInt aOff[],
							      const PylithInt aOff_x[],
							      const PylithScalar a[],
							      const PylithScalar a_t[],
							      const PylithScalar a_x[],
							      const PylithReal t,
							      const PylithScalar x[],
							      PylithScalar stress[]);

/** f1 entry function for 3-D incompressible isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [4].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param x Coordinates of point evaluation.
 * @param f0 Result [dim].
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_f1_IncomprUIntegral3D(const PylithInt dim,
				       const PylithInt numS,
				       const PylithInt numA,
				       const PylithInt sOff[],
				       const PylithInt sOff_x[],
				       const PylithScalar s[],
				       const PylithScalar s_t[],
				       const PylithScalar s_x[],
				       const PylithInt aOff[],
				       const PylithInt aOff_x[],
				       const PylithScalar a[],
				       const PylithScalar a_t[],
				       const PylithScalar a_x[],
				       const PylithReal t,
				       const PylithScalar x[],
				       PylithScalar f1[]);

/** Calculate deviatoric stress for 3-D incompressible isotropic linear
 *  elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [3].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param x Coordinates of point evaluation.
 * @param stress Result [dim*dim].
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticityIncompr3D(
						const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar stress[]);


/** g3_uu entry function for 3-D incompressible isotropic linear elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [2].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param utshift Coefficient for dF/ds_t term in Jacobian.
 * @param x Coordinates of point evaluation.
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_g3_uu_IncomprIsotropicLinearElasticity3D(
						const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithReal utshift,
						const PylithScalar x[],
						PylithScalar g3[]);

/** g2_uv entry function for 3-D incompressible isotropic linear elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [2].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param utshift Coefficient for dF/ds_t term in Jacobian.
 * @param x Coordinates of point evaluation.
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_g2_uv_IncomprIsotropicLinearElasticity3D(
						const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithReal utshift,
						const PylithScalar x[],
						PylithScalar g2[]);

#endif /* pylith_fekernels_elasticity_h */
