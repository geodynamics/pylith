// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/faults/KinSrc.hh"

// KinSrcBrune ------------------------------------------------------
/** @brief Slip function time history corresponding to the integral of Brune's (1970) far-field time function.
 *
 * slip = 0 for t < t0.
 * slip = final_slip * (1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau))
 * tau = 0.21081916 * rise_time
 *
 * slip_rate = 0 for t < t0.
 * slip_rate = final_slip * (1.0/tau)*(exp(-(t-t0)/tau) * (2.0 - (t-t0)/tau))
 */
class pylith::faults::KinSrcBrune : public KinSrc {
    friend class TestKinSrcBrune; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcBrune(void);

    /// Destructor.
    ~KinSrcBrune(void);

    /** Slip time function kernel.
     *
     * The "solution" field s is ignored.
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
     * @param[out] slip [dim].
     */
    static
    void slipFn(const PylithInt dim,
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
                const PylithInt numConstants,
                const PylithScalar constants[],
                PylithScalar slip[]);

    /** Slip rate time function kernel.
     *
     * The "solution" field s is ignored.
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
     * @param[out] slipRate [dim].
     */
    static
    void slipRateFn(const PylithInt dim,
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
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    PylithScalar slipRate[]);

    /** Slip acceleration time function kernel.
     *
     * The "solution" field s is ignored.
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
     * @param[out] slipAcc [dim].
     */
    static
    void slipAccFn(const PylithInt dim,
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
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar slipAcc[]);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * @param[in] scales Scales for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void _auxiliaryFieldSetup(const spatialdata::units::Scales& scales,
                              const spatialdata::geocoords::CoordSys* cs);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    KinSrcBrune(const KinSrcBrune&); ///< Not implemented
    const KinSrcBrune& operator=(const KinSrcBrune&); ///< Not implemented

}; // class KinSrcBrune

// End of file
