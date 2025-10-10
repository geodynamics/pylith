// =================================================================================================
// This code is part of SpatialData, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "scalesfwd.hh"

/// C++ object for scales related to elasticity.
class pylith::scales::ElasticityScales {
    friend class TestElasticityScales; // Unit testing

public:

    // PUBLIC METHODS /////////////////////////////////////////////////////

    /** Set defaults scales for quasi-static elasticity.
     *
     * @param[inout] Scales for nondimensionalization.
     * @param[in] lengthScale Default length scale in meters.
     * @param[in] timeScale Default time scale in seconds.
     *
     */
    static
    void setQuasistaticElasticity(pylith::scales::Scales* scales,
                                  const double lengthScale=100.0e+3,
                                  const double timeScale=31557600.0e+2);

    /** Set defaults scales for dynamic elasticity.
     *
     * @param[inout] Scales for nondimensionalization.
     * @param[in] lengthScale Default length scale in meters.
     * @param[in] velocityScale Default velocity scale in m/s.
     *
     */
    static
    void setDynamicElasticity(pylith::scales::Scales* scales,
                              const double lengthScale=100.0e+3,
                              const double velocityScale=3.0e+3);

    /** Set defaults scales for quasi-static poroelasticity.
     *
     * @param[inout] Scales for nondimensionalization.
     * @param[in] lengthScale Default length scale in meters.
     * @param[in] permeability Default permeability scale in m^2.
     * @param[in] viscosity Default viscosity scale in Pa*s.
     * @param[in] rigidity Default rigidity scale in Pa.
     */
    static
    void setQuasistaticPoroelasticity(pylith::scales::Scales* scales,
                                      const double lengthScale=100.0e+3,
                                      const double permeability=1.0e-12,
                                      const double viscosity=1.0e-3,
                                      const double rigidity=25.0e+9);

    /** Set time scale for poroelasticity.
     *
     * @param[in] lengthScale Default length scale in meters.
     * @param[in] permeability Default permeability scale in m^2.
     * @param[in] viscosity Default viscosity scale in Pa*s.
     * @param[in] rigidity Default rigidity scale in Pa.
     *
     * @returns Time scale in seconds.
     */
    static
    double computePoroelasticityTimeScale(const double viscosity,
                                          const double permeability,
                                          const double length,
                                          const double rigidity);

    /** Get value to nondimensionalize stress.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Stress scale in Pa (SI units).
     */
    static
    double getStressScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize pressure.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Fluid pressure scale in Pascals (SI units).
     */
    static
    double getFluidPressureScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize strain.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Strain scale.
     */
    static
    double getStrainScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize body force.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Body force scale in Pa/m (SI units).
     */
    static
    double getBodyForceScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize density.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Density scale in kg/m^3 (SI units).
     */
    static
    double getDensityScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize velocity.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Velocity scale in m/s (SI units).
     */
    static
    double getVelocityScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize acceleration.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Acceleration scale in m/s^2 (SI units).
     */
    static
    double getAccelerationScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize fluid viscosity.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Fluid viscosity scale in Pa*s (SI units).
     */
    static
    double getViscosityScale(const pylith::scales::Scales& scales);

    /** Get value to nondimensionalize permeability.
     *
     * @param[in] Scales for nondimensionalization.
     * @returns Permeability scale in m/s^2 (SI units).
     */
    static
    double getPermeabilityScale(const pylith::scales::Scales& scales);

}; // class ElasticityScales

// End of file
