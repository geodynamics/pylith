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

#include <stddef.h> // USES size_t

/// C++ object for managing parameters for nondimensionalization.
class pylith::scales::Scales {
    friend class TestScales; // Unit testing

public:

    // PUBLIC METHODS /////////////////////////////////////////////////////

    /// Default constructor
    Scales(void);

    /// Default destructor
    ~Scales(void);

    /** Copy constructor.
     *
     * @param dim Object to copy.
     */
    Scales(const Scales& dim);

    /** Assignment operator.
     *
     * @param dim Object to copy.
     * @returns Copy of this.
     */
    const Scales& operator = (const Scales& dim);

    /** Set value to nondimensionalize position.
     *
     * @param value Length scale in meters (SI units).
     */
    void setLengthScale(const double value);

    /** Get value to nondimensionalize position.
     *
     * @returns Length scale in meters (SI units).
     */
    double getLengthScale(void) const;

    /** Set value to nondimensionalize displacement.
     *
     * @param value Displacement scale in meters (SI units).
     */
    void setDisplacementScale(const double value);

    /** Get value to nondimensionalize displacement.
     *
     * @returns Displacement scale in meters (SI units).
     */
    double getDisplacementScale(void) const;

    /** Set value to nondimensionalize rigidity (e.g., elastic moduli).
     *
     * @param value Rigidity scale in Pascals (SI units).
     */
    void setRigidityScale(const double value);

    /** Get value to nondimensionalize rigidity (e.g., elastic moduli).
     *
     * @returns Rigidity scale in Pascals (SI units).
     */
    double getRigidityScale(void) const;

    /** Set value to nondimensionalize time.
     *
     * @param value Time scale in seconds (SI units).
     */
    void setTimeScale(const double value);

    /** Get value to nondimensionalize time.
     *
     * @returns Time scale in seconds (SI units).
     */
    double getTimeScale(void) const;

    /** Set value to nondimensionalize temperature.
     *
     * @param value Temperature scale in Kelvin (SI units).
     */
    void setTemperatureScale(const double value);

    /** Get value to nondimensionalize temperature.
     *
     * @returns Temperature scale in Kelvin (SI units).
     */
    double getTemperatureScale(void) const;

    /** Make value dimensionless.
     *
     * @param value Value with dimensions in SI units.
     * @param scale Scale used to nondimensionalize value.
     * @returns Dimensionless value.
     */
    double nondimensionalize(const double value,
                             const double scale) const;

    /** Make value dimensionless.
     *
     * @param value Dimensionless value.
     * @param value Value with dimensions in SI units.
     * @returns Scale used to nondimensionalize value.
     */
    double dimensionalize(const double value,
                          const double scale) const;

    /** Make values dimensionless.
     *
     * @param values Array of values with dimensions in SI units.
     * @param nvalues Number of values.
     * @param scale Scale used to nondimensionalize value.
     */
    void nondimensionalize(double* const values,
                           const size_t nvalues,
                           const double scale) const;

    /** Make values dimensionless.
     *
     * @param values Array of values with dimensions in SI units.
     * @param nvalues Number of values.
     * @param scale Scale used to nondimensionalize value.
     */
    void nondimensionalize(float* const values,
                           const size_t nvalues,
                           const double scale) const;

    /** Make value dimensionless.
     *
     * @param values Array of dimensionless values.
     * @param nvalues Number of values.
     * @param scale Scale used to nondimensionalize value.
     */
    void dimensionalize(double* const values,
                        const size_t nvalues,
                        const double scale) const;

    /** Make value dimensionless.
     *
     * @param values Array of dimensionless values.
     * @param nvalues Number of values.
     * @param scale Scale used to nondimensionalize value.
     */
    void dimensionalize(float* const values,
                        const size_t nvalues,
                        const double scale) const;

private:

    // PRIVATE MEMBERS ////////////////////////////////////////////////////

    double _length; ///< Length scale
    double _displacement; ///< Displacement scale
    double _rigidity; ///< Rigidity scale
    double _time; ///< Time scale
    double _temperature; ///< Temperature scale

}; // class Scales

#include "Scales.icc" // inline methods

// End of file
