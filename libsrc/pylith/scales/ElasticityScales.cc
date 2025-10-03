// =================================================================================================
// This code is part of SpatialData, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "ElasticityScales.hh" // implementation of class methods

#include "Scales.hh" // USES Scales

// ------------------------------------------------------------------------------------------------
// Set defaults scales for quasi-static elasticity.
void
pylith::scales::ElasticityScales::setQuasistaticElasticity(pylith::scales::Scales* scales,
                                                               const double lengthScale,
                                                               const double timeScale) {
    const double length = lengthScale;
    const double displacement = 1.0;
    const double rigidity = 2.5e+10;

    scales->setLengthScale(length);
    scales->setDisplacementScale(displacement);
    scales->setRigidityScale(rigidity);
    scales->setTimeScale(timeScale);
}


// ------------------------------------------------------------------------------------------------
// Set defaults scales for dynamic elasticity.
void
pylith::scales::ElasticityScales::setDynamicElasticity(pylith::scales::Scales* scales,
                                                           const double lengthScale,
                                                           const double velocityScale) {
    const double displacement = 1.0;
    const double rigidity = 2.25e+10;
    const double time = lengthScale / velocityScale;

    scales->setLengthScale(lengthScale);
    scales->setDisplacementScale(displacement);
    scales->setRigidityScale(rigidity);
    scales->setTimeScale(time);
}


// ------------------------------------------------------------------------------------------------
// Set defaults scales for quasi-static poroelasticity.
void
pylith::scales::ElasticityScales::setQuasistaticPoroelasticity(pylith::scales::Scales* scales,
                                                                   const double lengthScale,
                                                                   const double permeability,
                                                                   const double viscosity,
                                                                   const double rigidity) {
    const double length = lengthScale;
    const double displacement = 1.0;
    const double time = computePoroelasticityTimeScale(viscosity, permeability, length, rigidity);

    scales->setLengthScale(length);
    scales->setDisplacementScale(displacement);
    scales->setRigidityScale(rigidity);
    scales->setTimeScale(time);
}


// ------------------------------------------------------------------------------------------------
// Compute time scale for poroelasticity.
double
pylith::scales::ElasticityScales::computePoroelasticityTimeScale(const double viscosity,
                                                                     const double permeability,
                                                                     const double length,
                                                                     const double rigidity) {
    return (viscosity / permeability) * (length * length / rigidity);
} // computePoroelasticityTimeScale


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize stress.
double
pylith::scales::ElasticityScales::getStressScale(const pylith::scales::Scales& scales) {
    const double rigidity = scales.getRigidityScale();
    const double displacement = scales.getDisplacementScale();
    const double length = scales.getLengthScale();
    return rigidity * displacement / length;
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize pressure.
double
pylith::scales::ElasticityScales::getFluidPressureScale(const pylith::scales::Scales& scales) {
    return getStressScale(scales);
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize strain.
double
pylith::scales::ElasticityScales::getStrainScale(const pylith::scales::Scales& scales) {
    const double displacement = scales.getDisplacementScale();
    const double length = scales.getLengthScale();
    return displacement / length;
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize body force.
double
pylith::scales::ElasticityScales::getBodyForceScale(const pylith::scales::Scales& scales) {
    const double stress = getStressScale(scales);
    const double length = scales.getLengthScale();
    return stress / length;
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize density.
double
pylith::scales::ElasticityScales::getDensityScale(const pylith::scales::Scales& scales) {
    const double rigidity = scales.getRigidityScale();
    const double length = scales.getLengthScale();
    const double time = scales.getTimeScale();
    return (rigidity * time * time) / (length * length);
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize velocity.
double
pylith::scales::ElasticityScales::getVelocityScale(const pylith::scales::Scales& scales) {
    const double time = scales.getTimeScale();
    const double displacement = scales.getDisplacementScale();
    return displacement / time;
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize acceleration.
double
pylith::scales::ElasticityScales::getAccelerationScale(const pylith::scales::Scales& scales) {
    const double time = scales.getTimeScale();
    const double displacement = scales.getDisplacementScale();
    return displacement / (time * time);
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize viscosity.
double
pylith::scales::ElasticityScales::getViscosityScale(const pylith::scales::Scales& scales) {
    const double rigidityScale = scales.getRigidityScale();
    const double time = scales.getTimeScale();
    return rigidityScale * time;
}


// ------------------------------------------------------------------------------------------------
// Get value to nondimensionalize permability.
double
pylith::scales::ElasticityScales::getPermeabilityScale(const pylith::scales::Scales& scales) {
    const double length = scales.getLengthScale();
    return length * length;
}


// End of file
