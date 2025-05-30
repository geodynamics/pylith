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

namespace pylith {
    namespace fekernels {
        class Tensor;
        class TensorOps;

        class Solution;
        class DispVel;

        class Elasticity;
        class ElasticityPlaneStrain;
        class Elasticity3D;

        class IsotropicLinearElasticity;
        class IsotropicLinearElasticityPlaneStrain;
        class IsotropicLinearElasticity3D;

        class IsotropicLinearMaxwell;
        class IsotropicLinearMaxwellPlaneStrain;
        class IsotropicLinearMaxwell3D;

        class IsotropicLinearGenMaxwell;
        class IsotropicLinearGenMaxwellPlaneStrain;
        class IsotropicLinearGenMaxwell3D;

        class IsotropicPowerLaw;
        class IsotropicPowerLawPlaneStrain;
        class IsotropicPowerLaw3D;

        class IncompressibleElasticity;
        class IsotropicLinearIncompElasticity;
        class IsotropicLinearIncompElasticityPlaneStrain;
        class IsotropicLinearIncompElasticity3D;

        class Poroelasticity;
        class PoroelasticityPlaneStrain;
        class Poroelasticity3D;

        class IsotropicLinearPoroelasticity;
        class IsotropicLinearPoroelasticityPlaneStrain;
        class IsotropicLinearPoroelasticity3D;

        class TimeDependentFn;
        class NeumannTimeDependent;
        class AbsorbingDampers;

        class FaultCohesive;
        class FaultCohesiveKin;

        class BoundaryDirections;
    } // fekernels
} // pylith

// End of file
