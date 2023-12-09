// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/fekernels/fekernelsfwd.hh
 *
 * @brief Forward declarations for PyLith fekernels objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_fekernels_fekernelsfwd_hh)
#define pylith_fekernels_fekernelsfwd_hh

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

#endif // pylith_fekernels_fekernelsfwd_hh

// End of file
