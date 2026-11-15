// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

namespace pylith {
    namespace bc {
        class BoundaryCondition;
        class DiagnosticFieldFactory;
        class TimeDependentAuxiliaryFactory;

        class ConstraintBoundary;
        class Dirichlet;
        class DirichletTimeDependent;
        class DirichletUserFn;

        class NeumannTimeDependent;
        class NeumannUserFn;
        class AbsorbingDampers;

        class AbsorbingDampersAuxiliaryFactory;

    } // namespace bc
} // namespace pylith

// End of file
