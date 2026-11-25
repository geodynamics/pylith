// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#pragma once

namespace pylith {
    namespace materials {
        // New stuff
        class Material;

        class Elasticity;
        class AuxiliaryFactoryElasticity;
        class DerivedFactoryElasticity;

        class RheologyElasticity;
        class IsotropicLinearElasticity;
        class IsotropicLinearMaxwell;
        class IsotropicLinearGenMaxwell;
        class IsotropicPowerLaw;
        class IsotropicLinearIncompElasticity;
        class AuxiliaryFactoryElastic;
        class AuxiliaryFactoryViscoelastic;

        class IncompressibleElasticity;

        class RheologyIncompressibleElasticity;
        class IsotropicLinearIncompElasticity;

        class Poroelasticity;
        class AuxiliaryFactoryPoroelasticity;
        class DerivedFactoryPoroelasticity;

        class RheologyPoroelasticity;
        class IsotropicLinearPoroelasticity;
        class AuxiliaryFactoryPoroelastic;

        class Query;

    } // namespace materials
} // namespace pylith


// End of file
