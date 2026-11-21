// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#pragma once

#include "petsc.h"

typedef PetscScalar PylithScalar;
typedef PetscReal PylithReal;
typedef PetscInt PylithInt;

typedef PetscErrorCode (*PetscUserFieldFunc)(PetscInt,
                                             PetscReal,
                                             const PetscReal x[],
                                             PetscInt,
                                             PetscScalar* u,
                                             void* ctx);


// End of file
