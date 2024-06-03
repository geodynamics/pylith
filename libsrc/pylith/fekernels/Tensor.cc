/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */

#include <portinfo>

#include "pylith/fekernels/Tensor.hh"

const pylith::fekernels::TensorOps pylith::fekernels::Tensor::ops2D =
    pylith::fekernels::TensorOps::_create2D();

const pylith::fekernels::TensorOps pylith::fekernels::Tensor::ops3D =
    pylith::fekernels::TensorOps::_create3D();

// End of file
