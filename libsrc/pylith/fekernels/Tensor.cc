/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2022 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/Tensor.hh"

const pylith::fekernels::TensorOps pylith::fekernels::Tensor::ops2D =
    pylith::fekernels::TensorOps::_create2D();

const pylith::fekernels::TensorOps pylith::fekernels::Tensor::ops3D =
    pylith::fekernels::TensorOps::_create3D();

// End of file
