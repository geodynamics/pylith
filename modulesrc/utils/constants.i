// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


/**
 * @file modulesrc/utils/constants.i
 *
 * @brief PyLith constants.
 */

// ----------------------------------------------------------------------
// pylith::g_acc
%inline %{
  double
  g_acc(void) {
    return pylith::g_acc;
  }
%}


// ----------------------------------------------------------------------
// pylith::max_double
%inline %{
  double
  max_double(void) {
    return pylith::max_double;
  }
%}


// ----------------------------------------------------------------------
// pylith::max_float
%inline %{
  float
  max_float(void) {
    return pylith::max_float;
  }
%}


// End of file
