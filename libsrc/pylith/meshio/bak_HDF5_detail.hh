// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "hdf5.h"

#if H5_VERSION_GE(1,12,0)
#define PYLITH_HDF5_USE_API_112
#endif

#if H5_VERSION_GE(1,8,0)
#define PYLITH_HDF5_USE_API_18
#endif

namespace pylith {
    namespace meshio {
        namespace _HDF5 {
            /// RAII wrapper for HDF5 object handles.
            class H5Handle {
                using closer_t = herr_t (*)(hid_t);
public:

                hid_t id{-1};
                closer_t closer{nullptr};

                H5Handle() = default;
                H5Handle(hid_t h,
                         closer_t c) : id(h), closer(c) {}


                ~H5Handle() {
                    if ((id >= 0) && closer) { closer(id); }
                }

                H5Handle(const H5Handle&) = delete;
                H5Handle& operator=(const H5Handle&) = delete;

                H5Handle(H5Handle&& o) noexcept : id(o.id), closer(o.closer) {
                    o.id = -1;
                }

                H5Handle& operator=(H5Handle&& o) noexcept {
                    if (this != &o) {
                        if ((id >= 0) && closer) { closer(id); }
                        id = o.id;
                        closer = o.closer;
                        o.id = -1;
                    }
                    return *this;
                }

            };

        } // _HDF5
    } // meshio
} // pylith
