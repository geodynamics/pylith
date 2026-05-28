# Notes

At run time, a user specifies:

1. Governing equation (elasticity, incompressible elasticity, poroelasticity, ...)
   1. Elasticity
      1. Body force: yes/no (affects auxiliary field layout)
      2. Gravitational acceleration: yes/no (affects auxiliary field layout)
      3. Bulk rheology (isotropic linear, isotropic maxwell, isotropic powerlaw)
         1. IsotropicLinear
            1. Reference stress: yes/not (affects auxiliary field layout)
2. Inertia: yes/no (affects kernels and solution field layout)
3. Fault(s): yes/no (affects solution field layout); kernels added for each fault
4. Strain model: infinitesimal or finite (depends on dim, solution layout)

## PETSc interface

`∫Ω 𝜓^𝑢_trial ⋅ 𝑓0⁡(𝑡,\dot{𝑠},𝑠) + ∇𝜓^𝑢_trial : 𝒇1⁡(𝑡,\dot{𝑠},𝑠) ⁢𝑑⁢Ω = ∫Ω 𝜓^𝑢_trial ⋅ 𝑔0⁡(𝑡,𝑠) + ∇𝜓^𝑢_trial : 𝒈1⁡(𝑡,𝑠) ⁢𝑑⁢Ω`

PETSc has cell geometry, basis functions, and quadrature information.
PyLith supplies kernels for integration, e.g., f0, f1.

### f0 kernel for domain

```c++
void f0(const pylith::integer cellDim,
        const pylith::integer numS,
        const pylith::integer numA,
        const pylith::integer sOff[],
        const pylith::integer sOff_x[],
        const pylith::scalar s[],
        const pylith::scalar s_t[],
        const pylith::scalar s_x[],
        const pylith::integer aOff[],
        const pylith::integer aOff_x[],
        const pylith::scalar a[],
        const pylith::scalar a_t[],
        const pylith::scalar a_x[],
        const pylith::real t,
        const pylith::real x[],
        const pylith::integer numConstants,
        const pylith::scalar constants[],
        pylith::scalar f0[]);
```

## Momentum equation

Integration of momentum PDE requiries:

1. Body force term `f0 = f_i + ρ g_i`
   1. Spatial dimension
   2. Layout of auxiliary field
2. Divergence of stress term `f1 = -σ`
   1. Spatial dimension
   2. Strain model
   3. Stress model
   4. Layout of solution field
   5. Layout of auxiliary field
3. Jacobian for divergence of stress term `Jf3uu = j(f,g,df,dg) = C(f,df,g,dg)`
   1. Spatial dimension
   2. Strain model (maybe)
   3. Stress model
   4. Layout of solution field (maybe)
   5. Layout of auxiliary field

## Status

1. Implemented a static map to hold kernels to facilitate retrieval.
   1. Each rheology creates all kernels for a given PDE+rheology.
   2. Retrieval for given rheology is based on bitwise flags.

## Questions

1. Any obvious issues?
2. Strategies for organizing header files
   1. When is it reasonable to put multiple classes/structs/enums in one header file?
      1. Use `class` with `public` (new GPU compilers can handle C++).
      2. Avoid multiple classes in files.
      3. Move all inline functions to `.icc`
      4. Gather all Python exceptions in exceptions.py and would do same thing in C++. Strong exception hierarchy.
      5. Enums organize as necessary. Include in appropriate files.
      6. What about `fekernels` with just header files?
   2. Differentiating between the internal and public interface (implementation and conveying to users)
      1. Should public interface be `#include "pylith.hh"`?
         1. How to bridge gap between`pylith.hh` and header files in all of the directories?
      2. `lib/pyre/grid.h`
         1. What are the roles of `public.h`, `externals.h`, `forward.h`, `api.h`?
            1. `grid.h` will give you access to public interface.
               1. Long description of what package does.
               2. Pull in `public.h` and publish public headers.
               3. grid is a "bad" example because it is heavily templated.
               4. gcc bug with `#pragma once`: foo.h with same contents are not differentiated. Can use comment to differentiate.
            2. forward.h - sets namespace and classes.
            3. externals.h - Intent was to put all includes from external packages. Precompiled headers. Not used much anymore because link step is the bottleneck.
            4. api.h - templates with default arguments. More meaningful names.
            5. Packages will use public interface from sister subpackage.
      3. Helpers in own files and not exposed to user.
      4. concepts
         1. If kernel flags were types
         2. requires `class` class expansion for constraint is . If false at compile time, no expansion.
         3. TODO: Look up C++ concepts.
