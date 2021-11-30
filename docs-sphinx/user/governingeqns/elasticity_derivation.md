# Derivation of Elasticity Equation

For completeness we start our discussing of the governing equations with a derivation of the elasticity equation.
Consider domain {math}`\Omega` bounded by boundary {math}`\Gamma`.
Applying a Lagrangian description of the conservation of momentum gives
```{math} \frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\ d\Omega=\int_{\Omega}\vec{f}(\vec{x},t)\ d\ + \int_{\Gamma}\vec{\tau}(\vec{x},t)\ d\Gamma.
---
label: eqn:momentum:vec
---
```
The traction vector field is related to the stress tensor through
```{math} \vec{\tau}(\vec{x},t) = \boldsymbol{\sigma}(\vec{u}) \cdot \vec{n},
```
where {math}`\vec{n}` is the outward normal vector to {math}`\Gamma`.
Substituting into equation {math:numref}`eqn:momentum:vec` yields
```{math} \frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\ d\Omega=\int_{\Omega}\vec{f}(\vec{x},t)\ d\Omega+\int_{\Gamma}\boldsymbol{\sigma}(\vec{u})\cdot\vec{n}\ d\Gamma.
```
Applying the divergence theorem,
```{math} \int_{\Omega}\boldsymbol{\nabla}\cdot\vec{a}\ d\Omega=\int_{\Gamma}\vec{a}\cdot\vec{n}\ d\Gamma,
```
to the boundary integral results in
```{math} \frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\ d\Omega=\int_{\Omega}\vec{f}(\vec{x},t)\ d\Omega+\int_{\Omega}\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u})\ d\Omega,
```
which we can rewrite as
```{math} \int_{\Omega}\left(\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u})\right)\  d\Omega=\vec{0}.
```
Because the domain {math}`\Omega` is arbitrary, the integrand must be the zero vector at every location in the domain, so that we end up with
```{math} \rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}=\vec{0}\text{ in }\Omega,
```
```{math}
\boldsymbol{\sigma}(\vec{u})\cdot\vec{n}=\vec{\tau}(\vec{x},t)\text{ on }\Gamma_{\tau}\text{,}
```
```{math}
\vec{u}=\vec{u}_0(\vec{x},t)\text{ on }\Gamma_{u},\text{ and}
```
```{math}
\vec{u}^{+}-\vec{u}^{-}=\vec{d}\text{ on }\Gamma_{f}.
```
We specify tractions, {math}`\vec{\tau}`, on boundary {math}`\Gamma_{f}`, displacements,
{math}`\vec{u^{o}}`, on boundary {math}`\Gamma_{u}`, and slip, {math}`\vec{d}`, on fault
interface {math}`\Gamma_{f}`.
