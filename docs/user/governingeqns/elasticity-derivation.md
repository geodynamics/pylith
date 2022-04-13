# Derivation of Elasticity Equation

For completeness we start our discussion of the governing equations with a derivation of the elasticity equation.
Consider domain $\Omega$ bounded by boundary $\Gamma$.
Applying a Lagrangian description of the conservation of momentum gives
%
```{math}
:label: eqn:momentum:vec
\frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\, d\Omega=\int_{\Omega}\vec{f}(\vec{x},t)\, d\ + \int_{\Gamma}\vec{\tau}(\vec{x},t)\, d\Gamma.
```
%
The traction vector field is related to the stress tensor through
%
\begin{equation}
\vec{\tau}(\vec{x},t) = \boldsymbol{\sigma}(\vec{u}) \cdot \vec{n},
\end{equation}
%
where $\vec{n}$ is the outward normal vector to $\Gamma$.
Substituting into equation {math:numref}`eqn:momentum:vec` yields
%
\begin{equation}
\frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\, d\Omega = \int_{\Omega}\vec{f}(\vec{x},t)\, d\Omega+\int_{\Gamma}\boldsymbol{\sigma}(\vec{u})\cdot\vec{n}\, d\Gamma.
\end{equation}
%
Applying the divergence theorem,
%
\begin{equation}
\int_{\Omega}\boldsymbol{\nabla}\cdot\vec{a}\: d\Omega=\int_{\Gamma}\vec{a}\cdot\vec{n}\: d\Gamma,
\end{equation}
%
to the boundary integral results in
%
\begin{equation}
\frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\, d\Omega=\int_{\Omega}\vec{f}(\vec{x},t)\, d\Omega+\int_{\Omega}\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u})\, d\Omega,
\end{equation}
%
which we can rewrite as
%
\begin{equation}
\int_{\Omega}\left(\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u})\right)\, d\Omega=\vec{0}.
\end{equation}
%
Because the domain $\Omega$ is arbitrary, the integrand must be the zero vector at every location in the domain, so that we end up with
%
\begin{gather}
\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}=\vec{0}\text{ in }\Omega,\\
\boldsymbol{\sigma}(\vec{u})\cdot\vec{n}=\vec{\tau}(\vec{x},t)\text{ on }\Gamma_{\tau}\text{,}\\
\vec{u}=\vec{u}_0(\vec{x},t)\text{ on }\Gamma_{u},\text{ and}\\
\vec{u}^{+}-\vec{u}^{-}=\vec{d}\text{ on }\Gamma_{f}.
\end{gather}
%
We specify tractions, $\vec{\tau}$, on boundary $\Gamma_{f}$, displacements, $\vec{u^{o}}$, on boundary $\Gamma_{u}$, and slip, $\vec{d}$, on fault interface $\Gamma_{f}$.
