(sec-user-physics-absorbing-dampers)=
# Absorbing Boundary Conditions

You can use the `AbsorbingDampers` boundary condition to prevent seismic waves reflecting off of a boundary.
Normally incident dilatational and shear waves are perfectly absorbed.
Waves incident at other angles are only partially absorbed.
This boundary condition is simpler than a perfectly matched layer (PML) boundary condition but does not perform quite as well, especially for surface waves.
If the waves arriving at the absorbing boundary are relatively small in amplitude compared to the amplitudes of primary interest, this boundary condition gives reasonable results.

The auxiliary field spatial database contains the bulk rheology properties for an isotrpoic, linear elastic material (density, Vs (S-wave speed), and Vp (P-wave speed).
You can simply use the same spatial database that was used to specify the elastic properties of the material.

:::{seealso}
[`AbsorbingDampers` Component](../components/bc/AbsorbingDampers.md) for the Pyre properties and facilities and configuration examples.
:::

## Finite-Element Implementation of Absorbing Boundary

:::{admonition} TODO
:class: error

Move this to the multiphysics implementation section.

Consider a plane wave propagating at a velocity $c$.
We can write the displacement field as
%
\begin{equation}
\vec{u}(\vec{x},t)=\vec{u^{t}}(t-\frac{\vec{x}}{c}),
\end{equation}
%
where $\vec{x}$ is position, $t$ is time, and $\vec{u^{t}}$ is the shape of the propagating wave.
For an absorbing boundary we want the traction on the boundary to be equal to the traction associated with the wave propagating out of the domain.
Starting with the expression for the traction on a boundary, $T_{i}=\sigma_{ij}n_{j},$ and using the local coordinate system for the boundary $s_{h}s_{v}n,$ where $\vec{n}$ is the direction normal to the boundary, $\vec{s}_h$ is the horizontal direction tangent to the boundary, and $\vec{s}_v$ is the vertical direction tangent to the boundary, the tractions on the boundary are
%
\begin{gather}
T_{s_{h}}=\sigma_{s_{h}n}\\
T_{s_{v}}=\sigma_{s_{v}n}\\
T_{n}=\sigma_{nn}.
\end{gather}

In the case of a horizontal boundary, we can define an auxiliary direction in order to assign unique tangential directions.
For a linear elastic isotropic material, $\sigma_{ij}=\lambda\epsilon_{kk}\delta_{ij}+2\mu\epsilon_{ij},$ and we can write the tractions as
%
\begin{gather}
T_{s_{h}}=2\mu\epsilon_{s_{h}n}\\
T_{s_{v}}=2\epsilon_{s_{v}n}\\
T_{n}=(\lambda+2\mu)\epsilon_{nn}+\lambda(\epsilon_{s_{h}s_{h}}+\epsilon_{s_{v}s_{v}}).
\end{gather}
%
For infinitesimal strains, $\epsilon_{ij}=\frac{1}{2}(u_{i,j}+u_{j,i})$ and we have
%
\begin{gather}
\epsilon_{s_{h}n}=\frac{1}{2}(u_{s_{h},n}+u_{n,s_{h}})\\
\epsilon_{s_{v}n}=\frac{1}{2}(u_{s_{v},n}+u_{n,s_{v}})\\
\epsilon_{nn}=u_{n,n}.
\end{gather}
%
For our propagating plane wave, we recognize that
%
\begin{equation}
\frac{\partial\vec{u^{t}}(t-\frac{\vec{x}}{c})}{\partial x_{i}}=-\frac{1}{c}\frac{\partial\vec{u^{t}}(t-\frac{\vec{x}}{c})}{\partial t},
\end{equation}
%
so that our expressions for the tractions become
%
\begin{gather}
T_{s_{h}}=-\frac{\mu}{c}\left(\frac{\partial u_{s_{h}}^{t}(t-\frac{\vec{x}}{c})}{\partial t}+\frac{\partial u_{n}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\right),\\
T_{s_{v}}=-\frac{\mu}{c}\left(\frac{\partial u_{s_{v}}^{t}(t-\frac{\vec{x}}{c})}{\partial t}+\frac{\partial u_{n}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\right).
\end{gather}
%
For the normal traction, consider a dilatational wave propagating normal to the boundary at speed $v_p$; in this case $u_{s_{h}}=u_{s_{v}}=0$ and $c=v_{p}$.
For the shear tractions, consider a shear wave propagating normal to the boundary at speed $v_s$; we can decompose this into one case where $u_{n}=u_{s_{v}}=0$ and another case where $u_{n}=u_{s_{h}}=0$, with $c=v_{s}$ in both cases.
We also recognize that $\mu=\rho v_{s}^{2}$ and $\lambda+2\mu=\rho v_{p}^{2}$.
This leads to the following expressions for the tractions:
%
\begin{gather}
T_{s_{h}}=-\rho v_{s}\frac{\partial u_{s_{h}}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\\
T_{s_{v}}=-\rho v_{s}\frac{\partial u_{v}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\\
T_{n}=-\rho v_{p}\frac{\partial u_{n}^{t}(t-\frac{\vec{x}}{c})}{\partial t}
\end{gather}
%
We write the weak form of the boundary condition as
%
\begin{equation}
\int_{S_{T}}T_{i}\phi_{i}\, dS=\int_{S_{T}}-\rho c_{i}\frac{\partial u_{i}}{\partial t}\phi_{i}\, dS,
\end{equation}
%
where $c_{i}$ equals $v_{p}$ for the normal traction and $v_{s}$ for the shear tractions, and $\phi_{i}$ is our weighting function.
We express the trial solution and weighting function as linear combinations of basis functions,
%
\begin{gather}
u_{i}=\sum_{m}a_{i}^{m}N^{m},\\
\phi_{i}=\sum_{n}c_{i}^{n}N^{n}.
\end{gather}
%
Substituting into our integral over the absorbing boundaries yields
%
\begin{equation}
\int_{S_{T}}T_{i}\phi_{i}\, dS=\int_{S_{T}}-\rho c_{i}\sum_{m}\dot{a}_{i}^{m}N^{m}\sum_{n}c_{i}^{n}N^{n}\, dS.
\end{equation}
%
In the derivation of the governing equations, we recognized that the weighting function is arbitrary, so we form the residual by setting the terms associated with the coefficients $c_{i}^{n}$ to zero,
%
\begin{equation}
r_{i}^{n}=\sum_{\text{tract cells}}\sum_{\text{quad pts}}-\rho(x_{q})c_{i}(x_{q})\sum_{m}\dot{a}_{i}^{m}N^{m}(x_{q})N^{n}(x_{q})w_{q}|J_{cell}(x_{q})|,
\end{equation}
%
where $x_{q}$ are the coordinates of the quadrature points, $w_{q}$ are the weights of the quadrature points, and $|J_{cell}(x_{q})|$ is the determinant of the Jacobian matrix evaluated at the quadrature points associated with mapping the reference cell to the actual cell.

The appearance of velocity in the expression for the residual means that the absorbing dampers also contribute to the system Jacobian matrix.
Using the central difference method, the velocity is written in terms of the displacements,
%
\begin{equation}
\dot{u}_{i}(t)=\frac{1}{2\Delta t}(u_{i}(t+\Delta t)-u_{i}(t-\Delta t)).
\end{equation}
%
Expressing the displacement at time $t+\Delta t$ in terms of the displacement at time $t$ ($u_{i}(t)$) and the increment in the displacement at time $t$ ($du_{i}(t)$) leads to
%
\begin{equation}
\dot{u}_{i}(t)=\frac{1}{2\Delta t}(du_{i}(t)+u_{i}(t)-u_{i}(t-\Delta t))
\end{equation}
%
The terms contributing to the system Jacobian are associated with the increment in the displacement at time $t$.
Substituting into the governing equations and isolating the term associated with the increment in the displacement at time t yields
%
\begin{equation}
A_{ij}^{nm}=\sum_{\text{tract cells}}\sum_{\text{quad pts}}\delta_{ij}\frac{1}{2\Delta t}\rho(x_{q})v_{i}(x_{q})N^{m}(x_{q})N^{n}(x_{q})w_{q}|J_{cells}(x_{q})|,
\end{equation}
%
where $A_{ij}^{mn}$ is an $nd$ by $md$ matrix ($d$ is the dimension of the vector space), $m$ and $n$ refer to the basis functions and $i$ and $j$ are vector space components.
:::
