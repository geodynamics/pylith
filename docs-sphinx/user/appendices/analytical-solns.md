(sec-user-appendices-analytical-solns)=
# Analytical Solutions

(sec-user-appendices-airy-stress-function)=
## Airy Stress Functions

Airy stress functions provide a simple approach for solving 2D elastoplastic problems with uniform isotropic linearly elastic material properties.
They can be useful in creating tests using the Method of Manufactured Solutions.

We start with the equilibrium equation for static elasticity in Cartesian coordinates
%
```{math}
:label: eqn:stress:fn:equilibrium
\begin{gather}
\frac{\partial\sigma_{xx}}{\partial x} + \frac{\partial\sigma_{xy}}{\partial y} + f_x = 0 \\
\frac{\partial\sigma_{yy}}{\partial y} + \frac{\partial\sigma_{xy}}{\partial x} + f_y = 0,
\end{gather}
```
%
where $f_x$ and $f_y$ are the body force components in the $x$ and $y$ directions, respectively.
We assume the body forces can be derived from a potential $\psi$
%
\begin{align}
f_x &= -\frac{\partial \psi}{\partial x}, \\
f_y &= -\frac{\partial \psi}{\partial y}.
\end{align}
We choose an Airy stress function $\phi(x,y)$ to trivially satisfy the equilibrium equations so that
%
```{math}
:label: eqn:stress:function:components
\begin{align}
\sigma_{xx} &= \frac{\partial^{2}\phi}{\partial y^{2}} + \psi,\\
\sigma_{yy} &= \frac{\partial^{2}\phi}{\partial x^{2}} + \psi,\\
\sigma_{xy} &= -\frac{\partial^{2}\phi}{\partial x\partial y}.
\end{align}
```

We must also satisfy the compatibility equations.
For plane strain, we have
%
\begin{equation}
\left(\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}\right)\left(\sigma_{xx} + \sigma_{yy}\right) - \frac{1}{1-\nu}\left(\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2}\right) = 0.
\end{equation}
%
Substituting in the expressions for the stress components in terms of the Airy stress function leads to
%
```{math}
:label: eqn:stress:function:constraint
\frac{\partial^4 \phi}{\partial x^4} + \frac{\partial^4 \phi}{\partial x^2 \partial y^2} + \frac{\partial^4 \phi}{\partial y^4} + \frac{1-2\nu}{1-\nu}\left(\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2}\right)= 0.
```

### Example

We select a second order polynomial for the Airy stress function with the form
%
\begin{equation}
\phi = \frac{1}{2} a x^2 + b x y + \frac{1}{2} c y^2.
\end{equation}
%
and no body forces ($\psi=0$).
By inspection we see that this equation trivially satisfies the equilibrium equation {math:numref}`eqn:stress:function:constraint`.
Using equation {math:numref}`eqn:stress:function:components`, we have
\begin{align}
\sigma_{xx} &= c, \\
\sigma_{yy} &= a, \\
\sigma_{xy} &= -b.
\end{align}
%
Our stress function corresponds to a uniform stress field.

Let us now consider axial extension of a rectangular block with roller boundary conditions on two sides as shown in {numref}`fig:stress:function:axial:extension`.
We have
\begin{gather}
\tau_x = \tau_0 \text{ on } x=x_1,\\
u_x = 0 \text{ on } x=x_0,\\
u_y = 0 \text{ on } y=y_0.
\end{gather}
The other boundaries are free surfaces.
Because $\sigma_{yy} = \sigma_{xy} = 0$ and $\sigma_{xx} = \tau_0$, we have $a = b = 0$ and $c = \tau_0$.
For plane strain the out of plane stress is given by
\begin{equation}
\sigma_{zz} = \frac{\lambda}{2\lambda+2\mu}\left(\sigma_{xx} + \sigma_{yy}\right) = \frac{\lambda}{2\lambda+2\mu} \tau_0.
\end{equation}

The strain components for plane strain are
\begin{align}
\epsilon_{xx} &= \frac{1}{2\mu}\left(\sigma_{xx} - \frac{\lambda}{3\lambda+2\mu}\left(\sigma_{xx}+\sigma_{yy}+\sigma_{zz}\right)\right), \\
%
\epsilon_{yy} &= \frac{1}{2\mu}\left(\sigma_{yy} - \frac{\lambda}{3\lambda+2\mu}\left(\sigma_{xx}+\sigma_{yy}+\sigma_{zz}\right)\right), \\
%
\epsilon_{xy} &= \frac{1}{2\mu} \sigma_{xy}.
\end{align}
%
Substituting in our expressions for the stress components leads to
\begin{align}
\epsilon_{xx} &= \frac{\tau_0}{4\mu}\frac{\lambda+2\mu}{\lambda+\mu}, \\
\epsilon_{yy} &= -\frac{\tau_0}{4\mu}\frac{\lambda}{\lambda+\mu}, \\
\epsilon_{xy} &= 0.
\end{align}
The corresponding displacement field is
\begin{align}
u_x &= \frac{\tau_0}{4\mu}\frac{\lambda+2\mu}{\lambda+\mu} x, \\
u_y &= -\frac{\tau_0}{4\mu}\frac{\lambda}{\lambda+\mu} y.
\end{align}

Putting everything together the solution to our boundary value problem is
\begin{align}
\sigma_{xx} &= \tau_0, \\
\sigma_{yy} &= 0, \\
\sigma_{xy} &= 0, \\
\sigma_{zz} &= \frac{\lambda}{2\lambda+2\mu} \tau_0, \\
\epsilon_{xx} &= \frac{\tau_0}{4\mu}\frac{\lambda+2\mu}{\lambda+\mu}, \\
\epsilon_{yy} &= -\frac{\tau_0}{4\mu}\frac{\lambda}{\lambda+\mu}, \\
\epsilon_{xy} &= 0, \\
u_x &= \frac{\tau_0}{4\mu}\frac{\lambda+2\mu}{\lambda+\mu} x, \\
u_y &= -\frac{\tau_0}{4\mu}\frac{\lambda}{\lambda+\mu} y.
\end{align}

:::{figure-md} fig:stress:function:axial:extension
<img src="figs/airy_stress_axial_extension.*" alt="Geometry for Airy stress function axial extension example." width="400px" />

Diagram of axial extension example for the Airy stress function.
We constraint the normal degree of freedom on the $x=x_0$ and $y=y_0$ boundaries and apply a uniform normal traction $\tau_0$ on $x=x_1$.
:::
