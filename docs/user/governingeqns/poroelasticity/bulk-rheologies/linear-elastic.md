(sec-user-governingeqns-poroelasticity-isolinearelasticity)=
# Linear Isotropic

We assume linear elasticity for the solid phase, so we have the following formulation for the stress tensor:
%
\begin{equation}
  \boldsymbol{\sigma}(\vec{u},p) = \boldsymbol{C}:\boldsymbol{\epsilon} - \alpha p \boldsymbol{I}
  = \lambda \boldsymbol{I} \epsilon_{v} + 2 \mu  \boldsymbol{\epsilon} - \alpha \boldsymbol{I} p
\end{equation}
%
where $\lambda$ and $\mu$ are Lam&eacute;'s parameters, $\lambda = K_{d} - \frac{2 \mu}{3}$, $\mu$ is the shear modulus, and the volumetric strain is defined as $\epsilon_{v} = \nabla \cdot \vec{u}$.
For isotropic linear elasticity, all components of $C_{ikjl}$ are zero except for:
%
```{math}
\begin{gathered}
C_{1111} = C_{2222} = C_{3333} = \lambda + 2 \mu, \\
C_{1122} = C_{1133} = C_{2233} = \lambda, \\
C_{1212} = C_{2323} = C_{1313} = \mu. \\
\end{gathered}
```
%
The deviatoric elastic constants are:
%
```{math}
\begin{gathered}
C_{1111}^{dev} = C_{2222}^{dev} = C_{3333}^{dev} = \frac{4}{3} \mu, \\
C_{1122}^{dev} = C_{1133}^{dev} = C_{2233}^{dev} = -\frac{2}{3} \mu, \\
C_{1212}^{dev} = C_{2323}^{dev} = C_{1313}^{dev} = \mu. \\
\end{gathered}
```
