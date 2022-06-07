(sec-user-governingeqns-bulk-elasticity-constitutive)=
# Elasticity Constitutive Models

The Jacobian for the elasticity equation {math:numref}`eqn:elasticity:quasistatic:jacobian:pointwise` is
%
```{math}
:label: eqn:elasticity:jacobian
J_{f3}^{uu} = \frac{\partial F^{u_{i}}}{\partial u_{j}}.
```
%
In computing the derivative, we consider the linearized form:
%
```{math}
:label: eqn:linearized:elasticity
\begin{gather}
\sigma_{ik} = C_{ikjl} \epsilon_{jl} \\
\sigma_{ik} = C_{ikjl} \frac{1}{2}\left(u_{j,l} + u_{l,j}\right) \\
\sigma_{ik} = \frac{1}{2}\left(C_{ikjl} + C_{iklj}\right)u_{j,l} \\
\sigma_{ik} = C_{ikjl} u_{j,l}. \\
\end{gather}
```
%
In computing the Jacobian, we take the derivative of the stress tensor with respect to the displacement field,
%
```{math}
:label: eqn:elasticity:jacobian:deriv
\frac{\partial}{\partial u_{j}}\sigma_{ik} = C_{ikjl} {\psi_\mathit{basis^{}}^{u}}_{j,l},
```
%
so we have
%
```{math}
:label: eqn:elasticity:jacobian:cikjl
J_{f3}^{uu}(i,j,k,l) = -C_{ikjl}.
```
%
For many elasticity constitutive models we prefer to separate the stress into the mean stress and deviatoric stress:
%
```{math}
:label: eqn:elasticity:devstress
\begin{gather}
\boldsymbol{\sigma} = \sigma^{mean} \boldsymbol{I} + \boldsymbol{\sigma}^{dev}, \text{where} \\
\sigma^{mean} = \frac{1}{3}\mathrm{Tr}(\boldsymbol{\sigma}) = \frac{1}{3}\left(\sigma_{11} + \sigma_{22} + \sigma_{33}\right). \\
\end{gather}
```
%
Sometimes it is convenient to use pressure (positive pressure corresponds to compression) instead of the mean stress:
%
```{math}
:label: eqn:elasticity:devstress:pressure
\begin{gather}
\boldsymbol{\sigma} = -p \boldsymbol{I} + \boldsymbol{\sigma}^{dev}, \text{where} \\
p = -\frac{1}{3}\mathrm{Tr}(\boldsymbol{\sigma}). \\
\end{gather}
```
%
The Jacobian with respect to the deviatoric stress is
%
```{math}
:label: eqn:elasticity:devstress:jacobian
\begin{gather}
\frac{\partial \sigma_{ik}^{dev}}{\partial u_{j}} = \frac{\partial}{\partial u_{j}}\left(\sigma_{ik} - \frac{1}{3}\sigma_{mm}\delta_{ik}\right) \\
\frac{\partial \sigma_{ik}^{dev}}{\partial u_{j}} = C_{ikjl} {\psi_\mathit{basis^{}}^{u}}_{j,l} - \frac{1}{3}C_{mmjl} \delta_{ik}{\psi_\mathit{basis^{}}^{u}}_{j,l}.
\end{gather}
```
%
We call these modified elastic constants $C_{ijkl}^{dev}$, so that we have
%
```{math}
:label: eqn:elasticity:devstress:jacobian2
C_{ikjl}^{dev} = C_{ikjl} - \frac{1}{3}C_{mmjl} \delta_{ik}.
```
