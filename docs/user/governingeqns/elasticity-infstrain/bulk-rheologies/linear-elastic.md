(sec-user-governingeqns-elasticity-isolinearelasticity)=
# Linear Isotropic Elastic Models

We implement isotropic linear elasticity both with and without a reference stress-strain state.
With a linear elastic material it is often convenient to compute the deformation relative to an unknown initial stress-strain state.
Furthermore, when we use an initial undeformed configuration with zero stress and strain, the reference stress and strain are zero, so this presents a simplification of the more general case of the stress-strain state relative to the reference stress-strain state.

Without a reference stress-strain state, we have
%
```{math}
:label: eqn:elasticity:hookeslaw:noref
\sigma_{ij} = \lambda \epsilon_{kk} \delta{ij} + 2 \mu \epsilon_{ij}
```
%
and with a reference stress-strain state, we have
%
```{math}
:label: eqn:elasticity:hookeslaw:ref
\sigma_{ij} = \sigma_{ij}^{ref} + \lambda\left(\epsilon_{kk} - \epsilon_{kk}^{ref}\right) \delta{ij} + 2 \mu \left(\epsilon_{ij} - \epsilon_{ij}^{ref}\right).
```
%
The mean stress is
%
```{math}
:label: eqn:elasticity:meanstress
\begin{gather}
\sigma^{mean} = \frac{1}{3} \sigma_{kk}, \\
\sigma^{mean} = \frac{1}{3} \sigma_{kk}^{ref} + K\left(\epsilon_{kk} - \epsilon_{kk}^{ref}\right), \\
\end{gather}
```
%
where $K = \lambda + 2 \mu/3$ is the bulk modulus.
If the reference stress and reference strain are both zero, then this reduces to
%
```{math}
:label: eqn:elasticity:meanstress:noref
\sigma^{mean} = K \epsilon_{kk}.
```
%
The deviatoric stress is
%
```{math}
:label: eqn:elasticity:devstress
\begin{gather}
\sigma_{ij}^{dev} = \sigma_{ij} - \sigma^{mean} \delta_{ij}, \\
\sigma_{ij}^{dev} = \sigma_{ij}^{ref} + \lambda \left(\epsilon_{kk} -
\epsilon_{kk}^{ref}\right) \delta_{ij} + 2 \mu \left(\epsilon_{ij} -
\epsilon_{ij}^{ref}\right) - \left[\frac{1}{3}\sigma_{kk}^{ref} +
\left(\lambda + \frac{2}{3}\mu\right)\left(\epsilon_{kk}-
\epsilon_{kk}^{ref}\right)\right] \delta_{ij}, \\
\sigma_{ij}^{dev} = \sigma_{ij}^{ref} - \frac{1}{3}\sigma_{kk}^{ref}\delta_{ij} + 2 \mu \left(\epsilon_{ij} - \epsilon_{ij}^{ref}\right) - \frac{2}{3} \mu\left(\epsilon_{kk}- \epsilon_{kk}^{ref}\right) \delta_{ij}. \\
\end{gather}
```
%
For isotropic linear elasticity, all components of $C_{ikjl}$ are zero except for:
%
```{math}
:label: eqn:elasticity:nonzero:cikjl
\begin{gather}
C_{1111} = C_{2222} = C_{3333} = \lambda + 2 \mu, \\
C_{1122} = C_{1133} = C_{2233} = \lambda, \\
C_{1212} = C_{2323} = C_{1313} = \mu. \\
\end{gather}
```
%
The deviatoric elastic constants are:
%
```{math}
:label: eqn:elasticity:nonzero:cikjl:deviatoric
\begin{gather}
C_{1111}^{dev} = C_{2222}^{dev} = C_{3333}^{dev} = \frac{4}{3} \mu, \\
C_{1122}^{dev} = C_{1133}^{dev} = C_{2233}^{dev} = -\frac{2}{3} \mu, \\
C_{1212}^{dev} = C_{2323}^{dev} = C_{1313}^{dev} = \mu. \\
\end{gather}
```
%