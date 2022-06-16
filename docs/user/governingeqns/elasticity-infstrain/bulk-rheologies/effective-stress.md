(sec-user-governingeqns-elasticity-effective-stress)=
# Effective Stress Formulation for Viscoelastic Materials

An alternative technique is used for power-law viscoelastic materials, as well as for Drucker-Prager elasto-plastic materials.
This approach, known as the effective stress function formulation {cite}`Kojic:Bathe:1987` is more suitable for such nonlinear behavior, as it simplifies the equations that must be solved.
As for our linear viscoelastic models, the viscous volumetric strains are zero (incompressible flow), and we separate the general stress-strain relationship at time $t + \Delta t$ into deviatoric and volumetric parts:
%
```{math}
:label: eqn:effstress:stress-strain
\begin{gathered}
\boldsymbol{\sigma}^{dev}(t + \Delta t) = 2\mu \left[\boldsymbol{\epsilon}^{dev}(t + \Delta t) - \boldsymbol{\epsilon}^{creepdev}(t + \Delta t) - \boldsymbol{\epsilon}^{refdev}\right] + \boldsymbol{\sigma}^{refdev} \\
\boldsymbol{\sigma}^{dev}(t + \Delta t) = \frac{1}{a_{E}} \left[\boldsymbol{\epsilon}^{dev}(t + \Delta t) - \boldsymbol{\epsilon}^{creepdev}(t + \Delta t) - \boldsymbol{\epsilon}^{refdev}\right] + \boldsymbol{\sigma}^{refdev} \\
P(t + \Delta t) = 3K\left[\theta(t + \Delta t) - \theta^{ref}\right] + P^{ref} = \frac{1}{a_{m}}\left[\theta(t + \Delta t) - \theta^{ref}\right] + P^{ref}, \\
\end{gathered}
```
%
where $\boldsymbol{\epsilon}^{dev}(t + \Delta t)$ is the total deviatoric strain, $\boldsymbol{\epsilon}^{creepdev}(t + \Delta t)$ is the total viscous strain, $\boldsymbol{\epsilon}^{refdev}$ is the reference deviatoric strain, $P(t + \Delta t)$ is the total pressure, $\theta (t + \Delta t)$ is the mean strain evaluated at time $t + \Delta t$ and $\theta^{ref}$ is the reference mean strain.
The reference deviatoric stress and reference pressure are given by $\boldsymbol{\sigma}^{refdev}$ and $P^{ref}$, respectively.
The middle equation in {math:numref}`eqn:effstress:stress-strain` may also be written as
%
```{math}
:label: eqn:effstress:stress-strain2
\boldsymbol{\sigma}^{dev}(t + \Delta t) = \frac{1}{a_{E}} \left[\boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) - \boldsymbol{\Delta \epsilon}^{creepdev} \right] + \boldsymbol{\sigma}^{refdev},
```
%
where
%
```{math}
:label: eqn:effstress:devstrain-prime
\begin{gathered}
\boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) = \boldsymbol{\epsilon}^{dev}(t + \Delta t) - \boldsymbol{\epsilon}^{creepdev}(t) - \boldsymbol{\epsilon}^{refdev}, \\
\boldsymbol{\Delta \epsilon}^{creepdev} = \boldsymbol{\epsilon}^{creepdev}(t + \Delta t) - \boldsymbol{\epsilon}^{creepdev}(t). \\
\end{gathered}
```
%
The creep strain increment is approximated using
%
```{math}
:label: eqn:effstress:creep-incr
\boldsymbol{\Delta \epsilon}^{creepdev} = \Delta t \gamma(\tau)\boldsymbol{\sigma}^{dev}(\tau),
```
%
where, using the $\alpha$-method of time integration,
%
```{math}
:label: eqn:effstress:stress-alpha
\boldsymbol{\sigma}^{dev}(\tau) = (1 - \alpha) \boldsymbol{\sigma}^{dev}(ref, t) + \alpha \boldsymbol{\sigma}^{dev}(ref, t + \Delta t) + \boldsymbol{\sigma}^{refdev} = (1 - \alpha) \boldsymbol{\sigma}^{dev}(t) + \alpha \boldsymbol{\sigma}^{dev}(t + \Delta t),
```
%
and
%
```{math}
:label: eqn:effstress:gammatau
\gamma(\tau) = \frac{3 \Delta \overline{\epsilon}^{creepdev}}{2 \Delta t \overline{\sigma}(\tau)},
```
%
where
%
```{math}
:label: eqn:effstress:effcreepincr
\Delta \overline{\epsilon}^{creepdev} = \sqrt{\frac{2}{3} \boldsymbol{\Delta \epsilon}^{creepdev} : \boldsymbol{\Delta \epsilon}^{creepdev}}
```
%
and
%
```{math}
:label: eqn:effstress:effstresstau
\overline{\sigma}(\tau) = (1 - \alpha) \overline{\sigma}(ref, t) + \alpha \overline{\sigma}(ref, t + \Delta t) + \overline{\sigma}^{ref} = \sqrt{3 J_{2}^{\prime}(\tau)}.
```
