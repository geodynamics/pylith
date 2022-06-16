(sec-user-governingeqns-elasticity-powerlaw)=
# Power-law Viscoelastic Models

Laboratory results on rock rheology are typically performed using a triaxial experiment, and the creep data are fit to a power-law equation of the form (e.g., {cite}`Kirby:Kronenberg:1987`:
%
```{math}
:label: eqn:powerlaw:triax
\dot{\epsilon}_{11}^{creep} = A_{E} \exp\left(\frac{-Q}{RT}\right) \left(\sigma_{1} - \sigma_{3}\right)^{n} = A_{E} \exp\left(\frac{-Q}{RT}\right) \sigma_{d}^{n},
```
%
where $\dot{\epsilon}_{11}^{creep}$ is the strain rate in the direction of the maximum principal stress $(\sigma_{1})$, $A_{E}$ is the experimentally-derived pre-exponential constant, $Q$ is the activation enthalpy, $R$ is the universal gas constant, $T$ is the absolute temperature, $n$ is the power-law exponent, $\sigma_{3} (= \sigma_{2})$ is equal to the confining pressure, and $\sigma_{d}$ is the differential stress.
To properly formulate the flow law, it must be generalized so that the results are not influenced by the experiment type or the choice of coordinate systems (e.g., {cite}`Paterson:1994`).
The flow law may then be generalized in terms of the deviatoric stress and strain rate invariants:
%
```{math}
:label: eqn:powerlaw:flowlaw:invars
\sqrt{\dot{L}_{2}^{\prime creep}} = A_{M} \exp \left(\frac{-Q}{RT}\right)\sqrt{J_{2}^{\prime}}^{n},
```
%
where $A_{M}$ is now a pre-exponential constant used in the formulation for modeling.
In practice, it is necessary to compute each strain rate component using the flow law.
This is accomplished using:
%
```{math}
:label: eqn:powerlaw:flowlaw:components
\dot{\epsilon}_{ij}^{creep} = A_{M} \exp \left(\frac{-Q}{RT}\right)\sqrt{J_{2}^{\prime}}^{n-1}\sigma_{ij}^{dev}.
```
%
Note that {math:numref}`eqn:powerlaw:flowlaw:invars` and {math:numref}`eqn:powerlaw:flowlaw:components` are consistent, since {math:numref}`eqn:powerlaw:flowlaw:invars` may be obtained from {math:numref}`eqn:powerlaw:flowlaw:components` by taking the scalar inner product of both sides, multiplying by 1/2, and taking the square root.

In a triaxial experiment with confining pressure $P_{c}$, we have
%
```{math}
:label: eqn:powerlaw:triax2
\begin{gathered}
\sigma_{2} = \sigma_{3} = P_{c} \\
\sigma_{1} = \sigma_{1}^{app} \\
P = \frac{\sigma_{1} + 2P_{c}}{3}, \\
\end{gathered}
```
%
where $\sigma_{1}^{app}$ is the applied load.
The deviatoric stresses are then:
%
```{math}
:label: eqn:powerlaw:triax:devstress
\begin{gathered}
\sigma_{1}^{dev} = \frac{2}{3} \left(\sigma_{1} - P_{c}\right) \\
\sigma_{2}^{dev} = \sigma_{3}^{dev} = -\frac{1}{3} \left(\sigma_{1} - P_{c}\right). \\
\end{gathered}
```
%
This gives
%
```{math}
:label: eqn:powerlaw:triax:devstress2
\begin{gathered}
\sigma_{1}^{dev} = \frac{2}{3} \left(\sigma_{1} - \sigma_{3}\right) = \frac{2}{3} \sigma_{d} \\
\sigma_{2}^{dev} = \sigma_{3}^{dev} = -\frac{1}{3} \left(\sigma_{1} - \sigma_{3}\right) = -\frac{1}{3} \sigma_{d}. \\
\end{gathered}
```
%
In terms of the second deviatoric stress invariant, we then have
%
```{math}
:label: eqn:powerlaw:triax:invar
\sqrt{J_{2}^{\prime}} = \frac{\sigma_{d}}{\sqrt{3}}.
```

Under the assumption that the creep measured in the laboratory experiments is incompressible, we have
%
```{math}
:label: eqn:powerlaw:triax:strainrate
\begin{gathered}
\dot{\epsilon}_{11}^{creepdev} = \dot{\epsilon}_{11} \\
\dot{\epsilon}_{22}^{creepdev} = \dot{\epsilon}_{33} = -\frac{1}{2}\dot{\epsilon}_{11}. \\
\end{gathered}
```
%
In terms of the second deviatoric strain rate invariant we then have
%
```{math}
:label: eqn:powerlaw:triax:strainrate:invar
\sqrt{\dot{L}_{2}^{\prime creep}} = \frac{\sqrt{3}}{2} \dot{\epsilon}_{11}.
```
%
Substituting {math:numref}`eqn:powerlaw:triax:invar` and {math:numref}`eqn:powerlaw:triax:strainrate:invar` into {math:numref}`eqn:powerlaw:triax`, we obtain
%
```{math}
:label: eqn:powerlaw:triax:invars
\sqrt{\dot{L}_{2}^{\prime creep}} = A_{E} \frac{\sqrt{3}^{n+1}}{2} \exp\left(\frac{-Q}{RT}\right)\sqrt{J_{2}^{\prime}}^{n},
```
%
and therefore,
%
```{math}
:label: eqn:powerlaw:consts
A_{M} = \frac{\sqrt{3}^{n+1}}{2}A_{E}.
```

When the exponential factor is included, we define a new parameter:
%
```{math}
:label: eqn:powerlaw:const:exponential
A_{T} = A_{M}\exp\left(\frac{-Q}{RT}\right) = \frac{\sqrt{3}^{n+1}}{2}A_{E} \exp\left(\frac{-Q}{RT}\right).
```

There is a problem with the usage of parameters $A_{E}$, $A_{M}$, and $A_{T}$.
The dimensions of these parameters depend on the value of the power-law exponent; they are not really constants.
In addition to being logically inconsistent, this presents problems when specifying parameters for PyLith, since the power-law exponent must be known before the units can be determined.
An alternative way of writing the flow rule is (e.g., {cite}`Prentice:1968`):
%
```{math}
:label: eqn:powerlaw:flowrule:alternate
\frac{\sqrt{\dot{L}_{2}^{\prime creep}}}{\dot{\epsilon}_{0}^{dev}} = \left(\frac{\sqrt{J_{2}^{\prime}}}{\sigma_{0}^{dev}}\right)^{n},
```
%
where $\dot{\epsilon}_{0}^{dev}$ and $\sigma_{0}^{dev}$ are reference values for the deviatoric strain rate and stress.
This means that
%
```{math}
:label: eqn:powerlaw:flowrule:alternate:consts
\frac{\dot{\epsilon}_{0}^{dev}}{\sigma_{0}^{dev}} = A_{T}.
```

Users must therefore specify three parameters for a power-law material.
The properties `power_law_reference_strain_rate`, `power_law_reference_stress` and `power_law_exponent` in {ref}`tab:elasticity:auxiliary:subfields` refer to $\dot{\epsilon}_{0}^{dev}$, $\sigma_{0}^{dev}$, and $n$, respectively.
To specify the power-law properties for PyLith using laboratory results, the user must first compute $A_{T}$ using {math:numref}`eqn:powerlaw:const:exponential`.
Then, values for $\dot{\epsilon}_{0}^{dev}$ and $\sigma_{0}^{dev}$ must be provided.
The simplest method is probably to assume a reasonable value for the reference strain rate, and then compute $\sigma_{0}^{dev}$ as
%
```{math}
:label: eqn:powerlaw:flowrule:refstress
\sigma_{0}^{dev} = \left(\frac{\dot{\epsilon}_{0}^{dev}}{A_{T}}\right)^{\frac{1}{n}}.
```

We provide (`powerlaw_gendb.py`) to convert laboratory results to the properties used by PyLith.
To use the code, you must specify the spatial variation of $A_{E}$, $Q$, $n$, and $T$.
An additional parameter is given to define the units of $A_{E}$.
You must also specify either a reference stress or a reference strain rate.

The flow law in component form is
%
```{math}
:label: eqn:powerlaw:flowlaw:components2
\dot{\epsilon}_{ij}^{creepdev} = \frac{\dot{\epsilon}_{0}^{dev}\sqrt{J_{2}^{\prime}}^{n-1}\sigma_{ij}^{dev}}{\left(\sigma_{0}^{dev}\right)^{n}},
```
%
and the creep strain increment is approximated as
%
```{math}
:label: eqn:powerlaw:creepincr
\boldsymbol{\Delta \epsilon}^{creepdev} \approx \frac{\Delta t \dot{\epsilon}_{0}^{creepdev}\sqrt{J_{2}^{\prime}(\tau)}^{n-1}\boldsymbol{\sigma}^{dev}(\tau)}{\left(\sigma_{0}^{dev}\right)^{n}} = \frac{\Delta t \dot{\epsilon}_{0}^{creepdev}\sqrt{3}\overline{\sigma}(\tau)^{n-1}\boldsymbol{\sigma}^{dev}(\tau)}{\left(\sqrt{3}\sigma_{0}^{dev}\right)^{n}}.
```
%
Therefore,
%
```{math}
:label: eqn:powerlaw:creepincr:scalar
\Delta \overline{\epsilon}^{creepdev} \approx \frac{2 \Delta t \dot{\epsilon}_{0}^{creepdev}\sqrt{J_{2}^{\prime}(\tau)}^{n}}{\sqrt{3}\left(\sigma_{0}^{dev}\right)^{n}} = \frac{2 \Delta t \dot{\epsilon}_{0}^{creepdev}\overline{\sigma}(\tau)^{n}}{\sqrt{3}^{n+1}\left(\sigma_{0}^{dev}\right)^{n}},
```
%
and
%
```{math}
:label: eqn:powerlaw:gamma
\gamma(\tau) = \frac{\dot{\epsilon}_{0}^{creepdev}\sqrt{J_{2}^{\prime}(\tau)}^{n-1}}{\left(\sigma_{0}^{dev}\right)^{n}}.
```
%
Combining {math:numref}`eqn:effstress:stress-alpha`, {math:numref}`eqn:powerlaw:creepincr`, {math:numref}`eqn:powerlaw:creepincr:scalar`, {math:numref}`eqn:powerlaw:gamma`, and {math:numref}`eqn:effstress:stress-strain2`, we obtain:
%
```{math}
:label: eqn:powerlaw:stress-alpha
\boldsymbol{\sigma}^{dev}(t + \Delta t) = \frac{1}{a_{E}} \left\{\boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) - \Delta t \gamma(\tau)\left[(1-\alpha)\boldsymbol{\sigma}^{dev}(t) + \alpha \boldsymbol{\sigma}^{dev}(t + \Delta t)\right]\right\} + \boldsymbol{\sigma}^{refdev},
```
%
which we can rewrite as
%
```{math}
:label: eqn:powerlaw:stress-alpha:alternate
\boldsymbol{\sigma}^{dev}(t + \Delta t) \left[a_{E} + \alpha \Delta t \gamma(\tau)\right] = \boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) - \Delta t \gamma(\tau)(1-\alpha)\boldsymbol{\sigma}^{dev}(t) + a_{E} \boldsymbol{\sigma}^{refdev}.
```
%
Taking the scalar inner product of both sides we obtain:
%
```{math}
:label: eqn:powerlaw:effstressfn
a^{2}J_2^{\prime}(t + \Delta t) - b + c \gamma(\tau) - d^{2}\gamma(\tau)^{2} = F = 0,
```
%
where
%
```{math}
:label: eqn:powerlaw:effstressfn:variables
\begin{gathered}
a = a_{E} + \alpha \Delta t \gamma(\tau), \\
b = \frac{1}{2} \boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) : \boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) + a_{E} \boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) : \boldsymbol{\sigma}^{refdev} + a_{E}^{2} J_{2}^{\prime ref}, \\
c = \Delta t(1 - \alpha) \boldsymbol{\epsilon}^{\prime dev}(t + \Delta t) : \boldsymbol{\sigma}^{dev}(t) + \Delta t(1 - \alpha) a_{E} \boldsymbol{\sigma}^{dev}(t) : \boldsymbol{\sigma}^{refdev}, \\
d = \Delta t (1 - \alpha) \sqrt{J_{2}^{\prime}(t)}. \\
\end{gathered}
```

Equation {math:numref}`eqn:powerlaw:effstressfn` is a function of a single unknown -- the square root of the second deviatoric stress invariant at time $t + \Delta t$ -- which we can solve by bisection or by Newton’s method.
Then we compute the deviatoric stresses for the current time step from {math:numref}`eqn:effstress:effstresstau`, {math:numref}`eqn:powerlaw:gamma`, and {math:numref}`eqn:powerlaw:stress-alpha`.
To compute the total stress we combine the deviatoric and volumetric components from {math:numref}`eqn:effstress:stress-strain`.

To compute the tangent stress-strain relation, we first rewrite {math:numref}`eqn:powerlaw:stress-alpha:alternate` as
%
```{math}
:label: eqn:powerlaw:stress-alpha:function
\begin{gathered}
\sigma_{ij}^{dev}(t + \Delta t) \left[a_{E} + \alpha \Delta t \gamma(\tau)\right] - \epsilon_{ij}^{\prime dev}(t + \Delta t) + \sigma_{ij}^{dev}(t)\Delta t \gamma(\tau)(1-\alpha) - a_{E} \sigma_{ij}^{refdev} = 0, \\
\sigma_{ij}^{dev}(t + \Delta t)G  - \epsilon_{ij}^{\prime dev}(t + \Delta t) + \sigma_{ij}^{dev}(t)H - a_{E} \sigma_{ij}^{refdev} = 0. \\
\end{gathered}
```
%
Taking derivatives with respect to $\epsilon_{kl}$, we obtain
%
```{math}
:label: eqn:powerlaw:stress-alpha:derivative
G \frac{\partial \sigma_{ij}^{dev}}{\partial \epsilon_{kl}} + \sigma_{ij}^{dev} \frac{\partial G}{\partial \sigma_{ij}^{dev}} \frac{\partial \sigma_{ij}^{dev}}{\partial \epsilon_{kl}}  - \frac{\partial \epsilon_{ij}^{\prime dev}}{\epsilon_{kl}} + \sigma_{ij}^{dev}(t)\frac{\partial H}{\partial \sigma_{ij}^{dev}}\frac{\partial \sigma_{ij}^{dev}}{\partial \epsilon_{kl}} = 0,
```
%
where we have removed the time notation for quantities evaluated at $t + \Delta t$.
Rearranging, we compute the Jacobian for the deviatoric stress components using
%
```{math}
:label: eqn:powerlaw:stress-alpha:jacobian
\frac{\partial \sigma_{ij}^{dev}}{\partial \epsilon_{kl}} =
\frac{\frac{\partial \epsilon_{ij}^{\prime dev}}{\partial
\epsilon_{kl}}}{G + \sigma_{ij}^{dev} \frac{\partial G}{\partial \sigma_{ij}^{dev}} + \sigma_{ij}^{dev}(t)\frac{\partial H}{\partial \sigma_{ij}^{dev}}}.
```
