(sec-user-governingeqns-elasticity-genmaxwell)=
# Generalized Maxwell Viscoelastic Models

The most general form of linear viscoelastic model is the generalized Maxwell model, which can be represented by a spring in parallel with a number of Maxwell models (see {numref}`fig:viscoelasticity:models`).
With this model it is possible to represent a number of simpler viscoelastic models.
For example, a simple Maxwell model corresponds to setting the elastic constants of all springs to zero, with the exception of the spring contained in the first Maxwell model ($\mathit{\mu}_{1}$).
Similarly, the Kelvin-Voigt model corresponds to setting the elastic constants $\mathit{\mu}_{2} = \mathit{\mu}_{3} = 0$, and setting $\mathit{\mu}_{1} = \infty$ (or a very large number).

:::{admonition} TODO
:class: error

Add more information about using the generalized Maxwell model to implement various viscoelastic models (Burgers material).
:::

We follow formulations similar to those used by {cite}`Zienkiewicz:Taylor:2000` and {cite}`Taylor:2003`.
In this formulation, we specify the total shear modulus of the model ($\mu_{tot}$) and the bulk modulus ($K$).
We then provide the fractional shear modulus for each Maxwell element spring in the model.
It is not necessary to specify the fractional modulus for $\mu_{0}$, because we can find it by subtracting the sum of the other ratios from 1.
Note that the sum of all these fractions must equal 1.
We use a similar formulation for our linear Maxwell viscoelastic model, but in that case $\mu_{0}$ is always zero and we only use a single Maxwell model.
The parameters defining the both materials are listed in {numref}`tab:elasticity:auxiliary:subfields`.

As for all our viscoelastic models, the volumetric strain is completely elastic, and the viscoelastic deformation may be expressed purely in terms of the deviatoric components:
%
```{math}
:label: eqn:genmax:stress-strain
\boldsymbol{\sigma}^{dev} = 2\mu_{tot}\left[\mu_{0}\boldsymbol{\epsilon}^{dev} + \sum_{i=1}^{N}\mu_{i}\boldsymbol{q}^{i} - \boldsymbol{\epsilon}^{refdev}\right] + \boldsymbol{\sigma}^{refdev}; P = 3K(\theta - \theta^{ref}) + P^{ref},
```
%
where $K$ is the bulk modulus, $N$ is the number of Maxwell models, the terms with $ref$ superscripts refer to reference states, and the variable $q_{i}$ follows the evolution equations
%
```{math}
:label: eqn:genmax:statevar:evolution
\boldsymbol{\dot{q}}^{i} + \frac{1}{\tau_{i}}\boldsymbol{q}^{i} = \boldsymbol{\dot{\epsilon}}^{dev}.
```
%
The $\tau_{i}$ are the relaxation times for each Maxwell model:
%
```{math}
:label: eqn:genmax:relaxtime
\tau_{i} = \frac{\eta_{i}}{\mu_{tot}\mu_{i}}.
```

An alternative to the differential equation form above is an integral equation form expressed in terms of the relaxation modulus function.
This function is defined in terms of an idealized experiment in which, at time zero ($t = 0$), a specimen is subjected to a constant strain, $\boldsymbol{\epsilon}^{dev}_{0}$, and the stress response, $\boldsymbol{\sigma}^{dev}(t)$, is measured. For a linear material we obtain:
%
```{math}
:label: eqn:genmax:relaxmodulus:fn
\boldsymbol{\sigma}(t) = 2\mu(t)(\boldsymbol{\epsilon_{0}^{dev} - \boldsymbol{\epsilon}^{refdev}}) + \boldsymbol{\sigma}^{refdev},
```
%
where $\mu(t)$ is the shear relaxation modulus function. Using linearity and superposition for an arbitrary state of strain yields an integral equation:
%
```{math}
:label: eqn:genmax:relaxmodulus:integral
\boldsymbol{\sigma}^{dev}(t) = \intop_{-\infty}^{t}\mu(t - T)\boldsymbol{\dot{\epsilon}}^{dev} dT.
```

Writing the modulus function in Prony series form we obtain
%
```{math}
:label: eqn:genmax:relaxmodulus:prony
\mu(t) = \mu_{tot}\left(\mu_{0} + \sum_{i=1}^{N}\mu_{i}\exp\frac{-t}{\tau_{i}}\right),
```
%
where
%
```{math}
:label: eqn:genmax:shearfraction
\mu_{0} + \sum_{i=1}^{N}\mu_{i} = 1.
```
%
With the form in {math:numref}`eqn:genmax:relaxmodulus:prony`, the integral equation form is identical to the differential equation form.

If we assume the material is undisturbed until a strain is suddenly applied at time zero, we can divide the integral into
%
```{math}
:label: eqn:genmax:integral:timedivide
\intop_{-\infty}^{t} (\cdot) dT = \intop_{-\infty}^{0^{-}} (\cdot) dT + \intop_{0^{-}}^{0^{+}} (\cdot) dT + \intop_{0^{+}}^{t} (\cdot) dT .
```
%
The first term is zero, the second term includes a jump term associated with $\boldsymbol{\epsilon}^{dev}_{0}$ at time zero, and the last term covers the subsequent history of strain.
Applying this separation to {math:numref}`eqn:genmax:relaxmodulus:integral`,
%
```{math}
:label: eqn:genmax:stressintegral
\boldsymbol{\sigma}^{dev}(t) = 2\mu(t)(\boldsymbol{\epsilon}^{dev}_{0} - \boldsymbol{\epsilon}^{refdev}) + \boldsymbol{\sigma}^{refdev} + 2\intop_{0}^{t}\mu(t - T) \boldsymbol{\dot{\epsilon}}^{dev}(T) dT,
```
%
where we have left the sign off of the lower limit on the integral.
%
Substituting {math:numref}`eqn:genmax:relaxmodulus:prony` into {math:numref}`eqn:genmax:stressintegral`, we obtain
%
```{math}
:label: eqn:genmax:sum:integral
\boldsymbol{\sigma}^{dev}(t) = 2\mu_{tot}\left\{ \mu_{0}\boldsymbol{\epsilon}^{dev}(t) + \sum_{i=1}^{N}\left[\mu_{i}\exp\frac{-t}{\tau_{i}}\left(\boldsymbol{\epsilon}_{0}^{dev} + \intop_{0}^{t}\exp\frac{t}{\tau_{i}} \boldsymbol{\dot{\epsilon}}^dev(T) dT\right)\right] - \boldsymbol{\epsilon^{refdev}}\right\} + \boldsymbol{\sigma}^{refdev}.
```
%
We then split each integral into two ranges: from $0$ to $t_{n}$, and from $t_{n}$ to $t$, and define each integral as
%
```{math}
:label: eqn:genmax:integral:part
\boldsymbol{i}_{i}^{1}(t) = \intop_{0}^{t}\exp\frac{T}{\tau_{i}}\boldsymbol{\dot{\epsilon}}^{dev}(T) dT.
```
%
The integral then becomes
%
```{math}
:label: eqn:genmax:integral:part1
\boldsymbol{i}_{i}^{1}(t) = \boldsymbol{i}_{i}^{1}(t_{n}) + \intop_{t_{n}}^{t}\exp\frac{T}{\tau_{i}}\boldsymbol{\dot{\epsilon}}^{dev}(T) dT.
```
%
Including the negative exponential multiplier, we have
%
```{math}
:label: eqn:genmax:integral:negexp
\boldsymbol{h}_{i}^{1}(t) = \exp\frac{-t}{\tau_{i}}\boldsymbol{i}_{i}^{1}.
```
%
Then
%
```{math}
:label: eqn:genmax:integral:hdefine
\boldsymbol{h}_{i}^{1}(t) = \exp\frac{-\Delta t}{\tau_{i}}\boldsymbol{h}_{i}^{1}(t_{n}) + \Delta \boldsymbol{h}_{i}.
```
%
where
%
```{math}
:label: eqn:genmax:integral:deltahdefine
\Delta \boldsymbol{h}_{i} = \exp\frac{-t}{\tau_{i}}\intop_{t_{n}}^{t}\exp\frac{T}{\tau_{i}}\boldsymbol{\dot{\epsilon}}^{dev}(T) dT.
```
%
Approximating the strain rate as constant over each time step, the solution may be found as
%
```{math}
:label: eqn:genmax:integral:deltahsoln
\Delta \boldsymbol{h}_{i} = \frac{\tau_{i}}{\Delta t}\left(1 - \exp\frac{-\Delta t}{\tau_{i}}\right)\left(\boldsymbol{\epsilon}^{dev} - \boldsymbol{\epsilon}_{n}^{dev}\right) = \Delta h_{i}\left(\boldsymbol{\epsilon}^{dev} - \boldsymbol{\epsilon}_{n}^{dev}\right).
```

The approximation is singular for zero time steps, but a series expansion may be used for small time-step sizes:
%
```{math}
:label: eqn:genmax:integral:deltahapprox
\Delta h_{i}\approx1-\frac{1}{2}\left(\frac{\Delta t}{\tau_{i}}\right)+\frac{1}{3!}\left(\frac{\Delta t}{\tau_{i}}\right)^{2}-\frac{1}{4!}\left(\frac{\Delta t}{\tau_{i}}\right)^{3}+\cdots\,.
```
%
This converges with only a few terms.
With this formulation, the constitutive relation now has the simple form:
%
```{math}
:label: eqn:genmax:soln
\boldsymbol{\sigma}^{dev}(t) = 2\mu_{tot}\left( \mu_{0}\boldsymbol{\epsilon}^{dev}(t) + \sum_{i=1}^{N}\mu_{i}\boldsymbol{h}_{i}^{1}(t) - \boldsymbol{\epsilon^{refdev}}\right) + \boldsymbol{\sigma}^{refdev}.
```

We need to compute the tangent constitutive matrix when forming the Jacobian matrix.
In addition to the volumetric contribution to the tangent constitutive matrix, we require the deviatoric part:
%
```{math}
:label: eqn:genmax:devpart
\frac{\partial\sigma_{ij}^{dev}}{\partial\epsilon_{kl}} = \frac{\partial\sigma_{ij}^{dev}}{\partial\epsilon_{mn}^{dev}}\frac{\partial\epsilon_{mn}^{dev}}{\partial\epsilon_{kl}},
```
%
where the second derivative on the right may be easily deduced from the definition of deviatoric strain in {ref}`tab:viscoelasticity:notation`. The other derivative is given by
%
```{math}
:label: eqn:genmax:devjacobian
\frac{\partial\sigma_{ij}^{dev}}{\partial\epsilon_{mn}^{dev}} = 2\mu_{tot} \left[\mu_{0}\delta_{im}\delta_{jn} + \sum_{l=1}^{N}\mu_{l}\frac{\partial h_{lij}^{1}}{\partial \epsilon_{mn}^{dev}}\right].
```
%
From {math:numref}`eqn:genmax:integral:hdefine` and {math:numref}`eqn:genmax:integral:deltahsoln`, the derivative inside the brackets is
%
```{math}
:label: eqn:genmax:bracketderiv
\frac{\partial h_{lij}^{1}}{\partial\epsilon_{mn}^{dev}} = \Delta h_{l}(\Delta t)\delta_{im}\delta_{jn}.
```
%
The complete deviatoric tangent relation is then
%
```{math}
:label: eqn:genmax:jacobian
\frac{\partial\sigma_{ij}^{dev}}{\partial\epsilon_{mn}} = 2\mu_{tot} \left[\mu_{0} + \sum_{l=1}^{N}\mu_{l}\Delta h_{l}(\Delta t)\right]\frac{\partial \epsilon_{ij}^{dev}}{\partial \epsilon_{mn}}.
```

We use this formulation for both our Maxwell and generalized Maxwell viscoelastic models.
For the Maxwell model, $\mu_{0} = 0$ and $N = 1$.
For the generalized Maxwell model, $N = 3$.
The stable time step is equal to 1/5 of the minimum relaxation time for all of the Maxwell models {math:numref}`eqn:genmax:relaxtime`.
