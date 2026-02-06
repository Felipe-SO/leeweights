# leeweights
## Description
Calculates weights to implement Lee bounds (Lee, 2009) with fractional trimming and discrete conditional monotonicity (Semenova, 2025).

Fractional trimming assigns partial weight to the boundary observation(s) when the trimming proportion falls between discrete jumps of the empirical distribution, most notably in the presence of tied values at the trimming cutoff.

As in Semenova (2025), the direction of monotonicity is allowed to vary across the covariate space. Unlike Semenova's full aproach, however, this command only allows for a discrete partition of the covariate space. As a result, the procedure does not rely on regularization methods or orthogonalization corrections.

## Syntax
`varname`: Outcome for which effects will be estimated.

`select(varname)`: Variable that indicates selection into the sample.

`treatment(varname)`: Binary treatment indicator.

`tight(varlist)`: List of variables on which the covariate space will be partitioned through a full interaction. Used to define regions of conditional monotonicity and intensity of differential attrition.

`suffix(name)`: Suffix for the weight variables. Preceded by `wgt_[lb|ub]_`.
