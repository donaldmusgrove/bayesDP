# bayesDP 1.1.0
## Major new features
* Supports two-arm survival analysis via hazard rate comparisons
* Completely revamped summary and print methods to produce better formatted results
* Plot method allows users to specify a `type`
* Added vignettes for each of `bdpbinomial`, `bdpnormal`, and `bdpsurvival`
* Implemented the `fix_alpha` input which allows users to set the historical data weight at `alpha_max`

## Bug fixes and minor improvements
* Fixed error with two-arm analysis where models did not fit if either the current or historical control data were not input
* Changed `two_side` input to logical
* Consolidated several internal functions into a single function for computational efficiency gains


# bayesDP 1.0.3
* README update
* Added plot types
* Added Vignettes
* Added logo
* Improved documentation
* Updated print, summary, plot methods
* Refactored bdpnormal/bdpbinomial

# bayesDP 1.0.2
* Crucial bugfixes

# bayesDP 1.0.1
* User requested bugfixes

# bayesDP 1.0.0
* Initial CRAN release with normal, binomial and survival functions.
