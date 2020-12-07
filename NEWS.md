# bbreg 2.0.3

* Improved print and summary outputs when one or more covariates are factors.
* Improved initial guesses for the precision parameters.
* Fixed a bug on the predict argument of the bbreg function when using the model "beta" or when using discrimination but the "beta" model was chosen.
* Fixed a bug on vcov function
* Added 'start', 'fitted.values' and 'coefficients' to the returned list by the bbreg function. Also named the vectors and matrices (rows and columns) accordingly.

