# bbreg 2.0.5

* Added optimization controls
* Added Efron's pseudo R2 on summary and as an element of the returned list of bbreg
* Minor general fixes.

# bbreg 2.0.4

* Created the function bbreg.fit
* Improved documentation
* Improved file organization
* Improved the optimization algorithm

# bbreg 2.0.3

* Improved print and summary outputs when one or more covariates are factors.
* Improved initial guesses for the precision parameters.
* Fixed a bug on the predict argument of the bbreg function when using the model "beta" or when using discrimination but the "beta" model was chosen.
* Fixed a bug on vcov function
* Added 'start', 'fitted.values' and 'coefficients' to the returned list by the bbreg function. Also named the vectors and matrices (rows and columns) accordingly.
* Created pbessel, qbessel and rbessel functions. Also adjusted dbessel function to accept vectors as input.

# bbreg 2.0.2 

* Some minor bugs were fixed.

# bbreg 2.0.1

* Correction on the database SA.
* Adjust on plot method to allow additional graphical arguments.

# bbreg 2.0.0

* Several 'for' loops were rewritten in vectorized form.
* Several S3 methods implemented.
* Different link functions for the mean and precision parameters added.
* Functions 'bbsummary' and 'eplot' deprecated.
* Vignette written.
* Github repository created.
* Several improvements on the outputs.
* Initial guesses modified.
* Version 1.0.0 written by Vinicius Mayrink, from version 2.0.0 on, all the code modifications, additions, etc., were written by Alexandre B. Simas.

