## R CMD check results

0 errors | 0 warnings | 0 note

## Resubmission
This is a resubmission. In this version I have:

* used message()/warning() instead of of print()/cat() if possible,

* used if(verbose)cat(..) if I really had to write text to the console,

* ensured that I do not use more than 2 cores in your examples, vignettes, etc.
