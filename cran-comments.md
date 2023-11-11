## 2023.12.1 Release

Initial release - resubmission


Added \donttest{] to examples that have large CPU time.}
No more NOTES from local, win_release() or win_devel() checks other than a
NEW submission.

> If I understand correctly, your \dontrun{} example is wrapped in
> \dontrun{} for its long run time? If it is otherwise executable, please
> wrap it in \donttest{}. That way the example is run if the user calls
> examples() but won't be tested regularly. If you think \dontrun{} is the
> correct wrapper, please let us know why so we can publish your package.
> Otherwise, please fix and resubmit.

Ok...I removed \dontrun and\donttest from the MCMC examples (calls to BTSPAS package).
I've cut down the examples so they run in < 5 seconds on my Mac, 
but they take a litte longer on Windows machines:

2 notes
Examples with CPU (user + system) or elapsed time > 10s  <- win_release()
                       user system elapsed
LP_BTSPAS_fit_NonDiag 17.54   0.04   17.59
LP_BTSPAS_fit_Diag    10.60   0.02   10.83

or

Examples with CPU (user + system) or elapsed time > 10s  <- win_devel()
                       user system elapsed
LP_BTSPAS_fit_NonDiag 14.57   0.02   14.60
LP_BTSPAS_fit_Diag     9.61   0.01   10.28

I added some comments to the examples, indicating that a real analysis likely
needs different MCMC parameter values for the number of iterations etc to 
ensure proper mixing.

I hope this is now ok.




> Please add \value to .Rd files regarding exported methods and explain
> the functions results in the documentation. Please write about the
> structure of the output (class) and also what the output means. (If a
> function does not return a value, please document that too, e.g.
> \value{No return value, called for side effects} or similar)
> Missing Rd-tags:
>      fit_classes.Rd:  \value

Fixed, but this is NOT a useable function, but only documents how I use "classes"
to ensure that objects of different "types" are not mixed in AICc computations, etc.


> \dontrun{} should only be used if the example really cannot be executed
> (e.g. because of missing additional software, missing API keys, ...) by
> the user. That's why wrapping examples in \dontrun{} adds the comment
> ("# Not run:") as a warning for the user.
> Does not seem necessary.
> Please unwrap the examples if they are executable in < 5 sec, or replace
> \dontrun{} with \donttest{}.

Examples take about 30 seconds to run. Reduced number of MCMC iterations
to reduce execution time, but now results are not properly mixed etc.
Added a warning to the example about this. Dealt with problem of un-closed connections
from the BTSPAS package that also prevented examples from running.





Large documentation object has been moved to GitHub repository with link from Vignettes

> Please reduce the length of the title to less than 65 characters.

Done

> Please rather use the Authors@R field and declare Maintainer, Authors
> and Contributors with their appropriate roles with person() calls.
> e.g. something like:
> Authors@R: c(person("Alice", "Developer", role = c("aut", "cre","cph"),
> email = "alice.developer@some.domain.net"),
> person("Bob", "Dev", role = "aut") )

Done.

> If there are references describing the methods in your package, please
> add these in the description field of your DESCRIPTION file in the form
> authors (year) <doi:...>
> authors (year) <arXiv:...>
> authors (year, ISBN:...)
> or if those are not available: <https:...>
> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
> auto-linking. (If you want to add a title as well please put it in
> quotes: "Title")

This package is a consolidation of many, many authors and methods that are 
likely too numerous to mention in the Description file and best left for the individual
functions and/or the (large) vignette.

> Please add \value to .Rd files regarding exported methods and explain
> the functions results in the documentation. Please write about the
> structure of the output (class) and also what the output means. (If a
> function does not return a value, please document that too, e.g.
> \value{No return value, called for side effects} or similar)
> Missing Rd-tags:
>      LP_IS_print.Rd: \value

All functions now have @returns information provided.

> You write information messages to the console that cannot be easily
> suppressed.
> It is more R like to generate objects that can be used to extract the
> information a user is interested in, and then print() that object.
> Instead of print()/cat() rather use message()/warning() or
> if(verbose)cat(..) (or maybe stop()) if you really have to write text to
> the console. (except for print, summary, interactive functions)
> -> R/LP_IS_fit.R; R/LP_SPAS_fit.R

Corrected. In some cases, e.g. calls to SPAS or BTSPAS, these occur deep in the functions
in the other package. In these cases, added an argument *quietly* that temporarily
"sinks" console output and hides it from the user.


## Test environments
* local OS X install, R 4.3.2
* devtools::check_win_release()
* devtools::check_win_devel()
* devtools::check_rhub()

## R CMD check results
There were no ERRORs or WARNINGs. 

There were up to 1 Notes depending on R version/OS used: 

N Maintainer: ‘Carl James Schwarz <cschwarz.stat.sfu.ca@gmail.com>’
New submission

  Yes, this is a NEW package.


## Reverse dependencies

None.
