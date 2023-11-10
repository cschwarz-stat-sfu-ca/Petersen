## 2023.12.1 Release

Initial release - resubmission

> Please add \value to .Rd files regarding exported methods and explain
> the functions results in the documentation. Please write about the
> structure of the output (class) and also what the output means. (If a
> function does not return a value, please document that too, e.g.
> \value{No return value, called for side effects} or similar)
> Missing Rd-tags:
>      fit_classes.Rd:  \value

Fixed, but this is NOT a useable function, but only documents how I use "classes"
to ensure that objects of different "types" are not mixed in AICc computations, etc.
Is there a better way to add a helpful on a topic other than a vignette?


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
from the BTSPAS package.





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
