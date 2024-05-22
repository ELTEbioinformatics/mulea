## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## urlchecket::url_check() results

Error: SSL certificate problem: unable to get local issuer certificate
https://regulondb.ccg.unam.mx/

My understanding is that the website should update its certificate but I did not
find any direct contact to the website maintainers where I could request this.

Error: 403: Forbidden
https://doi.org/10.1073/pnas.0506580102

The other doi.org urls pass the checks. The url seems to work.

## Comments from CRAN Team

"Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'mulea', ...
Please note that package names are case sensitive."

mulea has been replaced with 'mulea' within the description section.

"Please always explain all acronyms in the description text. -> FDR"

While "FDR" is mentioned without explanation in the title, it is then explained
in the description section: "'mulea' employs an innovative empirical false
discovery rate (eFDR) correction method". If possible, I would like to keep the
abbreviation as is in the title. Other packages on CRAN follow similar approach
e.g. 'fftw' (https://cran.r-project.org/web/packages/fftw/index.html).

"\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user. Does not seem necessary.
Please replace \dontrun with \donttest.

Please unwrap the examples if they are executable in < 5 sec, or replace
dontrun{} with \donttest{}."

\dontrun{} has been replaced with \donttest{}
