# mulea 1.1.2

## Bug Fixes

- The resampling eFDR's parallel backend now propagates the host's library paths to
  the PSOCK worker processes before loading the package. This fixes vignette
  re-building under `R CMD check` on Windows, which failed with "there is no package
  called 'mulea'" because the workers could not locate the package in the temporary
  `*.Rcheck` library.

# mulea 1.1.1

- Updated citation to the BMC Bioinformatics article.

# mulea 1.1.0

## Minor Improvements

- Random seed argument was introduced to the ora function.
- Range in internal function `addExtraParamsToPlt()` was updated from 0 to 0.2.

# mulea 1.0.1

## Bug Fixes

- CRAN check results indicated some OS dependent behaviour. This is now fixed.

# mulea 1.0.0

- First major release

# Version 0.99.10 (2024-02-26)

- NEWS file was created.
