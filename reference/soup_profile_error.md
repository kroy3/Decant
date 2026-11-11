# Soup-profile recovery error.

Soup-profile recovery error: how close the estimated ambient profile is
to the true soup. This is the mechanism the structured method is meant
to fix, so report it directly. L1 distance between probability vectors
(0 = perfect).

## Usage

``` r
soup_profile_error(soup_hat, soup_true)
```

## Arguments

- soup_hat, soup_true:

  ambient profiles (normalised internally).
