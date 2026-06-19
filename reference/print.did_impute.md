# Print a did_impute result

Print a did_impute result

## Usage

``` r
# S3 method for class 'did_impute'
print(x, ...)
```

## Arguments

- x:

  a `did_impute` object.

- ...:

  ignored.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
panel <- expand.grid(i = 1:4, t = 1:6)
panel$Ei <- ifelse(panel$i <= 2, 4L, NA_integer_)
panel$y  <- 0.3 * (!is.na(panel$Ei) & panel$t >= panel$Ei) + rnorm(nrow(panel), sd = 0.1)
res <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei")
#> WARNING: suppressing wtr, consider lower minn or minn=0.
#> The number of treated entities is too small for some cohorts. Standard Errors may be wrong, consider using avgeffectsby option, averaging the the effect by treated X post variable.
print(res)
#> <didimpute did_impute result>
#>   Effects: tau_ate=0 
#>   SE     : tau_ate=0 
#>   N obs: 18 
```
