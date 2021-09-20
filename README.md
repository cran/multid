
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multid

<!-- badges: start -->

[![Travis build
status](https://app.travis-ci.com/vjilmari/multid.svg?branch=main)](https://app.travis-ci.com/vjilmari/multid/)
<!-- badges: end -->

The goal of multid is to provide tools for regularized measurement of
multivariate differences between two groups (e.g., sex differences).
Regularization via logistic regression variants enables inclusion of
large number of correlated variables in the multivariate set while
providing k-fold cross-validation and regularization to avoid
overfitting.

## Installation

You can install the released version of multid from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("multid")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vjilmari/multid")
```

## Example

This is a basic example which shows you how to measure standardized
multivariate (both Sepal and Petal dimensions, four variables in total)
distance between setosa and versicolor Species in iris dataset. Similar
example is provided with artificially generated multi group data which
are used as separate data folds in the procedure.

``` r
library(multid)
set.seed(91237)

D_regularized(
  data = iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
  mv.vars = c(
    "Sepal.Length", "Sepal.Width",
    "Petal.Length", "Petal.Width"
  ),
  group.var = "Species",
  group.values = c("setosa", "versicolor")
)$D
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       50           50 8.348364       -9.379  1.386586      2.247504
#>      pooled.sd     diff        D
#> [1,]  1.867337 17.72736 9.493393
# Use different partitions of data for regularization and estimation
D_regularized(
  data = iris[iris$Species == "setosa" |
    iris$Species == "versicolor", ],
  mv.vars = c(
    "Sepal.Length", "Sepal.Width",
    "Petal.Length", "Petal.Width"
  ),
  group.var = "Species",
  group.values = c("setosa", "versicolor"),
  size = 35,
  out = TRUE
)$D
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       15           15 7.614128    -8.973335   1.79109      2.353035
#>      pooled.sd     diff        D
#> [1,]  2.091026 16.58746 7.932692
# multigroup model where the group variable is defined as fold.var, and output is produced separately for each group/fold. Groups are also used as folds in the k-fold cross-validation procedure.

# separate sample folds
# generate data for 10 groups
set.seed(34246)
n1 <- 100
n2 <- 10
d <-
  data.frame(
    sex = sample(c("male", "female"), n1 * n2, replace = TRUE),
    fold = sample(x = LETTERS[1:n2], size = n1 * n2, replace = TRUE),
    x1 = rnorm(n1 * n2),
    x2 = rnorm(n1 * n2),
    x3 = rnorm(n1 * n2)
  )
#'
# Fit and predict with same data
D_regularized(
  data = d,
  mv.vars = c("x1", "x2", "x3"),
  group.var = "sex",
  group.values = c("female", "male"),
  fold.var = "fold",
  fold = TRUE,
  rename.output = TRUE
)$D
#>   n.1 n.2           m.1          m.2       sd.1       sd.2  pooled.sd
#> A  53  48  0.0147710988 -0.022086654 0.07263555 0.06499505 0.06911364
#> B  56  53 -0.0019900692 -0.019396118 0.07002969 0.05978853 0.06525374
#> C  39  56 -0.0214165020 -0.008847558 0.07734223 0.06494597 0.07027580
#> D  53  58 -0.0008188678 -0.011421247 0.06352770 0.06461516 0.06409867
#> E  44  52 -0.0030891361  0.005871956 0.06004801 0.06733488 0.06410440
#> F  60  31 -0.0123361550 -0.005902017 0.06477705 0.06218219 0.06391415
#> G  56  51 -0.0137985663 -0.025847991 0.06259123 0.05462635 0.05893284
#> H  47  51 -0.0040611163 -0.014056292 0.06945851 0.07320530 0.07143449
#> I  42  45  0.0176532701 -0.014835959 0.07031247 0.06559995 0.06791388
#> J  48  57 -0.0049288771 -0.011236030 0.05491690 0.06776115 0.06222991
#>           diff          D pooled.sd.1 pooled.sd.2 pooled.sd.total  d.sd.total
#> A  0.036857753  0.5332921  0.06661768   0.0647895      0.06570629  0.56094714
#> B  0.017406049  0.2667441  0.06661768   0.0647895      0.06570629  0.26490690
#> C -0.012568944 -0.1788517  0.06661768   0.0647895      0.06570629 -0.19128983
#> D  0.010602379  0.1654072  0.06661768   0.0647895      0.06570629  0.16136020
#> E -0.008961092 -0.1397890  0.06661768   0.0647895      0.06570629 -0.13638105
#> F -0.006434138 -0.1006685  0.06661768   0.0647895      0.06570629 -0.09792272
#> G  0.012049425  0.2044603  0.06661768   0.0647895      0.06570629  0.18338313
#> H  0.009995176  0.1399209  0.06661768   0.0647895      0.06570629  0.15211902
#> I  0.032489229  0.4783886  0.06661768   0.0647895      0.06570629  0.49446151
#> J  0.006307153  0.1013524  0.06661768   0.0647895      0.06570629  0.09599010
#'
# Different partitions for regularization and estimation for each data fold
D_regularized(
  data = d,
  mv.vars = c("x1", "x2", "x3"),
  group.var = "sex",
  group.values = c("female", "male"),
  fold.var = "fold",
  size = 17,
  out = TRUE,
  fold = TRUE,
  rename.output = TRUE
)$D
#>   n.1 n.2           m.1           m.2      sd.1      sd.2 pooled.sd
#> A  36  31 -0.0008117477 -0.0054231327 0.1783541 0.1495396 0.1656790
#> B  39  36 -0.0447770702 -0.0388386735 0.1615780 0.1473451 0.1549173
#> C  22  39 -0.0676320635  0.0135231815 0.1691657 0.1297412 0.1450075
#> D  36  41  0.0076243457  0.0088113020 0.1375257 0.1663247 0.1535588
#> E  27  35 -0.0719427999  0.0008396515 0.1315043 0.1126188 0.1211645
#> F  43  14 -0.0426691993  0.0379265035 0.1718121 0.1678070 0.1708739
#> G  39  34 -0.0297091831 -0.0394578004 0.1233644 0.1497425 0.1362613
#> H  30  34  0.0376714002  0.0090690130 0.1629470 0.1540525 0.1582751
#> I  25  28  0.0182168066 -0.0295980682 0.1372502 0.1361344 0.1366606
#> J  31  40  0.0042322250  0.0005736008 0.1539483 0.1513596 0.1524905
#>           diff            D pooled.sd.1 pooled.sd.2 pooled.sd.total
#> A  0.004611385  0.027833254   0.1542809   0.1462916       0.1503151
#> B -0.005938397 -0.038332698   0.1542809   0.1462916       0.1503151
#> C -0.081155245 -0.559662487   0.1542809   0.1462916       0.1503151
#> D -0.001186956 -0.007729654   0.1542809   0.1462916       0.1503151
#> E -0.072782451 -0.600691247   0.1542809   0.1462916       0.1503151
#> F -0.080595703 -0.471667611   0.1542809   0.1462916       0.1503151
#> G  0.009748617  0.071543569   0.1542809   0.1462916       0.1503151
#> H  0.028602387  0.180713173   0.1542809   0.1462916       0.1503151
#> I  0.047814875  0.349880403   0.1542809   0.1462916       0.1503151
#> J  0.003658624  0.023992473   0.1542809   0.1462916       0.1503151
#>     d.sd.total
#> A  0.030678130
#> B -0.039506332
#> C -0.539900950
#> D -0.007896456
#> E -0.484199322
#> F -0.536178487
#> G  0.064854561
#> H  0.190282908
#> I  0.318097694
#> J  0.024339704
```
