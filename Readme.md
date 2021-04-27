# blotIt3 - a framework for alignment of biological replicate data

The present package is a rewritten version of [blotIt2](https://github.com/dkaschek/blotIt2) by Daniel Kaschek. The aim of this toolbox is to scale biological replicate data to a common scale, making the quantitative data of different gels comparable.

## System preperation

blotIt3 requires the `R` packages `utils, MASS, data.table, ggplot2, rootSolve` and `trust`. Additionally, the package `devtools` is needed to install blotIt3 from github. If not already done, the required packages can be installed by executing

```r
install.packages(c("utils", "MASS", "data.table", "ggplot2", "rootSolve", "trust", "devtools"))
```
blotIt3 then is installed via `devtools`:
```r
devtools::install_github("SeverinBang/blotIt3")
```
## Usage

### Data import
First, the packe is imported
```r
library(blotIt3)
```
A .csv file is imported and is formated by the function `read_wide`. An example data file is supplied. It can be accessed by 
```r
example_data_path <- system.file(
                "extdata", "sim_data_wide.csv",
                package = "blotIt3"
            )
```
This reads out the provided example file, transferes it to a temporary location and stores the path to this temporary location in `example_data_path`.
The example file is structured as follows
|time|	condition|	ID|	pAKT|	pEPOR|	pJAK2|...|
|--- | --- | --- | --- | ---|--- | ---|
|0	|0Uml Epo	|1	|116.838271399017|	295.836863524109| |...
|5	|0Uml Epo	|1	|138.808500374087|	245.229971713582| |...
|...|...|...|...|...|...|...
0	|0Uml Epo	|2	|94.4670174938645|		|293.604761934545|	...
5	|0Uml Epo	|2	|	|	|398.958892340432|	...
|...|...|...|...|...|...|...

The first three columns contain descriptional data: timepoints, measurement conditions and IDs (e.g. the IDs of the different gels). All following columns contain the measurements of different targets, with the first row containing the names and the following the measurement values corresponding to the time, condition and ID stated in the first columns.

The information which columns contain discriptions has to be passed to `read_wide`:
```r
imported_data <- read_wide(
    file = example_data_path, # path to the example file
    description = seq(1,3), # Indices of columns containing the information
    sep = ",", # sign seperating the colums
    dec = "." # decimal sign
)
```
The result is then a long table of the form

|    |time| condition| ID|  name |    value|
|--- | --- | --- | --- | ---|--- |
pAKT1|       0|  0Uml Epo|  1|  pAKT| 116.83827
pAKT2|       5|  0Uml Epo|  1|  pAKT| 138.80850
pAKT3|      10|  0Uml Epo|  1|  pAKT|  99.09068
pAKT4|      20|  0Uml Epo|  1|  pAKT| 106.68584
pAKT5|      30|  0Uml Epo|  1|  pAKT| 115.02805
pAKT6|      60|  0Uml Epo|  1|  pAKT| 111.91323
pAKT7|     240|  0Uml Epo|  1|  pAKT| 132.56618
|...|...| ...| ...|...|...|

While the first (nameles) colums just contains (unique) row names. New are the columns `name` and `value`. While the column names of the original file are pasted in the former, the latter contains the respective values.
The data.frame `imported_data` can now be passed to the main function.
### Scale data
The full function call is
```r
scaled_data <- align_me(
  data = imported_data,
  model = "yi / sj",
  error_model = "value * sigmaR",
  distinguish = yi ~ name + time + condition,
  scaling = sj ~ name + ID,
  error = sigmaR ~ name + 1,
  input_scale = "linear",
  normalize = TRUE,
  average_techn_rep = FALSE,
  verbose = FALSE,
  normalize_input = TRUE
)
```
We will go now through the parameters individually:
- `data` A long table, usually the output of `read_wide`
- `model` A formula like deskribing the model used for aligning. The present one `yi / sj` means that the measured values `Y_i` are the real values `yi` scaled by scaling factors `sj`. The model therefore is the real value devided by the corresponding scaling factor.
- `error_model` A description of which errors affect the data. Here, only a relative error is present, where the parameter `sigmaR` is scaled by the respective `value`
- `distinguish` Discription of which parameter (left hand side of the tilde) represented by which columns (right hand side of the tilde) contain the "distinguisable effects". In the present example, the model states that the real value is represented by `yi` -- which is the left hand side of the present `distinguish` entry. The present right hand side is "name", "time" and "condition".
In short: we state that the entries "name", "time" and "condition" contain _real_, biological differences.
- `scaling` Same as above, but here is defined which colums contain identificators of different scaling. Here it is "name" and "ID", meaning that measurements with differ in this effects, (but have the same `distinguish` effects) are scaled uppod another.
- `error` Describes how the error affects the values individually. The present formulation means, that the error parameter is _not_ individually adjusted.
- `input_scale` Describes the scale of the present data. `align_me()` accepts "linear", "log", "log2" and "log10". Please be aware the that in the current alpha stage, only "linear" was tested.
- `average_techn_rep` A logical parameter that indicates, if technical replicates should be averaged before the scaling.
- `verbose` If set to `TRUE` additional information will be printed in the console.
- `normalize_input` If set to `TRUE`, the data will be scaled before the actual scaling. This means that the raw input will be scaled to a common order of magnitude before the scaling parameters will be calculated. This is only a computational aid, to eliminate a rare fail of convergence when the different values differ by many orders of magnitude. Setting this to `TRUE` makes only sense (and is only supported) for `input_scale = "linear"`.

The result of `align_me()` is a list with the entries
- `aligned`
