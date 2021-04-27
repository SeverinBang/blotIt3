pAKT5|      30|  0Uml Epo|  1|  pAKT| 115.02805
pAKT6|      60|  0Uml Epo|  1|  pAKT| 111.91323
pAKT7|     240|  0Uml Epo|  1|  pAKT| 132.56618
  ...|...| ...| ...|...|...|

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
