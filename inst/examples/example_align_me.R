## Get example dataset
example_data_file <- system.file(
  "extdata", "sim_data_wide.csv",
  package = "blotIt3"
)

## read in example data by use of 'read_wide()'
example_data <- read_wide(example_data_file, description = seq_len(3))

## execute align_me
out <- align_me(
  data = example_data,
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
