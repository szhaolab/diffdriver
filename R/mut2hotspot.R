

mut2hotspot <- function(data) {
  # Ensure data is a data.table
  data.table::setDT(data)

  # Create a shifted version of the data to compare positions
  data_shifted <- copy(data)
  data_shifted[, Position_shifted_right := Position + 3]
  data_shifted[, Position_shifted_left := Position -3]

  # Perform a non-equi join to find matching rows
  result <- data[data_shifted, on = .(Chromosome, Position >= Position_shifted_left, Position <= Position_shifted_right),
                 nomatch = 0, allow.cartesian = TRUE]

  # Filter out rows with the same SampleID
  result <- result[SampleID != i.SampleID]

  # Select unique rows based on the original data
  selected_rows <- unique(result[, .(Chromosome, Position)])

  return(selected_rows)
}
