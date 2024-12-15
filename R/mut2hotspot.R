

mut2hotspot <- function(data){
  # data is a data.table object

  # Initialize a vector to store the indices of rows to be selected
  selected_indices <- integer(0)

  # Iterate over each row in the data table
  for (i in 1:nrow(data)) {
    chrom_i <- data[i, Chromosome]
    start_i <- data[i, Position]
    sampleID_i <- data[i, SampleID]

    # Check for other rows with the same chrom, start difference within 3, and different sampleID
    for (j in 1:nrow(data)) {
      if (i != j) {
        chrom_j <- data[j, Chromosome]
        start_j <- data[j, Position]
        sampleID_j <- data[j, SampleID]
        # the criteria for defining hotspots
        if (chrom_i == chrom_j && abs(start_i - start_j) <= 3 && sampleID_i != sampleID_j) {
          selected_indices <- c(selected_indices, i)
          break
        }
      }
    }
  }

  # Select the rows based on the indices
  selected_rows <- data[selected_indices]

  return(selected_rows)

}
