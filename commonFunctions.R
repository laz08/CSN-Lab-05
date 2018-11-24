

saveNodesDegreeOnFile <- function(t, k, baseName) {
    
    filename = paste(baseName, sprintf("%05d", t), ".csv", sep = "")
    cat("Saving degree sequence", filename, "\n")
    write.csv(k, file=filename, row.names = FALSE)
    cat("Done.\n")
}
