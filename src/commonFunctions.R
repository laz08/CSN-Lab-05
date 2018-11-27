

saveNodesDegreeOnFile <- function(t, k, baseName) {
    
    filename = paste(baseName, sprintf("%05d", t), ".csv", sep = "")
    cat("Saving degree sequence", filename, "\n")
    write.csv(k, file=filename, row.names = FALSE)
    cat("Done.\n")
}

saveNodesSequenceOnFile <- function(t, k, baseName) {
    
    filename = paste(baseName, sprintf("%05d", t), ".csv", sep = "")
    cat("Saving degree sequence", filename, "\n")
    write.csv(k, file=filename, row.names = FALSE)
    cat("Done.\n")
}


plotVertexEvolution <- function(table, x, vid) {
    
    plot <- ggplot(data = table) +
        aes(x = sequ, y = table[, x]) +
        geom_point(color = '#0c4c8a') +
        labs(title = paste('Vertex',vid, 'degree evolution'),
             x = 't',
             y = 'k') +
        theme_calc()
    
    return(plot)
}


plotAllEvolutions <- function(table, title){
    p1 <- plotVertexEvolution(table, 2, 1)
    p2 <- plotVertexEvolution(table, 3, 10)
    p3 <- plotVertexEvolution(table, 4, 100)
    p4 <- plotVertexEvolution(table, 5, 1000)
    grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2, top=title)
}


loadModelExecutions <- function(files) {
    
    ks.1 = read.csv2(files[1]); ks.1 = ks.1$x
    ks.2 = read.csv2(files[2]); ks.2 = ks.2$x
    ks.3 = read.csv2(files[3]); ks.3 = ks.3$x
    ks.4 = read.csv2(files[4]); ks.4 = ks.4$x
    
    table = data.frame(sequ = seq(length(ks.1)), ks.1, ks.2, ks.3, ks.4)
    return(table)
}
