## Modified version of RNA-seq analysis with DESeq2 by Stephen Turner, @genetics_blog

# RNA-seq data from GSE52202
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse52202. All patients with
# ALS, 4 with C9 expansion ('exp'), 4 controls without expansion ('ctl')

# Import & pre-process ----------------------------------------------------

# Import data from featureCounts Previously ran at command line something like
# this: featureCounts -a genes.gtf -o counts.txt -T 12 -t exon -g gene_id

deseq2pipeline <- function(x, output_path, num_ctl, num_exp) {
    # dirn <- dirname(x)
    filen <- basename(x)

    dispersion_plot_file <- paste(output_path, substr(filen, 1, nchar(filen) - 4),
        "_dispersion.png", sep = "")
    ma_plot_file <- paste(output_path, substr(filen, 1, nchar(filen) - 4), "_ma.png",
        sep = "")
    volcano_plot_file <- paste(output_path, substr(filen, 1, nchar(filen) - 4), "_volcano.png",
        sep = "")
    heatmap_plot_file <- paste(output_path, substr(filen, 1, nchar(filen) - 4), "_sample_heatmap.png",
        sep = "")
    pca_plot_file <- paste(output_path, substr(filen, 1, nchar(filen) - 4), "_pca.png",
        sep = "")
    result_file <- paste(output_path, substr(filen, 1, nchar(filen) - 4), "_result.csv",
        sep = "")


    countdata <- read.csv(x, row.names = 1)

    # Remove first five columns (chr, start, end, strand, length)
    countdata <- countdata[, 6:ncol(countdata)]
    for (i in seq_len(ncol(countdata))) {
        # for-loop over columns
        countdata[, i] <- as.integer(countdata[, i])
    }

    countdata <- as.matrix(countdata)
    print(head(countdata))

    # if (length(colnames(countdata)) == (2 * num_replicates)) {
    (condition <- factor(c(rep("ctl", num_ctl), rep("exp", num_exp))))

    # Analysis with DESeq2 -----------------------------------------
    library(DESeq2)

    # Create a coldata frame and instantiate the DESeqDataSet.  See
    # ?DESeqDataSetFromMatrix
    (coldata <- data.frame(row.names = colnames(countdata), condition))
    dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
    dds

    # Run the DESeq pipeline
    dds <- DESeq(dds)

    # Plot dispersions
    png(dispersion_plot_file, 1000, 1000, pointsize = 20)
    plotDispEsts(dds, main = "Dispersion plot")
    dev.off()

    # Regularized log transformation for clustering/heatmaps, etc
    rld <- rlogTransformation(dds)
    head(assay(rld))
    hist(assay(rld))

    # Colors for plots below Ugly: Use RColorBrewer, better
    library(RColorBrewer)
    (mycols <- brewer.pal(8, "Dark2")[seq_len(length(unique(condition)))])

    # Sample distance heatmap
    sample_dists <- as.matrix(dist(t(assay(rld))))
    library(gplots)
    png(heatmap_plot_file, w = 1000, h = 1000, pointsize = 20)
    heatmap.2(as.matrix(sample_dists), key = F, trace = "none", col = colorpanel(100,
                                                                                 "black", "white"), ColSideColors = mycols[condition], RowSideColors = mycols[condition],
              margin = c(10, 10), main = "Sample Distance Matrix")
    dev.off()

    png(pca_plot_file, 1000, 1000, pointsize = 20)
    DESeq2::plotPCA(rld, intgroup = "condition")
    dev.off()

    # Get differential expression results
    res <- results(dds)
    table(res$padj < 0.05)
    ## Order by adjusted p-value
    res <- res[order(res$padj), ]
    ## Merge with normalized count data
    resdata <- merge(as.data.frame(res),
                     as.data.frame(counts(dds, normalized = TRUE)),
                     by = "row.names", sort = FALSE)
    names(resdata)[1] <- "Gene"
    head(resdata)
    ## Write results
    write.csv(resdata, file = result_file)

    ## Examine plot of p-values hist(res$pvalue, breaks=50, col='grey')

    ## Examine independent filtering attr(res, 'filterThreshold')
    ## plot(attr(res,'filterNumRej'), type='b', xlab='quantiles of
    ## baseMean', ylab='number of rejections')

    png(ma_plot_file, 1500, 1000, pointsize = 20)
    DESeq2::plotMA(dds, ylim = c(-1, 1))
    dev.off()

    ## Volcano plot with 'significant' genes labeled
    volcanoplot <- function(res, lfcthresh = 2, sigthresh = 0.05, main = "Volcano Plot",
                            legendpos = "bottomright", labelsig = TRUE, textcx = 1, ...) {
        with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = main,
                       ...))
        with(subset(res, padj < sigthresh), points(log2FoldChange, -log10(pvalue),
                                                   pch = 20, col = "red", ...))
        with(subset(res, abs(log2FoldChange) > lfcthresh), points(log2FoldChange,
                                                                  -log10(pvalue), pch = 20, col = "orange", ...))
        with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh),
             points(log2FoldChange, -log10(pvalue), pch = 20, col = "green", ...))
        if (labelsig) {
            require(calibrate)
            with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh),
                 textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = textcx,
                        ...))
        }
        legend(legendpos, xjust = 1, yjust = 1, legend = c(paste("FDR<", sigthresh,
                                                                 sep = ""), paste("|LogFC|>", lfcthresh, sep = ""), "both"), pch = 20,
               col = c("red", "orange", "green"))
    }
    png(volcano_plot_file, 1200, 1000, pointsize = 20)

    volcanoplot(resdata, lfcthresh = 1, sigthresh = 0.05, textcx = 0.8,
                xlim = c(floor(min(resdata$log2FoldChange, na.rm = TRUE)),
                         ceiling(max(resdata$log2FoldChange, na.rm = TRUE))),
                ylim = c(0,
                         max(-log10(resdata$pvalue), na.rm = TRUE) + 0.5)
    )
    dev.off()
}
