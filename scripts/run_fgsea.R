#!/usr/bin/env Rscript

library(argparse)
library(fgsea)
library(data.table)
library(ggplot2)

parser <- ArgumentParser(description="Run GSEA with fgsea on a ranked gene list from DESeq2.")

parser$add_argument("--deseq2_file", type="character", required=TRUE,
                    help="Path to the DESeq2 result TSV file.")
parser$add_argument("--geneset_file", type="character", required=TRUE,
                    help="Path to the gene set file in .gmt format.")
parser$add_argument("--outdir", type="character", required=TRUE,
                    help="Directory to save the GSEA results.")
parser$add_argument("--top_n", type="integer", default=20,
                    help="Number of top pathways to plot in the table.")

args <- parser$parse_args()

# --- 1. Load Data ---
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

de_results <- fread(args$deseq2_file)
pathways <- gmtPathways(args$geneset_file)

# --- 2. Create Ranked Gene List ---
# Use the 'stat' column for ranking, as it's a robust measure
# Remove NA values and sort in descending order
ranks <- de_results[!is.na(stat), .(gene, stat)]
setkey(ranks, gene)
ranks_vec <- ranks$stat
names(ranks_vec) <- ranks$gene

# --- 3. Run fgsea ---
message("Running fgsea...")
fgsea_res <- fgsea(pathways = pathways, 
                   stats    = ranks_vec,
                   minSize  = 15,
                   maxSize  = 500)

fwrite(fgsea_res, file.path(args$outdir, "fgsea_results.tsv"), sep = "	", sep2 = c("", " ", ""))

# --- 4. Generate Plots ---
message("Generating result plots...")

# Plot the top N pathways
top_pathways_up <- fgsea_res[ES > 0][head(order(pval), n=args$top_n), pathway]
top_pathways_down <- fgsea_res[ES < 0][head(order(pval), n=args$top_n), pathway]
top_pathways <- c(top_pathways_up, rev(top_pathways_down))

pdf(file.path(args$outdir, "top_pathways_table.pdf"), width=12, height=max(6, length(top_pathways) * 0.3))
plotGseaTable(pathways[top_pathways], ranks_vec, fgsea_res, gseaParam=0.5)
dev.off()

# Plot enrichment for a top pathway (if any are significant)
if (nrow(fgsea_res[padj < 0.05]) > 0) {
    top_pathway_name <- fgsea_res[order(padj), pathway[1]]
    
    pdf(file.path(args$outdir, "top_enrichment_plot.pdf"), width=6, height=4)
    plotEnrichment(pathways[[top_pathway_name]], ranks_vec) + 
        labs(title=top_pathway_name)
    dev.off()
}

message("GSEA analysis complete.")



