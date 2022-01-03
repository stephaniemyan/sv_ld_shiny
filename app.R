library(shiny)
library(rsconnect)
library(dplyr)
library(data.table)

#================================================================#
# Load in LD data. Pre-processing steps: ======================= #
# 1. Subset to just LD between SNPs and SVs ==================== #
# 2. Remove SVs that violate HWE in more than 50% of populations #
# 3. Remove SVs with genotyping rate < 50% ===================== #
# 4. Remove SVs with AF = 0 ==================================== #
# 5. Annotate with SV length =================================== #
# 6. Annotate SNPs with rsIDs from dbSNP ======================= #
#================================================================#

# read in and filter LD data
# all_pops_ld <- fread("all_pops_ld_20210114.txt")
# gt_callrate <- fread("genotyping_callrate.txt")
# hwe_pvals <- fread("HWE_pvals.txt")
# afs <- fread("eichlerSVs_af_allpops_unfolded.frq.strat.gz") %>%
#   group_by(SNP) %>%
#   summarize(max_af = max(MAF)) %>%
#   setDT()
# ld <- semi_join(all_pops_ld,
#                 gt_callrate[ call_rate >= 0.5 ],
#                 by = "SV") %>%
#   semi_join(., hwe_pvals[ non_hwe_counts <= 13 ],
#             by = "SV") %>%
#   semi_join(., afs[ max_af > 0 ],
#             by = c("SV" = "SNP")) %>%
#   as.data.table()

# # change column names
# setnames(ld, c("sv_id", "sv_chr", "sv_pos", "snp_rsid", "snp_chr", "snp_pos", "r2", "pop"))
# # rearrange column order
# ld <- ld[, c("sv_id", "sv_chr", "sv_pos", "snp_rsid", "snp_chr", "snp_pos", "r2", "pop")]
# # modify chromosome columns
# ld[, sv_chr := paste0("chr", sv_chr)]
# ld[, snp_chr := paste0("chr", snp_chr)]
# ld[snp_chr == "chr23", snp_chr := "chrX"]
# ld[sv_chr == "chr23", sv_chr := "chrX"]
# # add superpopulation annotations
# ld[pop %in% c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"), superpop := "AFR"]
# ld[pop %in% c("CLM", "MXL", "PEL", "PUR"), superpop := "AMR"]
# ld[pop %in% c("CDX", "CHB", "CHS", "JPT", "KHV"), superpop := "EAS"]
# ld[pop %in% c("CEU", "FIN", "GBR", "IBS", "TSI"), superpop := "EUR"]
# ld[pop %in% c("BEB", "GIH", "ITU", "PJL", "STU"), superpop := "SAS"]

# add in SV length data
# svlen <- fread("paragraph_sv_lengths.vcf",
#                col.names = c("CHROM", "POS", "ID", "sv_length"))
# ld <- merge(ld, svlen[, c("ID", "sv_length")],
#             by.x = c("sv_id"), by.y = c("ID"))

# add in rsIDs from dbsnp - this doesn't account for duplicate rsIDs at the same position
# dbsnp_path <- "/scratch4/mschatz1/rmccoy22/dbSNP/"
# get_rsids <- function(chrom) {
#   subset <- ld[snp_chr == chrom]
#   vcf <- fread(cmd = paste0("zcat ", dbsnp_path,
#                             "GCF_000001405.38.", chrom, ".corrected.vcf.gz",
#                             " | grep -v '##'"))
#   merged <- merge(ld, vcf[, c("#CHROM", "POS", "ID", "REF", "ALT")],
#                   by.x = c("snp_chr", "snp_pos"), by.y = c("#CHROM", "POS"),
#                   all.x = TRUE)
# }
# ld_rsid <- pbmclapply(as.list(c(1:22, "X")),
#                       function(x) get_rsids(paste0("chr", x)))
# ld <- ld_rsid[, c("sv_id", "sv_chr", "sv_pos",  "sv_length", "ID", "snp_chr",
#                   "snp_pos", "REF", "ALT", "r2", "pop", "superpop")]
# setnames(ld, c("ID", "REF", "ALT"), c("snp_rsid", "snp_ref", "snp_alt"))

# fwrite(ld, "ld_filtered_20210114.txt", sep = "\t")
ld <- fread(cmd = "zcat ld_filtered_20210114.txt")

#==========================#
# Customize user interface #
#==========================#

ui <- fluidPage(
  
  # title
  titlePanel("Structural variant LD browser"),
  
  # table subsetting options
  sidebarPanel(
    # data description
    p("LD within each 1000 Genomes population was calculated between SVs and SNPs/indels within a 100-variant and 100 Mb window, as described in",
      a("Yan et al.",
        href = "https://elifesciences.org/articles/67615#s4"),
      "Based on position, 1000 Genomes SNPs/indels were assigned rsIDs from",
      a("dbSNP build 155.",
        href = "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/")),
    # tags$div(tags$ul(
      # tags$li("LD within each 1000 Genomes population was calculated between SVs and SNPs/indels within a 100-variant and 100 Mb window, as described in",
      #         a("Yan et al.",
      #           href = "https://elifesciences.org/articles/67615#s4")),
      # tags$li("Based on position, 1000 Genomes SNPs/indels were assigned rsIDs from",
      #         a("dbSNP build 155",
      #           href = "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/")))),
    
    # responsive input of rsID or SV ID
    selectizeInput("snp_rsid", label = "SNP rsID",
                   choices = NULL, selected = NULL),
    selectizeInput("sv_id", label = "SV ID",
                   choices = NULL, selected = NULL),
    # other inputs
    selectInput("sv_chr", "Chromosome:",
                c("All", unique(as.character(ld$sv_chr)))),
    textInput("snp_pos", label = "SNP position", placeholder = "796338"),
    textInput("sv_pos", label = "SV position", placeholder = "814625"),
    selectInput("pop", "1000 Genomes population:",
                c("All", unique(as.character(ld$pop)))),
    selectInput("superpop", "1000 Genomes superpopulation:",
                c("All", unique(as.character(ld$superpop)))),
    # manual input of r2 threshold
    textInput("r2", label = c(list(HTML("r<sup>2</sup>")), " threshold"), placeholder = "0.5")
  ),
  
  # app description
  mainPanel(
    p("Find SNPs in linkage disequilibrium (LD) with structural variants (SVs), and vice versa, on GRCh38."),
    p("SVs were called from long-read sequencing data of",
      a("15 individuals",
        href = "https://doi.org/10.1016/j.cell.2018.12.019"),
      "and genotyped in",
      a("2,504 individuals",
        href = "https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2"),
      "from the 1000 Genomes Project. SV calling and genotyping, as well as LD analysis, was performed by",
      a("Yan et al. (2021).",
        href = "https://elifesciences.org/articles/67615")),
    p("If you have questions or comments, contact Stephanie Yan at ",
      a("syan@jhu.edu.",
        href = "mailto:syan@jhu.edu"),
      "Source code for the Shiny app is available",
      a("here.",
        href = "https://github.com/stephaniemyan/sv_ld_shiny")),
    br(),
    
    # display table output
    DT::dataTableOutput("table")
  ),
)

#=============================#
# Define server logic for app #
#=============================#

# Define server logic----
server <- function(input, output, session) {
  
  # update input selection for selectize inputs
  updateSelectizeInput(session, 'snp_rsid', choices = c("All", ld$snp_rsid), server = TRUE)
  updateSelectizeInput(session, 'sv_id', choices = c("All", ld$sv_id), server = TRUE)
  updateSelectizeInput(session, 'snp_pos', choices = c("All", ld$snp_pos), server = TRUE)
  
  # filter data based on selections
  output$table <- DT::renderDataTable(DT::datatable({
    if (input$snp_rsid != "All") {
      ld <- ld[ld$snp_rsid == input$snp_rsid,]
    }
    if (input$sv_id != "All") {
      ld <- ld[ld$sv_id == input$sv_id,]
    }
    if (input$sv_chr != "All") {
      ld <- ld[ld$sv_chr == input$sv_chr,]
    }
    if (input$snp_pos != "") {
      ld <- ld[ld$snp_pos == input$snp_pos,]
    }
    if (input$sv_pos != "") {
      ld <- ld[ld$sv_pos == input$sv_pos,]
    }
    if (input$pop != "All") {
      ld <- ld[ld$pop == input$pop,]
    }
    if (input$superpop != "All") {
      ld <- ld[ld$superpop == input$superpop,]
    }
    if (input$r2 != "") {
      ld <- ld[ld$r2 >= as.numeric(input$r2),]
    }
    ld
  }))
  
}

#===============#
# Call shinyApp #
#===============#

shinyApp(ui = ui, server = server)
