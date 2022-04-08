install.packages("remotes")
remotes::install_github("JimWhiting91/genotype_plot")

library(GenotypePlot)
library(vcfR)

my_vcf <- read.vcfR("SNP-only_simplified_chr7.sweep.recode.vcf.gz")

# popmap = two column data frame with column 1 for individual IDs as they appear in the VCF and column 2 for pop labels
our_popmap <- data.frame(ind = c("PV-AK-1_PV-AK-1","PV-AK-10_PV-AK-10","PV-AK-2_PV-AK-2","PV-AK-3_PV-AK-3","PV-AK-4_PV-AK-4","PV-AK-5_PV-AK-5","PV-AK-6_PV-AK-6","PV-AK-7_PV-AK-7","PV-AK-8_PV-AK-8","PV-AK-9_PV-AK-9","PV-MW-11_PV-MW-11","PV-MW-12_PV-MW-12","PV-MW-13_PV-MW-13","PV-MW-14_PV-MW-14","PV-MW-16_PV-MW-16","PV-MW-2_PV-MW-2","PV-MW-3_PV-MW-3","PV-MW-7_PV-MW-7","PV-MW-9_PV-MW-9","PV-NE-27_PV-NE-27","PV-NE-28_PV-NE-28","PV-NE-29_PV-NE-29","PV-NE-32_PV-NE-32","PV-NE-36_PV-NE-36","PV-NE-37_PV-NE-37","PV-NE-51_PV-NE-51","PV-NR-01_PV-NR-01","PV-NR-02_PV-NR-02","PV-NR-03_PV-NR-03","PV-NR-04_PV-NR-04","PV-NR-05_PV-NR-05","PV-NR-06_PV-NR-06","PV-NR-07_PV-NR-07","PV-NR-08_PV-NR-08","PV-NR-09_PV-NR-09","PV-NR-10_PV-NR-10","PV-NW-12_PV-NW-12","PV-NW-16_PV-NW-16","PV-NW-17_PV-NW-17","PV-NW-18_PV-NW-18","PV-NW-21_PV-NW-21","PV-NW-23_PV-NW-23","PV-NW-25_PV-NW-25","PV-NW-26_PV-NW-26","PV-NW-27_PV-NW-27","PV-NW-8_PV-NW-8","PV-SE-01_PV-SE-01","PV-SE-02_PV-SE-02","PV-SE-03_PV-SE-03","PV-SE-04_PV-SE-04","PV-SE-05_PV-SE-05","PV-SE-10_PV-SE-10","PV-SE-6_PV-SE-6","PV-SE-7_PV-SE-7","PV-SE-8_PV-SE-8","PV-SE-9_PV-SE-9","PV-SR-06_PV-SR-06","PV-SR-12_PV-SR-12","PV-SR-16_PV-SR-16","PV-SR-18_PV-SR-18","PV-SR-19_PV-SR-19","PV-SR-23_PV-SR-23","PV-SR-24_PV-SR-24","PV-SR-27_PV-SR-27","PV-SR-31_PV-SR-31"),
                         pop = c(rep("AK", 10),rep("MW", 9),rep("NE", 7),rep("NR", 10),rep("NW", 10),
                                 rep("SE", 10),rep("SR", 9)),
                         stringsAsFactors = FALSE)

new_plot <- genotype_plot(vcf_object  =  my_vcf,
                          popmap = our_popmap,                              
                          cluster        = FALSE,                           
                          snp_label_size = 5000,                          
                          colour_scheme=c("turquoise2","pink","red"))

pdf("Chr.7_sweep2.pdf")
cowplot::plot_grid(new_plot$positions, new_plot$genotypes, axis="tblr",
                   align="v", nrow=2, ncol=1, rel_heights=c(1,9))
dev.off()
