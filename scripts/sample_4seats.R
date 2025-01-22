a <- read.table("sample_4seats.tsv", header = TRUE)
pdf("sample_4seats_boxplots.pdf")
boxplot(ASN ~ Verified, a, main = "4 Seats", horizontal = TRUE, log = "x")
dev.off()
