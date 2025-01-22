a <- read.table("sample_3seats.tsv", header = TRUE)
pdf("sample_3seats_boxplots.pdf")
boxplot(ASN ~ Verified, a, main = "3 Seats", horizontal = TRUE, log = "x")
dev.off()
