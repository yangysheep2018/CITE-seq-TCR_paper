library(ineq)
library(immunarch)  # Load the package into R
library(stringr)
# data(immdata)  # Load the test dataset
# immdata$data <- top(immdata$data, 4000)

file_path = "/Users/yanyang/Desktop/Bomi/TCR/TCR_outs/outs"
immdata_10x <- repLoad(file_path)
df = as.data.frame(repDiversity(.data = immdata_10x$data, .method = "gini", .do.norm = NA, .laplace = 0))
df$V2 = str_extract(rownames(df), "[^_]+")
df$V3 = str_extract(df$V2, "\\D*(?=\\d)")

ggplot(df, aes(x=V3, y=V1,fill = V3)) + 
  geom_boxplot() +
  scale_fill_brewer(palette="Set1") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.background = element_rect(colour = NA))

immdata_10x$meta$s = str_extract(immdata_10x$meta$Sample, "[^_]+")
immdata_10x$meta$g = str_extract(immdata_10x$meta$s, "\\D*(?=\\d)")
data(immdata)  # Load the test dataset
immdata$data <- top(immdata$data, 4000)
imm_raref <- repDiversity(immdata$data, "raref", .verbose = F)
p1 <- vis(imm_raref,.by = "Status", .meta = immdata$meta)
immdata$meta
immdata_10x$meta


imm_raref <- repDiversity(immdata_10x$data, "raref", .verbose = F)

p1 <- vis(imm_raref,.by = "s", .meta = immdata_10x$meta)
