bio3 <- read.csv("BIOCLIM_3.matched_allvalues.txt")
richness <- read.csv("Astragalus_richness_allvalues.txt")

df <- data.frame(cbind(bio3, richness))
df <- df[df$bio3 != -9999,]
df <- df[df$richness != -9999,]

summary(lm(richness ~ bio3, data = df))
