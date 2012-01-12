
library(splicing)

iso3file <- system.file("test_files/iso3.gff", package="splicing")

iso <- readGFF3(iso3file)

iso1 <- selectIso(iso, 1)
iso2 <- selectIso(iso, 2)
iso3 <- selectIso(iso, 3)

iso12 <- selectIso(iso, 1:2)
iso13 <- selectIso(iso, c(1,3))
iso23 <- selectIso(iso, 2:3)

str(iso1)
str(iso2)
str(iso3)

str(iso12)
str(iso13)
str(iso23)
