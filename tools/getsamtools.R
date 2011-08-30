

url <- "http://voxel.dl.sourceforge.net/project/samtools/samtools/0.1.17/samtools-0.1.17.tar.bz2"

outfile <- paste(tempdir(), sep="/", "samtools.tar.bz2")
download.file(url, outfile)

samdir <- paste(tempdir(), sep="/", "samtools")
untar(outfile, compressed="bzip2", exdir=samdir)

samdir2 <- paste(samdir, sep="/", list.files(samdir))
files <- list.files(samdir2, full.names=TRUE)

cfiles <- grep("\\.c$", files, value=TRUE)
hfiles <- grep("\\.h$", files, value=TRUE)

file.copy(c(cfiles, hfiles), "splicing/src/")

dir.create("splicing/src/bcftools")
bcffiles <- list.files(paste(samdir2, sep="/", "bcftools"), full.names=TRUE)

bcf.cfiles <- grep("\\.c$", bcffiles, value=TRUE)
bcf.hfiles <- grep("\\.h$", bcffiles, value=TRUE)

file.copy(bcf.hfiles, "splicing/src/bcftools")
file.copy(bcf.cfiles, "splicing/src/")

unlink(c("splicing/src/bgzip.c", "splicing/src/bamtk.c",
         "splicing/src/main.c", "splicing/src/razip.c"))

