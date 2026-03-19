library(MSnbase)
library(dplyr)

# Let's load the test dataset of MSnbase
data(itraqdata)

######################################
# How many spectra do we have here?
length(itraqdata)
# [1] 55

# Which are those spectra?
featureNames(itraqdata)
#  [1] "X1"  "X10" "X11" "X12" "X13" "X14" "X15" "X16" "X17" "X18" "X19" "X2"
# [13] "X20" "X21" "X22" "X23" "X24" "X25" "X26" "X27" "X28" "X29" "X3"  "X30"
# [25] "X31" "X32" "X33" "X34" "X35" "X36" "X37" "X38" "X39" "X4"  "X40" "X41"
# [37] "X42" "X43" "X44" "X45" "X46" "X47" "X48" "X49" "X5"  "X50" "X51" "X52"
# [49] "X53" "X54" "X55" "X6"  "X7"  "X8"  "X9"

# Let's load spectrum X42
sp = itraqdata[["X42"]]
sp
# Object of class "Spectrum2"
#  Precursor: 803.9468
#  Retention time: 33:38
#  Charge: 2
#  MSn level: 2
#  Peaks count: 2001
#  Total ion count: 49352468

# Let's have a look (we need to have XMing installed for X11 forwarding)
plot(sp)

# Let's zoom into a region
plot(trimMz(sp, mzlim = c(450, 500)))

# Can we see the reporter ions?
plot(trimMz(sp, mzlim = c(110, 130)))

######################################
# Since this has been identified already, we can check the identification data
fData(itraqdata)["X42", c("spectrum", "ProteinAccession", "ProteinDescription", "PeptideSequence")]
#     spectrum ProteinAccession             ProteinDescription PeptideSequence
# X42       42          ECA0631 conserved hypothetical protein    LTDDDLTVIEGK

# Let's confirm whether we get the same identification with Mascot
# First we create a fuction to filter a bit the lot intensity peaks

showMgfData <- function(sp, t = 10000) {
  tf <- tempfile(fileext = ".mgf")
  on.exit(unlink(tf), add = TRUE)

  sp2 <- pickPeaks(sp)
  sp2 <- clean(removePeaks(sp2, t = t))

  writeMgfData(sp2, con = tf)
  cat(readLines(tf), sep = "\n")
  invisible(NULL)
}

# Then, we obtain the MGF text
showMgfData(sp)

# But this practical is about quantification
# So, let's quantify these spectra!
msnset <- quantify(itraqdata, reporters = iTRAQ4, method = "trap")
View(exprs(msnset))
msnset_log <- log(msnset, base = 2)

itraq <- as.data.frame(exprs(msnset))

# Let's check the channel distributions
boxplot(
  log2(itraq + 1),
  main = "Reporter intensities by iTRAQ channel",
  ylab = "log2 intensity"
)

# Let's have a look at X47, which is enolase, a spike in protein in this sample
sp = itraqdata[["X47"]]
plot(sp)
plot(trimMz(sp, mzlim = c(110, 130)))

######################################
# Let's normalise the dataset
msnset_n <- normalize(msnset_log, method = "center.median")

# It would be nice to have the identification data on the same table
itraq_table <- cbind(
  ProteinAccession = fData(msnset_n)$ProteinAccession,
  PeptideSequence  = fData(msnset_n)$PeptideSequence,
  as.data.frame(exprs(msnset_n))
)

View(itraq_table)

######################################
# We filter out some bad data
ok <- with(itraq_table,
           !is.na(iTRAQ4.114) & !is.na(iTRAQ4.115))

itraq_ok <- itraq_table[ok, ]

# And generate the MA plot of 115 vs 114
itraq_ok$logFC_115_114 <- itraq_ok$iTRAQ4.115 - itraq_ok$iTRAQ4.114
itraq_ok$A_115_114 <- (itraq_ok$iTRAQ4.115 + itraq_ok$iTRAQ4.114) / 2

plot(itraq_ok$A_115_114, itraq_ok$logFC_115_114,
     pch = 16, cex = 0.8,
     xlab = "A = average normalized log2-intensity (114,115)",
     ylab = "M = normalized log2-difference (115 - 114)",
     main = "iTRAQ MA plot: 115 vs 114")
abline(h = 0, col = "red")

# Labelling the MA plot
text(itraq_ok$A_115_114, itraq_ok$logFC_115_114,
     labels = itraq_ok$ProteinAccession,
     pos = 3, cex = 0.6)

# Repeat MA plot, but labelling with the protein ids
plot(itraq_ok$A_115_114, itraq_ok$logFC_115_114,
     pch = 16, cex = 0.8,
     xlab = "A = average normalized log2-intensity (114,115)",
     ylab = "M = normalized log2-difference (115 - 114)",
     main = "iTRAQ MA plot: 115 vs 114")
abline(h = 0, col = "red")
text(itraq_ok$A_115_114, itraq_ok$logFC_115_114,
     labels = itraq_ok$ProteinAccession,
     pos = 3, cex = 0.6)

######################################
# We will try to make a volcano plot at protein level
# In this case we will accept only one spectrum
# because this is a test dataset and we have almost nothing
prot_tab <- itraq_ok |>
  mutate(logratio = iTRAQ4.115 - iTRAQ4.114) |>
  group_by(ProteinAccession) |>
  summarize(
    n_spectra = n(),
    logFC = mean(logratio, na.rm = TRUE),
    pval = if (n() >= 2) t.test(logratio, mu = 0)$p.value else 1,
    .groups = "drop"
  ) |>
  mutate(neglog10p = -log10(pval))

sel <- prot_tab$n_spectra >= 2

plot(prot_tab$logFC, prot_tab$neglog10p,
     pch = 16, cex = 0.8,
     xlab = "Protein mean normalized log2-difference (115 - 114)",
     ylab = "-log10 p-value",
     main = "Protein-level volcano: 115 vs 114")
abline(v = 0, col = "red")

text(prot_tab$logFC[sel], prot_tab$neglog10p[sel],
     labels = prot_tab$ProteinAccession[sel],
     pos = 3, cex = 0.8)