' shuffled_joined_to_fastq.R

Usage:
    shuffled_joined_to_fastq.R --R1_fpath R1_FPATH --R2_fpath R2_FPATH --output_dir OUTPUT_DIR [ --sample_index SAMPLE_INDEX ]

Options:
    -r --R1_fpath R1_FPATH          Full path to shuffled and joined file for R1
    -s --R2_fpath R2_FPATH          Full path to shuffled and joined file for R2
    -o --output_dir OUTPUT_DIR      Output directory
    -d --sample_index SAMPLE_INDEX  Sample index in read header (e.g. AGTCAGTC:AGTCAGTC)
' -> doc

library(docopt)
args <- docopt(doc)

suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(GenomicFeatures))
suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg38))
# genome: hg38
# provider: UCSC
# release date: Feb 2019

# args = list()
## F19K16+F24B22 only
# TMP_dir = "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/make_toy_fastqs_outputs/outputs/output_allchrs/toy_sample2/TMP/"
# args$R1_fpath = paste0(TMP_dir, "CMP-01-02-cfDNA-03.all.shuf222.R1.joined.F19K16_F24B22")
# args$R2_fpath = paste0(TMP_dir, "CMP-01-02-cfDNA-03.all.shuf222.R2.joined.F19K16_F24B22")
# args$R1_fpath = paste0(TMP_dir, "CMP-01-02-cfDNA-04.all.shuf222.R1.joined.F19K16_F24B22")
# args$R2_fpath = paste0(TMP_dir, "CMP-01-02-cfDNA-04.all.shuf222.R2.joined.F19K16_F24B22")
# args$output_dir = TMP_dir
# args$sample_index = "EFGH:EFGH"

## hg38+F19K16+F24B22
# TMP_dir = "/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/make_toy_fastqs_outputs/outputs/output_allchrs/TMP/"
# args$R1_fpath = paste0(TMP_dir, "CMP-01-02-cfDNA-03.all.shuf111.R1.joined")
# args$R2_fpath = paste0(TMP_dir, "CMP-01-02-cfDNA-03.all.shuf111.R2.joined")
# args$output_dir = TMP_dir
# args$sample_index = "AAAA:AAAA"

R1.joined.fpath = args$R1_fpath
R2.joined.fpath = args$R2_fpath
out.dir = paste0(args$output_dir,"/")

#if (!is.null(args$fake_fragname_seed)){
#  fragname_seed = args$fake_fragname_seed
#} else {
#  fragname_seed = 42
#}

## should be same for all samples subsampled from into new test sample
if (!is.null(args$sample_index)){
  sample_index = args$sample_index
} else {
  sample_index = "AGTCAGTC:AGTCAGTC"
}
# sample_index

R1.out.fname = paste0(basename(R1.joined.fpath) %>% gsub("joined","hg38.",.), gsub(":","",sample_index), ".fastq")
# R1.out.fname
R2.out.fname = paste0(basename(R2.joined.fpath) %>% gsub("joined","hg38.",.), gsub(":","",sample_index), ".fastq")
# R2.out.fname

## R1 #########################################################################

# read in data
R1.joined = read.table(R1.joined.fpath, sep = '\t', header = F, comment.char="")
# R1.joined %>% head()

R1.joined.BACs =
  R1.joined %>%
  filter(V6 == "F19K16" | V6 == "F24B22")
# R1.joined.BACs$V6 %>% unique()
# R1.joined.BACs %>% head()

R1.joined.hg38 =
  R1.joined %>%
  filter(V6 != "F19K16" & V6 != "F24B22")
# R1.joined.hg38 %>% head()

## extend all to 75bp
## ( ideally in both directions, but for now just change 'end'='start'+75,
## but have to take into account 3nt UMI + linker )
read_len = R1.joined.hg38$V3 %>% nchar() %>% unique()
R1.newEnd =
  R1.joined.hg38 %>%
  mutate(end = if_else(nchar(V9) == 3, V7 + read_len - (3+1) - 1,
                       if_else(nchar(V9) == 4, V7 + read_len - (4+1) - 1, 0)))
# R1.newEnd %>% nrow()

## dataframe to granges so can use DNAStringSet()
if (nrow(R1.newEnd) != 0){
  R1.gr =
    makeGRangesFromDataFrame(R1.newEnd,
                             seqnames.field = "V6",
                             start.field = "V7",
                             end.field = "end",
                             keep.extra.columns = T)
  R1.gr$hg38seq = DNAStringSet(as.character(Views(Hsapiens, R1.gr)))
  # R1.gr

  ## back to dataframe
  R1.hg38seq =
    data.frame(chr = R1.gr %>% seqnames(),
               start = R1.gr %>% start(),
               end = R1.gr %>% end(),
               umi = R1.gr$V9,
               sequence = R1.gr$hg38seq,
               spacer = R1.gr$V4,
               phred = R1.gr$V5,
               strand = R1.gr$V11) %>%
    mutate(sequence_umi = if_else(R1.gr$V11 == "+", paste0(umi, "T", R1.gr$hg38seq),
                                  if_else(R1.gr$V11 == "-", paste0(umi, "T", reverseComplement(R1.gr$hg38seq)),
                                          "ERROR"))) ## { if -ve strand in bam, have to rev_comp to turn it into FASTQ sequence }
  # R1.hg38seq$sequence_umi %>% nchar() %>% unique() # [1] 75
  # R1.hg38seq %>% head()
}

# R1.joined.BACs %>% head()
R1.joined.BACs.rearranged =
  data.frame(chr = R1.joined.BACs$V6,
             start = R1.joined.BACs$V7,
             end = R1.joined.BACs$V8,
             umi = R1.joined.BACs$V9,
             sequence_umi = R1.joined.BACs$V3,
             spacer = "+",
             phred = R1.joined.BACs$V5,
             strand = R1.joined.BACs$V11)
# R1.joined.BACs.rearranged %>% head()

if (nrow(R1.newEnd) != 0 & nrow(R1.joined.BACs.rearranged) != 0){
  R1.hg38seq.BAC.rbound = bind_rows(R1.hg38seq, R1.joined.BACs.rearranged)
} else {
  R1.hg38seq.BAC.rbound = R1.joined.BACs.rearranged
}
# R1.hg38seq.BAC.rbound$chr %>% unique()
# R1.hg38seq.BAC.rbound %>% head()
# R1.hg38seq.BAC.rbound %>% filter(chr == "F19K16" | chr == "F24B22") %>% head()

## make dummy frag names
## @TEST123:123:TESTING:1:random_4to5:random_4to5:random_4to5 1:N:0:random_8+random_8
vec_len = nrow(R1.hg38seq.BAC.rbound)
message(paste0("Number of reads subsampled: ", vec_len))
instrument_flowcell = rep("@TEST123:123:TESTING:1", vec_len)
R1_sample_index = rep(paste0("1:N:0:", sample_index), vec_len)

fragname_seed = sample(1:10000, 1);
message(paste0("Seed for these R1/R2 file's fragment names is: ", fragname_seed))
set.seed(fragname_seed)
flowcell_coord = paste(sample(1000:99999, vec_len, replace = TRUE),
                       sample(1000:99999, vec_len, replace = TRUE),
                       sample(1000:99999, vec_len, replace = TRUE),
                       sep = ":")
R1_frag_names = paste0(instrument_flowcell, ":", flowcell_coord, " ", R1_sample_index)
# R1_frag_names %>% head()

## add frag names to R1
R1.hg38seq.frag_name =
  R1.hg38seq.BAC.rbound %>%
  dplyr::mutate(frag_names = R1_frag_names) %>%
  dplyr::select(frag_names, sequence_umi, spacer, phred)
# R1.hg38seq.frag_name %>% head()
# R1.hg38seq.frag_name %>% filter(sequence_umi == "GTCTCATGCTCCCAACCAGCGTCTCTTCCTCTGCTTCCAATGGCGTGTTCCTCTTGTTCCCCGCTCGATCTGCCT")

R1.hg38seq.fastq =
  R1.hg38seq.frag_name %>%
  pivot_longer(cols = everything(), names_to = "col1", values_to = "col2") %>%
  dplyr::select(col2)
# R1.hg38seq.fastq %>% head(12)

message(paste0("Writing fastq for ", R1.out.fname))
write.table(R1.hg38seq.fastq,
            file = paste0(out.dir, R1.out.fname),
            quote = F, sep = "", row.names = F, col.names = F)


## R2 #########################################################################

# read in data
R2.joined = read.table(R2.joined.fpath, sep = '\t', header = F, comment.char="")
# R2.joined %>% head()

R2.joined.BACs =
  R2.joined %>%
  filter(V6 == "F19K16" | V6 == "F24B22")
# R2.joined.BACs$V6 %>% unique()
# R2.joined.BACs %>% head()

R2.joined.hg38 =
  R2.joined %>%
  filter(V6 != "F19K16" & V6 != "F24B22")
# R2.joined.hg38 %>% head()

## extend all to 75bp
## ( ideally in both directions, but for now just change 'end'='start'+75,
## but have to take into account 3nt UMI + linker )
read_len = R2.joined.hg38$V3 %>% nchar() %>% unique()
R2.newEnd =
  R2.joined.hg38 %>%
  mutate(end = if_else(nchar(V10) == 3, V7 + read_len - (3+1) - 1,                 ## R2 is V10, uses the other UMI
                       if_else(nchar(V10) == 4, V7 + read_len - (4+1) - 1, 0)))

## dataframe to granges so can use DNAStringSet()
if (nrow(R2.newEnd) != 0){
  R2.gr =
    makeGRangesFromDataFrame(R2.newEnd,
                             seqnames.field = "V6",
                             start.field = "V7",
                             end.field = "end",
                             keep.extra.columns = T)
  R2.gr$hg38seq = DNAStringSet(as.character(Views(Hsapiens, R2.gr)))
  # R2.gr

  ## back to dataframe
  R2.hg38seq =
    data.frame(chr = R2.gr %>% seqnames(),
               start = R2.gr %>% start(),
               end = R2.gr %>% end(), ## does start end have to reverse for R2??
               umi = R2.gr$V10, ## R2 is V10, uses the other UMI
               sequence = R2.gr$hg38seq,
               spacer = R2.gr$V4,
               phred = R2.gr$V5,
               strand = R2.gr$V11) %>%
    mutate(sequence_umi = if_else(R2.gr$V11 == "+", paste0(umi, "T", R2.gr$hg38seq),
                                  if_else(R2.gr$V11 == "-", paste0(umi, "T", reverseComplement(R2.gr$hg38seq)),
                                          "ERROR"))) ## { if -ve strand in bam, have to rev_comp to turn it into FASTQ sequence }
  # R2.hg38seq$sequence_umi %>% nchar() %>% unique() # [1] 75
  # R2.hg38seq %>% head()
}

# R2.joined.BACs %>% head()
R2.joined.BACs.rearranged =
  data.frame(chr = R2.joined.BACs$V6,
             start = R2.joined.BACs$V7,
             end = R2.joined.BACs$V8,
             umi = R2.joined.BACs$V10,
             sequence_umi = R2.joined.BACs$V3,
             spacer = "+",
             phred = R2.joined.BACs$V5,
             strand = R2.joined.BACs$V11)
# R2.joined.BACs.rearranged %>% head()
# R2.joined.BACs.rearranged %>% filter(sequence_umi == "CCCTTACTTCCATTTTCTCTAAAGATCAGCTTAAGCATGTCAGCAATGATTTCTAGTTTGCTTCCAGCGTTGGTC")

if (nrow(R2.newEnd) != 0 & nrow(R2.joined.BACs.rearranged) != 0){
  R2.hg38seq.BAC.rbound = bind_rows(R2.hg38seq, R2.joined.BACs.rearranged)
} else {
  R2.hg38seq.BAC.rbound = R2.joined.BACs.rearranged
}
# R2.hg38seq.BAC.rbound$chr %>% unique()
# R2.hg38seq.BAC.rbound %>% head()
# R2.hg38seq.BAC.rbound %>% filter(chr == "F19K16" | chr == "F24B22") %>% head()
# R2.hg38seq.BAC.rbound %>% filter(sequence_umi == "GGAGTCCTGAGGCAACTATGAGTCCGTACTACATGAAGAGCATCACACCGAGTGCTTATGAAGCCACACTTCGAC")

## make dummy frag names
R2_sample_index = rep(paste0("2:N:0:", sample_index), vec_len)
R2_frag_names = paste0(instrument_flowcell, ":", flowcell_coord, " ", R2_sample_index)
# R2_frag_names %>% head()

## add frag names to R2
R2.hg38seq.frag_name =
  R2.hg38seq.BAC.rbound %>%
  dplyr::mutate(frag_names = R2_frag_names) %>%
  dplyr::select(frag_names, sequence_umi, spacer, phred)
# R2.hg38seq.frag_name %>% head()

R2.hg38seq.fastq =
  R2.hg38seq.frag_name %>%
  pivot_longer(cols = everything(), names_to = "col1", values_to = "col2") %>%
  dplyr::select(col2)
# R2.hg38seq.fastq %>% head(12)
# R2.hg38seq.fastq %>% tail(12)

message(paste0("Writing fastq for ", R2.out.fname))
write.table(R2.hg38seq.fastq,
            file = paste0(out.dir, R2.out.fname),
            quote = F, sep = "", row.names = F, col.names = F)

# EOF
