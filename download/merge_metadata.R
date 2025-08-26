#!/usr/bin/env Rscript

## CONTROLLED ACCESS DATA- these files with phenotypes can't be used outside the controlled access file space

subject = read.delim("data/phenotypic_tsvs_controlled_access/subject.tsv")
sample = read.delim("data/phenotypic_tsvs_controlled_access/sample.tsv")
sequencing = read.delim("data/phenotypic_tsvs_controlled_access/sequencing.tsv")
bamfiles <- list.files("data/rnaseq_aligned_reads/", pattern = '.bam$')
bam_sequencing <- unique(sequencing[sequencing$file_name %in% bamfiles,])

# Do some field parsing to facilitate joins

sample$subject_id <- sub('(GTEX-[^\\-]+)-[^-]+.*', '\\1', sample$submitter_id)
bam_sequencing$specimen_id <- sub("(.*).Aligned.*", "\\1", bam_sequencing$file_name)

# Do some joins
sample_sequencing = merge(sample, bam_sequencing, by = 'specimen_id', suffixes = c(".sample", ".file"))
subject_sample_sequencing <- merge(subject, sample_sequencing, by.x = "submitter_id", by.y="subject_id", suffixes = c(".subject",".sample"))

# Write the joined output
write.table(subject_sample_sequencing, sep = '\t', file = "data/phenotypic_tsvs_controlled_access/subject_sample_sequencing.tsv", quote = FALSE, row.names = FALSE)

## PUBLIC METADATA- these were retrieved from the public GTEX so metadata can be displayed publicly

sample <- read.delim("data/phenotypic_tsvs_public/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
subject <- read.delim("data/phenotypic_tsvs_public/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

sample$SUBJID <- sub('(GTEX-[^\\-]+)-[^-]+.*', '\\1', sample$SAMPID)
subject_sample <- merge(subject, sample, by="SUBJID")
subject_sample_sequencing <- merge(subject_sample, bam_sequencing, by.x = "SAMPID", by.y = "specimen_id")

# Write the joined output
write.table(subject_sample_sequencing, sep = '\t', file = "data/phenotypic_tsvs_public/subject_sample_sequencing.tsv", quote = FALSE, row.names = FALSE)

