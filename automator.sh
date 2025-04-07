#!/bin/bash
export PATH="/usr/local/bin:/opt/homebrew/bin:$PATH"

set -e
set -u
set -o pipefail

# Adjust these if needed:
THREADS=16
STAR_RAM=128000000000   # 128GB for STAR genome generation (in bytes)

# Base directories
INPUT_DIR="$HOME/Desktop/RNAseq"                   # Folder where FASTQ files are placed by user
PIPELINE_DIR="$HOME/Desktop/rnaseq_pipeline"       # Main output folder
RAW_DIR="$PIPELINE_DIR/raw_data"
TRIM_DIR="$PIPELINE_DIR/trimmed_data"
FASTQC_DIR="$PIPELINE_DIR/fastqc_results"
ALIGN_DIR="$PIPELINE_DIR/alignment"
COUNTS_DIR="$PIPELINE_DIR/counts"
ANALYSIS_DIR="$PIPELINE_DIR/analysis"
REF_DIR="$PIPELINE_DIR/reference"                  # To store reference genome and GTF
STAR_INDEX_DIR="$PIPELINE_DIR/star_index"

# 1. Create output directories
echo "Setting up output directories..."
mkdir -p "$RAW_DIR" "$TRIM_DIR" "$FASTQC_DIR" "$ALIGN_DIR" "$COUNTS_DIR" "$ANALYSIS_DIR" "$REF_DIR" "$STAR_INDEX_DIR"

# 2. Move or copy FASTQ files into raw_data
if compgen -G "$INPUT_DIR/*_R1.fastq.gz" > /dev/null; then
    echo "Importing FASTQ files from $INPUT_DIR to $RAW_DIR..."
    # Copy (or move) FASTQ files; using copy to preserve originals
    cp -v "$INPUT_DIR"/*_R{1,2}.fastq.gz "$RAW_DIR/" 2>/dev/null || true
else
    echo "No FASTQ files found in $INPUT_DIR. Please add paired FASTQ.gz files and re-run."
    exit 1
fi

# 3. Download mouse GRCm39 genome and GENCODE annotation if not already present
GENOME_FA="$REF_DIR/GRCm39.primary_assembly.genome.fa"
GTF_FILE="$REF_DIR/gencode.vM36.primary_assembly.annotation.gtf"

if [ ! -f "$GENOME_FA" ]; then
    echo "Downloading GRCm39 (mm39) genome FASTA..."
    curl -s -S -L -o "${GENOME_FA}.gz" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz"
    gunzip "${GENOME_FA}.gz"
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "Downloading GENCODE vM36 annotation GTF..."
    curl -s -S -L -o "${GTF_FILE}.gz" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.annotation.gtf.gz"
    gunzip "${GTF_FILE}.gz"
fi

# 4. Build STAR index if not present
if [ ! -f "$STAR_INDEX_DIR/Genome" ]; then
    echo "Building STAR index for GRCm39 genome (this may take some time)..."
    STAR --runThreadN $THREADS \
         --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$GENOME_FA" \
         --sjdbGTFfile "$GTF_FILE" \
         --sjdbOverhang 100 \
         --limitGenomeGenerateRAM $STAR_RAM
fi

# Download Trimmomatic adapter file (Illumina TruSeq adapters) if not present
ADAPTER_FA="$PIPELINE_DIR/TruSeq3-PE.fa"
if [ ! -f "$ADAPTER_FA" ]; then
    echo "Downloading Trimmomatic adapter sequences..."
    curl -s -S -L -o "$ADAPTER_FA" "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa"
fi

# 5. Process each FASTQ pair: FastQC -> Trimmomatic -> STAR -> Samtools -> FeatureCounts
echo "Starting per-sample processing..."
cd "$RAW_DIR"
for R1_FILE in *_R1.fastq.gz; do
    # Derive sample name by stripping the _R1.fastq.gz suffix
    SAMPLE="${R1_FILE%_R1.fastq.gz}"
    R2_FILE="${SAMPLE}_R2.fastq.gz"
    echo "Processing sample: $SAMPLE"

    # Quality check on raw reads
    echo "  Running FastQC on raw reads..."
    fastqc -t $THREADS -o "$FASTQC_DIR" "$R1_FILE" "$R2_FILE"

    # Trim adapters and low-quality sequences with Trimmomatic
    echo "  Trimming adapters and low-quality bases..."
    TRIM_R1_PAIRED="$TRIM_DIR/${SAMPLE}_R1.trimmed.fastq.gz"
    TRIM_R1_UNPAIRED="$TRIM_DIR/${SAMPLE}_R1.unpaired.fastq.gz"
    TRIM_R2_PAIRED="$TRIM_DIR/${SAMPLE}_R2.trimmed.fastq.gz"
    TRIM_R2_UNPAIRED="$TRIM_DIR/${SAMPLE}_R2.unpaired.fastq.gz"
    trimmomatic PE -threads $THREADS -phred33 \
        "$R1_FILE" "$R2_FILE" \
        "$TRIM_R1_PAIRED" "$TRIM_R1_UNPAIRED" \
        "$TRIM_R2_PAIRED" "$TRIM_R2_UNPAIRED" \
        ILLUMINACLIP:"$ADAPTER_FA":2:30:10:8:TRUE SLIDINGWINDOW:4:20 MINLEN:36

    # (Optional) Run FastQC on trimmed reads as well
    echo "  Running FastQC on trimmed reads..."
    fastqc -t $THREADS -o "$FASTQC_DIR" "$TRIM_R1_PAIRED" "$TRIM_R2_PAIRED"

    # Align with STAR
    echo "  Aligning reads to genome with STAR..."
    STAR --runThreadN $THREADS \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$TRIM_R1_PAIRED" "$TRIM_R2_PAIRED" \
         --readFilesCommand "gunzip -c" \
         --outFileNamePrefix "$ALIGN_DIR/${SAMPLE}." \
         --outSAMtype BAM Unsorted

    # Sort alignments and index
    echo "  Sorting and indexing alignments..."
    samtools sort -@ $THREADS -m 16G -o "$ALIGN_DIR/${SAMPLE}.bam" "$ALIGN_DIR/${SAMPLE}.Aligned.out.bam"
    samtools index "$ALIGN_DIR/${SAMPLE}.bam"
    rm "$ALIGN_DIR/${SAMPLE}.Aligned.out.bam"  # remove unsorted BAM to save space
done

# FeatureCounts: count reads for all samples in one run
echo "Running featureCounts on all alignment files..."
featureCounts -T $THREADS -p -a "$GTF_FILE" -o "$COUNTS_DIR/counts.txt" "$ALIGN_DIR"/*.bam  \
    && echo "Gene count matrix saved to $COUNTS_DIR/counts.txt"

# 6. Differential expression with DESeq2 (via Rscript)
echo "Running DESeq2 differential expression analysis in R..."
# Create an R script dynamically:
cat > "$ANALYSIS_DIR/deseq2_analysis.R" << 'EOF'
#!/usr/bin/env Rscript
suppressMessages(library(DESeq2))
counts <- read.table(file.path(Sys.getenv("PIPELINE_DIR"), "counts/counts.txt"), 
                     header=TRUE, row.names=1, comment.char="#")
counts <- counts[ , grep(".bam$", colnames(counts)) ]
colnames(counts) <- sub(".bam$", "", basename(colnames(counts)))
sampleNames <- colnames(counts)
conditions <- ifelse(grepl("control", sampleNames, ignore.case=TRUE), "control", "treatment")
colData <- data.frame(row.names = sampleNames, condition = factor(conditions))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue), ]
write.csv(as.data.frame(resOrdered), file=file.path(Sys.getenv("PIPELINE_DIR"), "analysis/deseq2_results.csv"))
EOF

# Run the DESeq2 R script
PIPELINE_DIR="$PIPELINE_DIR" Rscript "$ANALYSIS_DIR/deseq2_analysis.R"

echo "DESeq2 analysis complete. Results saved in $ANALYSIS_DIR (e.g., 'deseq2_results.csv')."
echo "RNA-seq pipeline finished."
