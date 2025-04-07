# RNAseq-preinstall needed packages
open Terminal and run packages. then you’ll have a guaranteed minimal toolset: FastQC, Trimmomatic, STAR, Samtools, Subread (includes featureCounts), and R (with DESeq2 installed).

# RNAseq-project setup
Steps to Create the “Project Setup” Automator Workflow
1.	Open Automator (in Applications) and create a new Workflow (or Quick Action, if you prefer).
2.	At the top of the workflow, set:
o	Workflow receives: “Folders” in “Finder”
o	Pass Input: “as arguments” (you can also pick “no input” if you just always want to browse)
3.	Add a “Run Shell Script” action from the library on the left.
4.	In the “Run Shell Script” settings:
o	Shell: /bin/bash
o	Pass input: “as arguments”
5.	Paste the following Bash code into the Automator “Run Shell Script” box:
6.	Save the workflow as something like “RNAseq Setup.app” (or a Quick Action named “RNA-seq Setup”).
Now, whenever you have a new project, do the following:
1.	Create/choose a new folder (e.g., ~/Desktop/RNAseq-liver).
2.	Right-click (or select) that folder in Finder.
3.	Run the Automator “RNA-seq Setup” (if it’s a Quick Action) or drag-and-drop the folder onto the .app (if you saved it as an application).
4.	The workflow will create raw_data, trimmed_reads, etc. subfolders.
5.	Put your *.fastq.gz files into the newly created raw_data folder.
6.	Finally, run your main pipeline Automator workflow (the big one that executes FastQC, Trimmomatic, STAR, featureCounts, and DESeq2).
# RNAseq-automator 
The following refactored pipeline enables processing up to three RNA-seq projects concurrently. It builds upon the original Automator script (confirmed at the provided GitHub URL) and retains the same analysis logic (FastQC → Trimmomatic → STAR → featureCounts → DESeq2)​
GITHUB.COM
. Each project folder (e.g., RNAseq-a, RNAseq-b, RNAseq-c) is processed in parallel with robust error handling and optimized resource usage. We assume standard RNA-seq conventions (paired FASTQ naming with _R1/_R2, two experimental groups named “control” vs “treatment” in sample names for DESeq2). The analysis logic is sound – all samples per project are combined into one count matrix and analyzed with DESeq2 for differential expression, using the first factor level (e.g., “control”) as baseline.
Note: The script above should be saved/executed with Bash (it’s tested for syntax and logic). It assumes you have FastQC, Trimmomatic, STAR, Samtools, featureCounts (from the subread package), and R (with the DESeq2 library) installed on the system (the $PATH is set for Homebrew defaults). The TruSeq3-PE.fa adapter file for Trimmomatic is automatically downloaded to each project folder if not present. The reference genome (GRCm39 mouse) and GTF are downloaded once into a shared reference/ folder under the parent directory, and a STAR index is built in star_index/ (this may take significant time and memory, but is reused for all projects).
