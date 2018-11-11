# WDL Workflow for standard UofA WGBS pipeline

# Jimmy Breen (jimmymbreen@gmail.com)

## Runs:
##  - Trimming using trim_galore
##  - Mapping using Bismark
##  - Sort BAM and mark duplicates

## Validated the workflow using wdltool
## - And made acyclic graph using:
# dot <(wdltool graph wgbs.wdl) -Tpdf > wgbs.pdf

workflow WGBS {

    # Raw Sequencing Reads
    Array[File] filesR1
    Array[File] filesR2

    # Reference Genome Files
    File genomeFile
    String tempDir
    String index
    File PileOMeth
    File bismark
    File debismark

    # Read Trimming
    Int clipR1 = 6
    Int clipR2 = 6

    Int cpu

    scatter (pairedFiles in zip(filesR1, filesR2)) {
        call TrimGalorePaired {
            input:
                pairedFiles = pairedFiles,
                clipR1 = clipR1,
                clipR2 = clipR2
    }

        call BismarkPaired {
            input:
                index = index,
                bismark = bismark,
                pairedFiles = TrimGalorePaired.trimFiles,
                tempDir = tempDir,
                cpu = cpu
    }

        call Deduplicate {
            input:
                file = BismarkPaired.outputFile,
                debismark = debismark
    }

        call Sort_Bam {
            input:
                file = Deduplicate.outputFile,
                sambamba = sambamba,
                cpu = cpu
    }

        call PileOMeth_extract {
            input:
                PileOMeth = PileOMeth,
                genomeFile = genomeFile,
                file = Sort_Bam.outputFile
    }
}

    output {
        Array[File] outputFile = PileOMeth_extract.outputFile

        # Report files
        Array[Pair[File, File]] trimReportFiles = TrimGalorePaired.statsFiles
        Array[File] alignReportFile = BismarkPaired.reportFile
        Array[File] deduplicatedReportFiles = Deduplicate.reportFile
    }
}

## Tasks
# Trim reads using AdapterRemoval
task TrimGalorePaired {
    # General options
    Pair[File, File] pairedFiles
    String? adapter
    String? adapter2
    Boolean gzip = true
    Int? clipR1
    Int? clipR2

    command {
            trim_galore ${pairedFiles.left} ${pairedFiles.right} \
                    --phred33 --gzip --paired \
                    ${'--adapter ' + adapter} \
                    ${'--adapter2 ' + adapter2} \
                    ${'--clip_R1 ' + clipR1} \
                    ${'--clip_R2 ' + clipR2} \
                    ${true='--gzip' false='--dont_gzip' gzip}
    }

    output {
            Pair[File, File] trimFiles = (
                sub(basename(pairedFiles.left), ".fastq(.gz)?", "_val_1.fq") + if (gzip) then ".gz" else "",
                sub(basename(pairedFiles.right), ".fastq(.gz)?", "_val_2.fq") + if (gzip) then ".gz" else "")
            Pair[File, File] statsFiles = (
                basename(pairedFiles.left) + "_trimming_report.txt",
                basename(pairedFiles.right) + "_trimming_report.txt")
        }
}

task BismarkPaired {
    File index
    File bismark
    Pair[File, File] pairedFiles
    Int cpu
    String tempDir

  command {
      ${bismark} ${index} -1 ${pairedFiles.left} -2 ${pairedFiles.right} \
        --fastq --bowtie2 -p ${cpu} --temp_dir ${tempDir}
  }

  output {
      File outputFile = glob("*_pe.bam")[0]
      File reportFile = glob("*_PE_report.txt")[0]
      Pair[File, File] ambiguousFiles = (
          "${basename(pairedFiles.left)}_ambiguous_reads_1.fq.gz",
          "${basename(pairedFiles.right)}_ambiguous_reads_2.fq.gz")
      Pair[File, File] unmappedFiles = (
          "${basename(pairedFiles.left)}_unmapped_reads_1.fq.gz",
          "${basename(pairedFiles.right)}_unmapped_reads_2.fq.gz")
  }
}

task Deduplicate {
    File debismark
    File file

    command {
      ${debismark} ${file} --paired --bam
  }

  output {
    File outputFile = sub(basename(file), ".bam$", ".deduplicated.bam")
    File reportFile = sub(basename(file), ".bam$", ".deduplication_report.txt")
  }
}

task Sort_Bam {
    File sambamba
    File file
    Int cpu

    command {
      ${sambamba} sort -t ${cpu} ${file}
  }

  output {
    File outputFile = sub(basename(file), ".bam$", ".sorted.bam")
  }
}

task PileOMeth_extract {
    File PileOMeth
    File file
    File genomeFile

  command {
      ${PileOMeth} extract -d 10 ${genomeFile} ${file}
        }

  output {
      File outputFile = sub(basename(file), ".sorted.bam$", ".CpG.bedGraph.gz")
  }
}
