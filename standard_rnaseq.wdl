# WDL Workflow for standard UofA RNAseq pipeline

# Jimmy Breen (jimmymbreen@gmail.com)

## Runs:
##  - Trimming using AdapterRemoval
##  - Mapping using STAR
##  - Pseudomapping using kallisto
##  - Quant using featureCounts

## Validated the workflow using wdltool
## - And made acyclic graph using:
# dot <(wdltool graph standard_rnaseq.wdl) -Tpdf > standard_rnaseq_wdl_validate.pdf

workflow standard_rnaseq_quant {

    # Tools
    File AdapterRemoval
    File STAR
    File featureCounts

    # Annotation and indexes
    File annotation_file
    String STARindexDir

    # Data
    File sampleinfo
    String datadir
    Array[String] sample_name = read_lines(sampleinfo)

    scatter(idx in range(length(sample_name))) {
        call AdapterRemoval_Trim {
            input:
                AdapterRemoval = AdapterRemoval,
                sample_name = sample_name[idx],
                rawfastq_read1 = datadir+'/'+sample_name[idx]+"_R1.fastq.gz"#,
            }

        call rnaseq_star_map {
            input:
                STAR = STAR,
                STARindexDir = STARindexDir,
                sample_name = sample_name[idx],
                trim_fastqR1 = AdapterRemoval_Trim.out_trim_R1
            }

        call rnaseq_featureCounts_quant {
            input:
                featureCounts = featureCounts,
                annotation_file = annotation_file,
                sample_name = sample_name[idx],
                in_bam = rnaseq_star_map.out_bam
            }
    }
}

## Tasks
# Trim reads using AdapterRemoval
task AdapterRemoval_Trim {
    File AdapterRemoval
    String sample_name
    File rawfastq_read1
    Int cpu=12

    command {
        module load AdapterRemoval/2.2.0-foss-2016a
        ${AdapterRemoval} --file1 ${rawfastq_read1} \
            --output1 ${sample_name}_T1.fastq.gz \
            --trimns --trimqualities --gzip
        }

    output {
        # nice little hack from the ENCODE guys
        File out_trim_R1 = "${sample_name}_T1.fastq.gz"
        }

    runtime {
                cpu: cpu
        }
}

# Map RNAseq reads:
task rnaseq_star_map {
    File STAR
    String STARindexDir
    String sample_name
    File trim_fastqR1
    String genome="hg19"
    Int cpu=12

    command {
        ${STAR} --genomeDir ${STARindexDir} --readFilesIn ${trim_fastqR1} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample_name}.${genome}. --outSAMattrRGline ID:${sample_name} LB:library PL:illumina PU:machine SM:${genome} --outSAMmapqUnique 60 --runThreadN ${cpu}
        }

    output {
        File out_bam = "${sample_name}.${genome}.Aligned.sortedByCoord.out.bam"
    }
    runtime {
        cpu: cpu
  }
}

## featureCounts
task rnaseq_featureCounts_quant {
    File featureCounts
    File annotation_file
    String sample_name
    File in_bam
    Int cpu=6

    command {
        ${featureCounts} -a ${annotation_file} \
            -T ${cpu} \
            -s 1 \
            -o ${sample_name}.geneCounts.tsv \
            ${in_bam} > ${sample_name}.featureCounts.log
        }

    output {
        File out_counts = "${sample_name}.geneCounts.tsv"
        String out_log = "${sample_name}.featureCounts.log"
    }

    runtime {
        cpu: cpu
    }
}            
