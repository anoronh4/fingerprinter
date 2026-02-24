/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_fingerprinter_pipeline'
include {
        FINGERPRINT_GBCMS as FINGERPRINT_GBCMS_GRCH38 ;
        FINGERPRINT_GBCMS as FINGERPRINT_GBCMS_GRCH37 ;
        FINGERPRINT_GBCMS as FINGERPRINT_GBCMS_GRCH37_CHR
    } from '../subworkflows/msk/fingerprint_gbcms/main'
include { FINGERPRINT_GBCMS_BATCH } from '../subworkflows/msk/fingerprint_gbcms_batch/main'
include { FASTAREMOVECHRPREFIX    } from '../modules/local/fastaremovechrprefix/main'
include { FASTAADDCHRPREFIX       } from '../modules/local/fastaaddchrprefix/main'
include {
    SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GRCH38 ;
    SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GRCH37_CHR
} from '../modules/nf-core/samtools/faidx/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FINGERPRINTER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    FASTAREMOVECHRPREFIX(Channel.of([[:],params.fasta_grch38]))
    SAMTOOLS_FAIDX_GRCH38(FASTAREMOVECHRPREFIX.out.fasta,[[:],[]],true)

    FASTAADDCHRPREFIX(Channel.of([[:],params.fasta_grch37]))
    SAMTOOLS_FAIDX_GRCH37_CHR(FASTAADDCHRPREFIX.out.fasta,[[:],[]],true)

    ch_is_bam = ch_samplesheet.filter{ meta, bam_or_fp, bai ->
        meta.is_bam
    }
    ch_bam = ch_is_bam.map{ meta, bam_or_fp, bai ->
        [meta, bam_or_fp]
    }
    ch_bai = ch_is_bam.map{ meta, bam_or_fp, bai ->
        [meta, bai]
    }

    ch_fp_tsv = ch_samplesheet.filter{ meta, bam_or_fp, bai ->
        ! meta.is_bam
    }.map{ meta, bam_or_fp, bai ->
        [ meta, bam_or_fp ]
    }

    FINGERPRINT_GBCMS_GRCH38(
        ch_bam.filter{ meta, bam -> meta.genome == "GRCh38" },
        ch_bai.filter{ meta, bai -> meta.genome == "GRCh38" },
        ch_fp_tsv.filter{ meta, fp_tsv -> meta.genome == "GRCh38" },
        Channel.of([[:],file(workflow.projectDir + "/assets/combined.dbsnp.assayFP.18517.grch38.vcf")]),
        [],
        FASTAREMOVECHRPREFIX.out.fasta.map{it[1]},
        SAMTOOLS_FAIDX_GRCH38.out.fai.map{it[1]},
        "GRCh38",
        false
    )

    FINGERPRINT_GBCMS_GRCH37(
        ch_bam.filter{ meta, bam -> meta.genome == "GRCh37" },
        ch_bai.filter{ meta, bai -> meta.genome == "GRCh37" },
        ch_fp_tsv.filter{ meta, fp_tsv -> meta.genome == "GRCh37" },
        Channel.of([[:],file(workflow.projectDir + "/assets/combined.dbsnp.assayFP.18557.grch37.vcf")]),
        [],
        Channel.of(params.fasta_grch37),
        Channel.of(params.fasta_grch37_fai),
        "GRCh37",
        false
    )

    FINGERPRINT_GBCMS_GRCH37_CHR(
        ch_bam.filter{ meta, bam -> meta.genome == "GRCh37chr" },
        ch_bai.filter{ meta, bai -> meta.genome == "GRCh37chr" },
        ch_fp_tsv.filter{ meta, fp_tsv -> meta.genome == "GRCh37chr" },
        Channel.of([[:],file(workflow.projectDir + "/assets/combined.dbsnp.assayFP.18557.grch37.chr.vcf")]),
        [],
        FASTAADDCHRPREFIX.out.fasta.map{it[1]},
        SAMTOOLS_FAIDX_GRCH37_CHR.out.fai.map{it[1]},
        "GRCh37chr",
        false
    )

    FINGERPRINT_GBCMS_BATCH(
        FINGERPRINT_GBCMS_GRCH38.out.fp_tsv
            .mix(FINGERPRINT_GBCMS_GRCH37.out.fp_tsv)
            .mix(FINGERPRINT_GBCMS_GRCH37_CHR.out.fp_tsv)
            .map{ meta, fp ->
                def meta2 = meta.clone()
                meta2.genome = meta.genome.replace("chr","")
                [meta2, fp]
            },
        Channel.of(file(workflow.projectDir + "/assets/liftover_mapping.tsv")),
        "GRCh37",
        params.pool ? channel.of(params.pool,"") : []


    )

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'fingerprinter_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
