/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FILTER_RFAM_CM         } from '../modules/local/filter_rfam_cm/main'
include { GENOME_CHUNK           } from '../modules/local/genome_chunk/main'
include { CMSEARCH               } from '../modules/local/cmsearch/main'
include { PARSE_RFAM             } from '../modules/local/parse_rfam/main'
include { RFAM_TO_FORMATS        } from '../modules/local/rfam_to_formats/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ncrnannotator_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NCRNANNOTATOR {

    take:
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    // -----------------------------------------------------------------------
    // Resolve chunk size (mode-dependent default, user override via params)
    // -----------------------------------------------------------------------
    def chunk_size = params.chunk_size
        ? params.chunk_size as Integer
        : (params.mode == 'ensembl-vertebrates' ? 1_000_000
        :  params.mode == 'full'                ? 1_000_000
        :  params.mode == 'mgnify-assembly'     ? 50_000_000
        :                                         100_000)

    // -----------------------------------------------------------------------
    // Build input channels
    // -----------------------------------------------------------------------
    ch_genome = channel.value([
        [id: file(params.fasta).baseName],
        file(params.fasta, checkIfExists: true)
    ])

    ch_rfam_cm   = file(params.rfam_cm,   checkIfExists: true)
    ch_rfam_seed = file(params.rfam_seed, checkIfExists: true)

    // -----------------------------------------------------------------------
    // Step 1: Filter Rfam.cm to clade-specific accessions
    //         (skipped in mgnify-assembly and full modes — full Rfam.cm used directly)
    // -----------------------------------------------------------------------
    if (params.mode in ['ensembl-vertebrates', 'ensembl-invertebrates']) {
        def accession_basename = params.mode == 'ensembl-vertebrates' ? 'vertebrates' : 'invertebrates'
        ch_accessions = params.rfam_accessions
            ? file(params.rfam_accessions, checkIfExists: true)
            : file("${projectDir}/assets/rfam_accessions/${accession_basename}.txt",
                   checkIfExists: true)

        FILTER_RFAM_CM(ch_rfam_cm, ch_accessions)
        ch_versions = ch_versions.mix(FILTER_RFAM_CM.out.versions)
        ch_filtered_cm = FILTER_RFAM_CM.out.filtered_cm
    } else {
        ch_filtered_cm = ch_rfam_cm
    }

    // -----------------------------------------------------------------------
    // Step 2: Chunk genome into windows
    // -----------------------------------------------------------------------
    GENOME_CHUNK(ch_genome, channel.value(chunk_size))
    ch_versions = ch_versions.mix(GENOME_CHUNK.out.versions)

    // Scatter: one (meta, chunk_file) entry per chunk
    ch_chunks = GENOME_CHUNK.out.chunks
        .flatten()
        .map { f -> [ [id: f.baseName], f ] }

    // -----------------------------------------------------------------------
    // Step 3: cmsearch — runs in parallel on each chunk
    // -----------------------------------------------------------------------
    CMSEARCH(ch_chunks, ch_filtered_cm)
    ch_versions = ch_versions.mix(CMSEARCH.out.versions.first())

    // -----------------------------------------------------------------------
    // Step 4: Parse all tblout files, remove overlaps, apply GA thresholds
    // -----------------------------------------------------------------------
    PARSE_RFAM(
        CMSEARCH.out.tblout.collect { it[1] },
        ch_filtered_cm,
        ch_rfam_seed
    )
    ch_versions = ch_versions.mix(PARSE_RFAM.out.versions)

    // -----------------------------------------------------------------------
    // Step 5: Convert hits to GTF / GFF3 / BED
    // -----------------------------------------------------------------------
    RFAM_TO_FORMATS(PARSE_RFAM.out.results)
    ch_versions = ch_versions.mix(RFAM_TO_FORMATS.out.versions)

    // -----------------------------------------------------------------------
    // Collate software versions
    // -----------------------------------------------------------------------
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
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'ncrnannotator_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'ncrnannotator'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    gtf            = RFAM_TO_FORMATS.out.gtf
    gff3           = RFAM_TO_FORMATS.out.gff3
    bed            = RFAM_TO_FORMATS.out.bed
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
