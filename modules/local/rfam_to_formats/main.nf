process RFAM_TO_FORMATS {
    tag "rfam_to_formats"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    path hits_tsv

    output:
    path "annotation.gtf",  emit: gtf
    path "annotation.gff3", emit: gff3
    path "annotation.bed",  emit: bed
    path "versions.yml",    emit: versions

    script:
    """
    rfam_to_formats.py \\
        --hits          ${hits_tsv} \\
        --output_prefix annotation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
