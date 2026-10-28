process FILTER_RFAM_CM {
    tag "rfam_filter"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    path rfam_cm
    path accessions

    output:
    path "rfam_filtered.cm", emit: filtered_cm
    path "versions.yml",     emit: versions

    script:
    """
    filter_rfam_cm.py \\
        --rfam_cm ${rfam_cm} \\
        --accessions ${accessions} \\
        --output rfam_filtered.cm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
