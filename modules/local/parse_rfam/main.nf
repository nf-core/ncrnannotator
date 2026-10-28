process PARSE_RFAM {
    tag "parse_rfam"
    label 'process_medium'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    path tblouts, stageAs: "tblouts/*"
    path rfam_cm
    path rfam_seed

    output:
    path "rfam_hits.tsv", emit: results
    path "versions.yml",  emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    parse_rfam_results.py \\
        --tblout_dir tblouts/ \\
        --rfam_cm    ${rfam_cm} \\
        --rfam_seed  ${rfam_seed} \\
        --output     rfam_hits.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
