process CMSEARCH {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::infernal=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h7b50bb2_4' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h7b50bb2_4' }"

    input:
    tuple val(meta), path(chunk)
    path rfam_cm

    output:
    tuple val(meta), path("*.tblout"), emit: tblout
    path "versions.yml",               emit: versions

    script:
    def args    = task.ext.args ?: '--rfam --nohmmonly --cut_ga'
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    cmsearch \\
        ${args} \\
        --cpu ${task.cpus} \\
        --tblout ${prefix}.tblout \\
        ${rfam_cm} \\
        ${chunk}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmsearch -h 2>&1 | grep -m 1 "# INFERNAL" | sed 's/.*INFERNAL //; s/ .*//')
    END_VERSIONS
    """
}
