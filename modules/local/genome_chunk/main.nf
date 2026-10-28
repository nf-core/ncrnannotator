process GENOME_CHUNK {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(fasta)
    val(chunk_size)

    output:
    path "chunks/*.fa", emit: chunks
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p chunks

    python3 - <<'PYEOF'
    import gzip
    import os
    import re

    chunk_size = int("${chunk_size}")
    fasta_path = "${fasta}"

    def open_fasta(path):
        if path.endswith('.gz'):
            return gzip.open(path, 'rt')
        return open(path)

    def sanitise(name):
        # Replace characters that are problematic in filenames
        return re.sub(r'[^A-Za-z0-9._\\-]', '_', name)

    with open_fasta(fasta_path) as fh:
        seqname = None
        seq_parts = []

        def flush_seq(seqname, sequence):
            if not seqname or not sequence:
                return
            seq_len = len(sequence)
            pos = 0
            while pos < seq_len:
                start = pos + 1          # 1-based
                end   = min(pos + chunk_size, seq_len)  # 1-based inclusive
                chunk = sequence[pos:end]
                safe  = sanitise(seqname)
                fname = f"chunks/{safe}.rs{start}.re{end}.fa"
                header = f">{seqname}.rs{start}.re{end}"
                with open(fname, 'w') as out:
                    out.write(header + '\\n')
                    # Write 60 bp per line
                    for i in range(0, len(chunk), 60):
                        out.write(chunk[i:i+60] + '\\n')
                pos += chunk_size

        for line in fh:
            line = line.rstrip('\\n')
            if line.startswith('>'):
                flush_seq(seqname, ''.join(seq_parts))
                seqname = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        flush_seq(seqname, ''.join(seq_parts))

    print(f"genome_chunk: created {len(os.listdir('chunks'))} chunk files")
    PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
