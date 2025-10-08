process SALMON_INDEX {
    tag "salmon_index"
    conda "${projectDir}/envs/salmon.yml"
    publishDir "${params.reference.salmon_index}", mode: 'copy', overwrite: true
    cpus params.salmon?.threads ?: 4

    // Enable Nextflow caching
    cache true
    maxRetries 1

    input:
        tuple path(transcripts), path(annotation), path(decoy, optional: true)

    output:
        path "salmon_index", emit: index
        path "tx2gene.tsv", emit: tx2gene

    script:
        """
        set -euo pipefail

        # Check custom cache if enabled
        if [ "${CACHE_ENABLED}" = "true" ]; then
            # Generate cache key for this process
            input_hash=\$(python3 -c "
import sys
sys.path.append('${projectDir}/../../scripts')
from cache_manager import get_cache_manager

cache_mgr = get_cache_manager('${CACHE_DIR}', True)
input_files = ['${transcripts}', '${annotation}']
${decoy ? "input_files.append('${decoy}')" : ""}
params = {'decoy': '${decoy}', 'threads': ${task.cpus}}

input_hash = cache_mgr.compute_input_hash([__import__('pathlib').Path(f) for f in input_files], params)
print(input_hash[:16])
            ")

            cache_key="salmon_index_\${input_hash}"
            echo "Checking cache for key: \${cache_key}"

            # Check if cached result exists
            cache_check=\$(python3 -c "
import sys
sys.path.append('${projectDir}/../../scripts')
from cache_manager import get_cache_manager

cache_mgr = get_cache_manager('${CACHE_DIR}', True)
cache_key = 'salmon_index_\${input_hash}'
output_files = ['${projectDir}/${params.reference.salmon_index}/salmon_index', '${projectDir}/${params.reference.salmon_index}/tx2gene.tsv']

if cache_mgr.is_cached(cache_key, [__import__('pathlib').Path(f) for f in output_files]):
    print('CACHE_HIT')
else:
    print('CACHE_MISS')
" 2>/dev/null || echo "CACHE_MISS")

            if [ "\$cache_check" = "CACHE_HIT" ]; then
                echo "Using cached Salmon index"
                mkdir -p '${projectDir}/${params.reference.salmon_index}'
                touch '${projectDir}/${params.reference.salmon_index}/salmon_index'
                touch '${projectDir}/${params.reference.salmon_index}/tx2gene.tsv'
            else
                echo "Building new Salmon index"
                mkdir -p salmon_index
                decoy_flag=""
                if [ -s "${decoy}" ]; then
                  decoy_flag="-d ${decoy}"
                fi
                bash ${projectDir}/scripts/build_salmon_index.sh \
                  -t ${transcripts} \
                  -a ${annotation} \
                  ${decoy_flag} \
                  -o salmon_index \
                  -p ${task.cpus}
                cp salmon_index/tx2gene.tsv tx2gene.tsv

                # Mark as complete in cache
                python3 -c "
import sys
sys.path.append('${projectDir}/../../scripts')
from cache_manager import get_cache_manager

cache_mgr = get_cache_manager('${CACHE_DIR}', True)
cache_key = 'salmon_index_\${input_hash}'
input_files = ['${transcripts}', '${annotation}']
${decoy ? "input_files.append('${decoy}')" : ""}
params = {'decoy': '${decoy}', 'threads': ${task.cpus}}

output_files = ['${projectDir}/${params.reference.salmon_index}/salmon_index', '${projectDir}/${params.reference.salmon_index}/tx2gene.tsv']
cache_mgr.store_cache(cache_key, cache_mgr.compute_input_hash([__import__('pathlib').Path(f) for f in input_files], params), [__import__('pathlib').Path(f) for f in output_files])
                "
            fi
        else
            # Run without custom caching, rely on Nextflow cache
            mkdir -p salmon_index
            decoy_flag=""
            if [ -s "${decoy}" ]; then
              decoy_flag="-d ${decoy}"
            fi
            bash ${projectDir}/scripts/build_salmon_index.sh \
              -t ${transcripts} \
              -a ${annotation} \
              ${decoy_flag} \
              -o salmon_index \
              -p ${task.cpus}
            cp salmon_index/tx2gene.tsv tx2gene.tsv
        fi
        """
}

process SALMON_QUANT {
    tag "${sample}"
    conda "${projectDir}/envs/salmon.yml"
    publishDir { "${params.paths.salmon}/${sample}" }, mode: 'copy', overwrite: true
    cpus params.salmon?.threads ?: 4

    // Enable Nextflow caching
    cache true
    maxRetries 1

    input:
        tuple val(sample), val(meta), path(fastq1), path(fastq2, optional: true)
        path index_dir

    output:
        tuple val(sample), path("quant.sf"), path("lib_format_counts.json")

    script:
        """
        set -euo pipefail

        # Check custom cache if enabled
        if [ "${CACHE_ENABLED}" = "true" ]; then
            # Generate cache key for this sample
            input_hash=\$(python3 -c "
import sys
sys.path.append('${projectDir}/../../scripts')
from cache_manager import get_cache_manager

cache_mgr = get_cache_manager('${CACHE_DIR}', True)
input_files = ['${fastq1}']
${fastq2 && fastq2 != 'null' ? "input_files.append('${fastq2}')" : ""}
input_files.append('${index_dir}')
params = {'libtype': '${params.salmon?.libtype ?: 'A'}', 'extra': '${params.salmon?.extra ?: ''}', 'threads': ${task.cpus}}

input_hash = cache_mgr.compute_input_hash([__import__('pathlib').Path(f) for f in input_files], params)
print(input_hash[:16])
            ")

            cache_key="salmon_quant_${sample}_\${input_hash}"
            echo "Checking cache for key: \${cache_key}"

            sample_output_dir="${params.paths.salmon}/${sample}"

            # Check if cached result exists
            cache_check=\$(python3 -c "
import sys
sys.path.append('${projectDir}/../../scripts')
from cache_manager import get_cache_manager

cache_mgr = get_cache_manager('${CACHE_DIR}', True)
cache_key = 'salmon_quant_${sample}_\${input_hash}'
output_files = ['${sample_output_dir}/quant.sf', '${sample_output_dir}/lib_format_counts.json']

if cache_mgr.is_cached(cache_key, [__import__('pathlib').Path(f) for f in output_files]):
    print('CACHE_HIT')
else:
    print('CACHE_MISS')
" 2>/dev/null || echo "CACHE_MISS")

            if [ "\$cache_check" = "CACHE_HIT" ]; then
                echo "Using cached Salmon quantification for ${sample}"
                mkdir -p '${sample_output_dir}'
                touch '${sample_output_dir}/quant.sf'
                touch '${sample_output_dir}/lib_format_counts.json'
            else
                echo "Running Salmon quantification for ${sample}"
                mkdir -p quant
                if [ -s "${fastq2}" ] && [ "${params.se}" != "true" ]; then
                  salmon quant -i ${index_dir} --libType ${params.salmon?.libtype ?: 'A'} \
                    -1 ${fastq1} -2 ${fastq2} --threads ${task.cpus} ${params.salmon?.extra ?: ''} \
                    -o quant
                else
                  salmon quant -i ${index_dir} --libType ${params.salmon?.libtype ?: 'A'} \
                    -r ${fastq1} --threads ${task.cpus} ${params.salmon?.extra ?: ''} \
                    -o quant
                fi
                mv quant/quant.sf .
                mv quant/lib_format_counts.json .

                # Mark as complete in cache
                python3 -c "
import sys
sys.path.append('${projectDir}/../../scripts')
from cache_manager import get_cache_manager

cache_mgr = get_cache_manager('${CACHE_DIR}', True)
cache_key = 'salmon_quant_${sample}_\${input_hash}'
input_files = ['${fastq1}']
${fastq2 && fastq2 != 'null' ? "input_files.append('${fastq2}')" : ""}
input_files.append('${index_dir}')
params = {'libtype': '${params.salmon?.libtype ?: 'A'}', 'extra': '${params.salmon?.extra ?: ''}', 'threads': ${task.cpus}}

sample_output_dir = '${params.paths.salmon}/${sample}'
output_files = ['${sample_output_dir}/quant.sf', '${sample_output_dir}/lib_format_counts.json']
cache_mgr.store_cache(cache_key, cache_mgr.compute_input_hash([__import__('pathlib').Path(f) for f in input_files], params), [__import__('pathlib').Path(f) for f in output_files])
                "
            fi
        else
            # Run without custom caching, rely on Nextflow cache
            mkdir -p quant
            if [ -s "${fastq2}" ] && [ "${params.se}" != "true" ]; then
              salmon quant -i ${index_dir} --libType ${params.salmon?.libtype ?: 'A'} \
                -1 ${fastq1} -2 ${fastq2} --threads ${task.cpus} ${params.salmon?.extra ?: ''} \
                -o quant
            else
              salmon quant -i ${index_dir} --libType ${params.salmon?.libtype ?: 'A'} \
                -r ${fastq1} --threads ${task.cpus} ${params.salmon?.extra ?: ''} \
                -o quant
            fi
            mv quant/quant.sf .
            mv quant/lib_format_counts.json .
        fi
        """
}
