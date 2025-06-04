#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
helpMessage = """
细菌基因组组装分析流程

使用方法:
    nextflow run main.nf --input <fastq文件路径> [选项]

必需参数:
    --input                 输入文件路径（支持逗号分隔的多个文件）

可选参数:
    --output_dir           输出目录 (默认: Output)
    --assembly_strategy    组装策略选择 (默认: flye_medaka)
                          可选值: flye_medaka, canu_racon, all

组装参数:
    --genome_size          基因组大小 (默认: 5m)
    
    Flye + Medaka参数:
    --flye_threads         Flye线程数 (默认: 12)
    --flye_iterations      Flye迭代次数 (默认: 5)
    --flye_min_overlap     Flye最小重叠长度 (默认: 1000)
    --medaka_threads       Medaka线程数 (默认: 20)
    
    Canu + Racon参数:
    --canu_threads         Canu线程数 (默认: 40)
    --minimap2_threads     Minimap2线程数 (默认: 40)
    --racon_threads        Racon线程数 (默认: 40)

质控参数:
    --nanocomp_threads     NanoComp线程数 (默认: 24)
    --fastplong_threads    Fastplong线程数 (默认: 12)
    --min_quality          最小质量值 (默认: 10)
    --min_length           最小长度 (默认: 1000)

质量评估参数:
    --quast_threads        QUAST线程数 (默认: 40)
    --busco_threads        BUSCO线程数 (默认: 40)
    --checkm_threads       CheckM线程数 (默认: 20)
    --pplacer_threads      Pplacer线程数 (默认: 8)

Nextflow常用参数:
    -profile              配置选择
    -resume               从上次运行中断的地方继续
    -with-trace          生成执行跟踪报告
    -with-timeline       生成执行时间线报告
    -with-dag           生成DAG图（需要Graphviz）
    -with-report        生成HTML格式的执行报告
    -log                指定日志文件
    -c                  指定配置文件
    -bg                 在后台运行


示例:
    nextflow run main.nf --input "sample1.fastq.gz"
    nextflow run main.nf --input "sample1.fastq.gz" --assembly_strategy all
    nextflow run main.nf --input "sample1.fastq.gz" --genome_size 4m --flye_threads 24 
    
    从断点继续运行:
    nextflow run main.nf -resume --input "sample1.fastq.gz"
    
    生成执行报告:
    nextflow run main.nf --input "sample1.fastq.gz" -with-trace -with-timeline -with-dag
""".stripIndent()

if (params.help) {
        log.info helpMessage
        exit 0
    }


process nanocomp {
    publishDir "${params.output_dir}/1_Nanocomp/", mode:"copy",saveAs: {file -> file.replace("Nanocomp/","")}

    input:
    val sample_ids
    path fastq_files

    output:
    path "Nanocomp/*"

    script:
    """
    mkdir -p Nanocomp
    NanoComp --fastq ${fastq_files.join(" ")} \
    --names ${sample_ids.join(" ")} \
    -t ${params.nanocomp_threads} ${params.nanocomp_params} \
    -o Nanocomp
    """
}

process fastplong {
    publishDir "${params.output_dir}/2_Fastplong/${sample_id}", mode:"copy",saveAs: {file -> file.replace("fastplong/","")}

    input:
    tuple val(sample_id),path(fastq_file)

    output:
    path "fastplong/*"
    path "fastplong/${sample_id}.clean.fastq.gz",emit:clean_fastq

    script:
    """
    mkdir -p fastplong
    fastplong -i ${fastq_file} \
    -o fastplong/${sample_id}.clean.fastq.gz \
    -h fastplong/${sample_id}.report.html \
    -j fastplong/${sample_id}.report.json \
    -m ${params.min_quality} -l ${params.min_length} \
    ${params.fastplong_params} -w ${params.fastplong_threads}
    """
}

process flye_medaka{
    publishDir "${params.output_dir}/3_Flye_Medaka/${sample_id}", mode:"copy"

    input:
    tuple val(sample_id),path(fastq_file)

    output:
    path "flye/*"
    path "medaka/*"
    path "medaka/${sample_id}.flye.medaka.assembly.fa",emit:assembly_file

    script:
    """
    mkdir -p flye
    mkdir -p medaka
    if file -L ${fastq_file}|grep -q "gzip compressed data";then
        gunzip -c ${fastq_file} > ${sample_id}.clean.fastq
    else
        ln -s ${fastq_file} ${sample_id}.clean.fastq
    fi

    flye --nano-raw ${sample_id}.clean.fastq --out-dir flye/ \
    --genome-size ${params.genome_size} \
    --iterations ${params.flye_iterations} \
    --min-overlap ${params.flye_min_overlap} \
    ${params.flye_params} --threads ${params.flye_threads} 

    medaka_consensus -i ${sample_id}.clean.fastq \
    -d flye/assembly.fasta \
    -o medaka/ \
    ${params.medaka_params} -t ${params.medaka_threads} \
    --bacteria 
    mv medaka/consensus.fasta medaka/${sample_id}.flye.medaka.assembly.fa
    """
}

process canu_racon {
    publishDir "${params.output_dir}/4_Canu_Racon/${sample_id}", mode:"copy"

    input:
    tuple val(sample_id),path(fastq_file)

    output:
    path "canu_racon/*"
    path "canu_racon/${sample_id}.canu.racon.assembly.fa",emit:assembly_file

    
    script:
    """
    canu -p ${sample_id} -d canu \
    genomeSize=${params.genome_size} \
    ${params.canu_params} maxThreads=${params.canu_threads} \
    -nanopore ${fastq_file}

    minimap2 -d canu/${sample_id}_ref.mmi canu/${sample_id}.contigs.fasta -t ${params.minimap2_threads}
    minimap2 -ax map-ont canu/${sample_id}_ref.mmi ${fastq_file} > ${sample_id}.sam

    mkdir -p canu_racon
    racon -t ${params.racon_threads} ${fastq_file} \
    ${sample_id}.sam ${params.racon_params} \
    canu/${sample_id}.contigs.fasta > canu_racon/${sample_id}.canu.racon.assembly.fa
    mv canu/${sample_id}.contigs.fasta canu_racon/${sample_id}.canu.contigs.fasta
    """
}

process quast {
    publishDir "${params.output_dir}/5_Quast/${sample_id}", mode:"copy",saveAs: {file -> file.replace("quast/","")}

    input:
    tuple val(sample_id),val(assembly_strategy),path(assembly_file)
    
    output:
    path "quast/*"

    script:
    """
    quast.py -o quast/${assembly_strategy} \
    -t ${params.quast_threads} ${params.quast_params} \
    ${assembly_file}
    """
}

process busco {
    publishDir "${params.output_dir}/6_Busco/${sample_id}", mode:"copy",saveAs: {file -> file.replace("busco/","")}

    input:
    tuple val(sample_id),val(assembly_strategy),path(assembly_file)

    output:
    path "busco/*"

    script:
    """
    busco -i ${assembly_file} \
    -l ${params.busco_odb_dir} \
    -o busco/${assembly_strategy} \
    -m genome ${params.busco_params} \
    -c ${params.busco_threads} \
    """
}

process checkm {
    publishDir "${params.output_dir}/7_Checkm/${sample_id}", mode:"copy",saveAs: {file -> file.replace("checkm/","")}

    input:
    tuple val(sample_id),val(assembly_strategy),path(assembly_file)

    output:
    path "checkm/*"

    script:
    """
    mkdir -p checkm
    checkm data setRoot ${params.checkm_databases}
    checkm lineage_wf ./ checkm/ --extension fa \
    -t ${params.checkm_threads} ${params.checkm_params} \
    --pplacer_threads ${params.pplacer_threads} \
    -f checkm/${sample_id}_${assembly_strategy}.checkm.txt \
    --tmpdir checkm
    """
}

workflow {

    if(params.input == null || params.input == "") {
        error "param input is required"
    }

    if(!(params.assembly_strategy in ["flye_medaka", "canu_racon","all"])) {
        error "Invalid assembly_strategy: '${params.assembly_strategy}'. Valid options: 'flye_medaka', 'canu_racon', 'all' "
    }

    Channel
        .from(params.input.split(','))
        .map { path -> 
            def fastq_file = file(path)
            def sampleName = fastq_file.getName().replaceAll('(\\_|\\.)(fastq|fq)(\\.gz)?$', '')
            return tuple(sampleName, fastq_file)
        }
        .set { input_ch }
    all_samples = input_ch.map { it[0] }.collect()
    all_fastq_files = input_ch.map { it[1] }.collect()

    nanocomp(all_samples, all_fastq_files)
    fastplong(input_ch)
    clean_fastq_ch = input_ch.combine(fastplong.out.clean_fastq).map {sample_id,raw_fastq,clean_fastq -> tuple(sample_id,clean_fastq)}

    flye_medaka_ch = Channel.of()
    canu_racon_ch = Channel.of()
    
    if (params.assembly_strategy in ["flye_medaka", "all"]) {
        flye_medaka(clean_fastq_ch)
        flye_medaka_ch = clean_fastq_ch
                            .combine(flye_medaka.out.assembly_file)
                            .map {sample_id,fastq_file,assembly_file ->
                                tuple(sample_id, "flye_medaka", assembly_file)
                            }
    }
    
    if (params.assembly_strategy in ["canu_racon", "all"]) {
        canu_racon(clean_fastq_ch)
        canu_racon_ch = clean_fastq_ch
                            .combine(canu_racon.out.assembly_file)
                            .map {sample_id,fastq_file,assembly_file ->
                                tuple(sample_id, "canu_racon", assembly_file)
                            }
    }
    assembly_ch = Channel.of().mix(flye_medaka_ch).mix(canu_racon_ch)
    quast(assembly_ch)
    busco(assembly_ch)
    checkm(assembly_ch)
}

