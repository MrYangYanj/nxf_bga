# 细菌基因组组装分析流程

该Nextflow流程基于三代纳米孔测序数据（序风测序），用于细菌完成图基因组组装和评估，支持两种组装策略和质量评估。

## 功能特点

- 原始数据质量评估 (NanoComp)
- 数据质控 (Fastplong)
- 组装策略：
  - Flye + Medaka
  - Canu + Racon
- 组装质量评估：
  - QUAST
  - BUSCO
  - CheckM

## 安装要求

- Nextflow (>= 20.04.0)
- Docker 且用户需要有docker组权限

## 使用方法

### 基本运行

```bash
nextflow run main.nf --input path/to/reads.fastq.gz
```

### 主要参数

- `--input`: 输入文件路径（必需，支持逗号分隔的多个文件）
- `--output_dir`: 输出目录（默认：Output）
- `--assembly_strategy`: 组装策略选择，有以下三种
  - flye_medaka（默认）
  - canu_racon
  - all（同时运行两种策略）

### 详细参数配置

#### 组装参数
- `--genome_size`: 基因组大小（默认：5m）

#### Flye + Medaka参数
- `--flye_threads`: Flye线程数（默认：20）
- `--flye_iterations`: 迭代次数（默认：5）
- `--flye_min_overlap`: 最小重叠长度（默认：1000）
- `--medaka_threads`: Medaka线程数（默认：20）
- `--flye_params`: Flye额外参数
- `--medaka_params`: Medaka额外参数

#### Canu + Racon参数
- `--canu_threads`: Canu线程数（默认：40）
- `--minimap2_threads`: Minimap2线程数（默认：40）
- `--racon_threads`: Racon线程数（默认：40）
- `--canu_params`: Canu额外参数
- `--racon_params`: Racon额外参数

#### 质控参数
- `--nanocomp_threads`: NanoComp线程数（默认：24）
- `--nanocomp_params`: NanoComp额外参数（默认：-f pdf --plot violin）
- `--fastplong_threads`: Fastplong线程数（默认：12）
- `--min_quality`: 最小质量值（默认：10）
- `--min_length`: 最小长度（默认：1000）
- `--fastplong_params`: Fastplong额外参数

#### 质量评估参数
- `--quast_threads`: QUAST线程数（默认：40）
- `--quast_params`: QUAST额外参数
- `--busco_threads`: BUSCO线程数（默认：40）
- `--busco_odb_dir`: BUSCO数据库路径（默认：${projectDir}/databases/bacteria_odb12）
- `--busco_params`: BUSCO额外参数
- `--checkm_threads`: CheckM线程数（默认：20）
- `--pplacer_threads`: Pplacer线程数（默认：8）
- `--checkm_databases`: CheckM数据库路径（默认：${projectDir}/databases/checkm_data）
- `--checkm_params`: CheckM额外参数（默认：--tab_table）

### Nextflow常用参数
- `-profile`: 配置选择，目前没有设置不同的配置
- `-resume`: 从上次运行中断的地方继续
- `-with-trace`: 生成执行跟踪报告
- `-with-timeline`: 生成执行时间线报告
- `-with-dag`: 生成DAG图（需要Graphviz）
- `-with-report`: 生成HTML格式的执行报告
- `-log`: 指定日志文件
- `-c`: 指定配置文件
- `-bg`: 在后台运行

## 输出目录结构

```
Output/
├── 1_Nanocomp/          # 原始数据质量评估
├── 2_Fastplong/         # 数据质控结果
├── 3_Flye_Medaka/       # Flye+Medaka组装结果
├── 4_Canu_Racon/        # Canu+Racon组装结果
├── 5_Quast/            # QUAST评估结果
├── 6_Busco/            # BUSCO评估结果
└── 7_Checkm/           # CheckM评估结果
```

## 示例命令

1. 基本运行：
```bash
nextflow run main.nf --input sample1.fastq.gz
```

2. 多样本运行：
```bash
nextflow run main.nf --input sample1.fastq.gz,sample2.fastq.gz
```

3. 使用Flye+Medaka策略：
```bash
nextflow run main.nf --input sample1.fastq.gz --assembly_strategy flye_medaka
```

4. 使用Canu+Racon策略：
```bash
nextflow run main.nf --input sample1.fastq.gz --assembly_strategy canu_racon
```

5. 同时使用两种策略：
```bash
nextflow run main.nf --input sample1.fastq.gz --assembly_strategy all
```

6. 自定义参数：
```bash
nextflow run main.nf \
  --input sample1.fastq.gz \
  --genome_size 4m \
  --flye_threads 24 \
  --medaka_threads 32 \
  --output_dir "results"
```

7. 指定自定义配置文件：
```bash
nextflow run main.nf \
  --input sample1.fastq.gz \
  --genome_size 4m \
  --output_dir "results" \
  -c custom_nextflow.config
```

8. 从断点继续运行：
```bash
nextflow run main.nf -resume --input sample1.fastq.gz
```

9. 生成执行报告：
```bash
nextflow run main.nf --input sample1.fastq.gz -with-trace -with-timeline -with-dag
```

## 注意事项

1. 当输入文件指定了多个样本时，除NanoComp外其余流程会并行处理，请相应调节单样本运行线程
2. 根据实际数据量以及硬件配置调整线程数和内存
3. 确保有足够的磁盘空间存储中间文件和结果
4. 建议使用绝对路径指定输入文件
5. 程序运行中的样本名约定为fastq文件的前缀，请尽量不要使用特殊符号
6. 运行完成后，如不需要再进行断点继续运行，请清理nextflow生成的work文件夹（命令:nextflow clean -f 或rm -rf ./work）