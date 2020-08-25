version 1.0

struct Sample {
    String sample_name
    String read_type
    String genome
    Int genome_size
    String raw_bams
    Int raw_size_mb
}


task test_task {
    input{
        String output_dir
        String sample_name
    }
    String sample_dir = "~{output_dir}/~{sample_name}"
    File sample_tsv  = "~{output_dir}/config_files/~{sample_name}.tsv"
    Map[String, String] sample_map = read_map(sample_tsv)
    String sample_raws = sample_map['raw_bams']

    command{
        echo "~{sample_dir}"
        echo "~{sample_raws}"
    }
    output{
        String my_out = stdout()
    }
}

task bowtie2_align {
    input{
        String output_dir
        String config_dir
        String sample_name
        String bowtie2_index
        String? adapter_fasta

        Int cpus = 8
    }

    String sample_dir = "~{output_dir}/~{sample_name}"
    File sample_tsv  = "~{config_dir}/~{sample_name}.tsv"
    Map[String, String] sample_map = read_map(sample_tsv)
    String raw_bams = sample_map['raw_bams']
    Int raw_size_mb = sample_map['raw_size_mb']
    Int memory = 6000 + raw_size_mb
    String interleaved_in = if(sample_map['read_type'] == 'paired') then '--interleaved_in' else ' '
    String interleaved = if(sample_map['read_type'] == 'paired') then '--interleaved' else ' '

    command {
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{sample_dir}" ] && mkdir -p ~{sample_dir};

        for i in ~{raw_bams}; do samtools fastq $i 2>> "~{sample_dir}/~{sample_name}.samtools.log" ; done | \
            fastp --stdin ~{interleaved_in} --stdout --html "~{sample_dir}/~{sample_name}.fastp.html" --json "~{sample_dir}/~{sample_name}.fastp.json" 2> "~{sample_dir}/~{sample_name}.fastp.log" | \
            bowtie2 --very-sensitive --no-discordant -p ~{cpus} --maxins 2000 -x ~{bowtie2_index} --met-file "~{sample_dir}/~{sample_name}.bowtie2.met" ~{interleaved} - 2> "~{sample_dir}/~{sample_name}.txt" | \
            samblaster --addMateTags 2> "~{sample_dir}/~{sample_name}.samblaster.log" | \
            samtools sort -o "~{sample_dir}/~{sample_name}.bam" - 2>> "~{sample_dir}/~{sample_name}.samtools.log";
            samtools index "~{sample_dir}/~{sample_name}.bam" 2>> "~{sample_dir}/~{sample_name}.samtools.log";
            samtools flagstat "~{sample_dir}/~{sample_name}.bam" > "~{sample_dir}/~{sample_name}.samtools_flagstat.log";
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: "mediumq"
        rt_time: "2-00:00:00"
    }

    output {
        File output_bam = "~{sample_dir}/~{sample_name}.bam"
        File output_bai = "~{sample_dir}/~{sample_name}.bam.bai"
    }
}

workflow atacseq {
    input{
        String project_path
        String project_name
        String genome
        String regulatory_regions
        String blacklisted_regions
        String chromosome_sizes
        String unique_tss
        Array[String] sample_list
        String bowtie2_index
        String adapter_fasta
    }

    String output_dir = "~{project_path}/atacseq_results"
    String config_dir = "~{project_path}/config_files"

    scatter(sample in sample_list) {

        call bowtie2_align {
            input:
                output_dir = output_dir,
                config_dir = config_dir,
                sample_name = sample,
                bowtie2_index = bowtie2_index,
                adapter_fasta = adapter_fasta
        }
    }
}