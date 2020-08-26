version 1.0

struct Sample {
    String sample_name
    String read_type
    String genome
    String genome_size
    String raw_bams
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
        Sample sample
        String bowtie2_index
        String? adapter_fasta

        Int cpus = 8
    }

    String sample_dir = "~{output_dir}/~{sample.sample_name}"
    String bam_dir = "~{output_dir}/~{sample.sample_name}/bam"
    String raw_bams = sample.raw_bams

    String interleaved_in = if(sample.read_type == 'paired') then '--interleaved_in' else ' '
    String interleaved = if(sample.read_type == 'paired') then '--interleaved' else ' '

    #Int raw_size_mb = sample.raw_size_mb #sample_map['raw_size_mb']
    Int memory = 16000

    command {
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{sample_dir}" ] && mkdir -p ~{sample_dir};
        [ ! -d "~{bam_dir}" ] && mkdir -p ~{bam_dir};

        for i in ~{raw_bams}; do samtools fastq $i 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" ; done | \
            fastp --stdin ~{interleaved_in} --stdout --html "~{bam_dir}/~{sample.sample_name}.fastp.html" --json "~{bam_dir}/~{sample.sample_name}.fastp.json" 2> "~{bam_dir}/~{sample.sample_name}.fastp.log" | \
            bowtie2 --very-sensitive --no-discordant -p ~{cpus} --maxins 2000 -x ~{bowtie2_index} --met-file "~{bam_dir}/~{sample.sample_name}.bowtie2.met" ~{interleaved} - 2> "~{bam_dir}/~{sample.sample_name}.txt" | \
            samblaster --addMateTags 2> "~{bam_dir}/~{sample.sample_name}.samblaster.log" | \
            samtools sort -o "~{bam_dir}/~{sample.sample_name}.bam" - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log";
            samtools index "~{bam_dir}/~{sample.sample_name}.bam" 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log";
        #if [ $? -eq "0" ]; then
            samtools flagstat "~{bam_dir}/~{sample.sample_name}.bam" > "~{bam_dir}/~{sample.sample_name}.samtools_flagstat.log";
        #else
        #    exit $?;
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: "mediumq"
        rt_time: "2-00:00:00"
    }

    output {
        File output_bam = "~{bam_dir}/~{sample.sample_name}.bam"
        File output_bai = "~{bam_dir}/~{sample.sample_name}.bam.bai"
        File flagstat = "~{bam_dir}/~{sample.sample_name}.samtools_flagstat.log"
    }
}

task macs2_peak_call {
    input {
        String output_dir
        Sample sample

        File input_bam
        File input_bai
    }

    String sample_dir = "~{output_dir}/~{sample.sample_name}"
    String peaks_dir = "~{output_dir}/~{sample.sample_name}/peaks"

    command {
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{sample_dir}" ] && mkdir -p ~{sample_dir};
        [ ! -d "~{peaks_dir}" ] && mkdir -p ~{peaks_dir};

        macs2 callpeak -t ~{input_bam} \
            --nomodel --extsize 147 -g ~{sample.genome_size} \
            -n ~{sample.sample_name} \
            --outdir ~{peaks_dir} > "~{peaks_dir}/~{sample.sample_name}.macs2.log" 2>&1;

    }

    runtime {
        rt_cpus: 2
        rt_mem: 8000
        rt_queue: "shortq"
        rt_time: "12:00:00"
    }
    output {
        File peak_calls = "~{peaks_dir}/~{sample.sample_name}.narrowPeak"
        File macs2_xls = "~{peaks_dir}/~{sample.sample_name}.xls"
        File summits_bed = "~{peaks_dir}/~{sample.sample_name}_summits.bed"
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

    scatter(sample_name in sample_list) {
        File sample_tsv  = "~{config_dir}/~{sample_name}.tsv"
        Map[String, String] sample_map = read_map(sample_tsv)
        Sample sample = { "sample_name": sample_name,
                            "read_type": sample_map["read_type"],
                            "raw_bams": sample_map["raw_bams"],
                            "genome": sample_map["genome"], "genome_size": sample_map["genome_size"]
                        }

        call bowtie2_align {
            input:
                output_dir = output_dir,
                sample = sample,
                bowtie2_index = bowtie2_index,
                adapter_fasta = adapter_fasta
        }
        call macs2_peak_call {
            input:
                output_dir = output_dir,
                sample = sample,
                input_bam = bowtie2_align.output_bam,
                input_bai = bowtie2_align.output_bai
        }
    }
}