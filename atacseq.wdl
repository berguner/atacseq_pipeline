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
        String whitelist
        String chrM

        Int cpus = 8
    }

    String sample_dir = "~{output_dir}/~{sample.sample_name}"
    String bam_dir = "~{output_dir}/~{sample.sample_name}/mapped"
    String raw_bams = sample.raw_bams

    String interleaved_in = if(sample.read_type == "paired") then "--interleaved_in" else " "
    String interleaved = if(sample.read_type == "paired") then "--interleaved" else " "
    String filter = if(sample.read_type == "paired") then "-q 30 -F 2316 -f 2 -L ~{whitelist}" else "-q 30 -F 2316 -L ~{whitelist}"

    #Int raw_size_mb = sample.raw_size_mb #sample_map['raw_size_mb']
    Int memory = 16000

    command <<<
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{sample_dir}" ] && mkdir -p ~{sample_dir};
        [ ! -d "~{bam_dir}" ] && mkdir -p ~{bam_dir};

        for i in ~{raw_bams}; do samtools fastq $i 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" ; done | \
            fastp --stdin ~{interleaved_in} --stdout --html "~{bam_dir}/~{sample.sample_name}.fastp.html" --json "~{bam_dir}/~{sample.sample_name}.fastp.json" 2> "~{bam_dir}/~{sample.sample_name}.fastp.log" | \
            bowtie2 --very-sensitive --no-discordant -p ~{cpus} --maxins 2000 -x ~{bowtie2_index} --met-file "~{bam_dir}/~{sample.sample_name}.bowtie2.met" ~{interleaved} - 2> "~{bam_dir}/~{sample.sample_name}.txt" | \
            samblaster --addMateTags 2> "~{bam_dir}/~{sample.sample_name}.samblaster.log" | \
            samtools sort -o "~{bam_dir}/~{sample.sample_name}.bam" - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log";

        samtools index "~{bam_dir}/~{sample.sample_name}.bam" 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log";
        samtools idxstats "~{bam_dir}/~{sample.sample_name}.bam" | awk '{ sum += $3 + $4; if($1 == "~{chrM}") { mito_count = $3; }}END{ print "mitochondrial_fraction\t"mito_count/sum }' > "~{sample_dir}/~{sample.sample_name}.stats.tsv";
        samtools flagstat "~{bam_dir}/~{sample.sample_name}.bam" > "~{bam_dir}/~{sample.sample_name}.samtools_flagstat.log";

        samtools view ~{filter} -o "~{bam_dir}/~{sample.sample_name}.filtered.bam" "~{bam_dir}/~{sample.sample_name}.bam";
        samtools index "~{bam_dir}/~{sample.sample_name}.filtered.bam";

    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: "mediumq"
        rt_time: "2-00:00:00"
    }

    output {
        File output_bam = "~{bam_dir}/~{sample.sample_name}.bam"
        File output_bai = "~{bam_dir}/~{sample.sample_name}.bam.bai"
        File filtered_bam = "~{bam_dir}/~{sample.sample_name}.filtered.bam"
        File filtered_bai = "~{bam_dir}/~{sample.sample_name}.filtered.bam.bai"
    }
}

task macs2_peak_call {
    input {
        String output_dir
        String regulatory_regions
        Sample sample

        File input_bam
        File input_bai
    }

    String sample_dir = "~{output_dir}/~{sample.sample_name}"
    String peaks_dir = "~{output_dir}/~{sample.sample_name}/peaks"
    String format = if(sample.read_type == 'paired') then '--format BAMPE' else '--format BAM'

    command <<<
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{sample_dir}" ] && mkdir -p ~{sample_dir};
        [ ! -d "~{peaks_dir}" ] && mkdir -p ~{peaks_dir};

        macs2 callpeak -t ~{input_bam} ~{format} \
            --nomodel --keep-dup auto --extsize 147 -g ~{sample.genome_size} \
            -n ~{sample.sample_name} \
            --outdir ~{peaks_dir} > "~{peaks_dir}/~{sample.sample_name}.macs2.log" 2>&1;

        cat ~{peaks_dir}/~{sample.sample_name}_peaks.narrowPeak | wc -l | \
            awk '{print "peaks\t" $1}' >> "~{sample_dir}/~{sample.sample_name}.stats.tsv"

        TOTAL_READS=`samtools idxstats ~{input_bam} | awk '{sum += $3}END{print sum}'`;
        samtools view -c -L "~{peaks_dir}/~{sample.sample_name}_peaks.narrowPeak" ~{input_bam} | \
            awk -v total=$TOTAL_READS '{print "frip\t" $1/total}' >> "~{sample_dir}/~{sample.sample_name}.stats.tsv";

        samtools view -c -L ~{regulatory_regions} ~{input_bam} | \
            awk -v total=$TOTAL_READS '{print "regulatory_fraction\t" $1/total}' >> "~{sample_dir}/~{sample.sample_name}.stats.tsv";
    >>>

    runtime {
        rt_cpus: 2
        rt_mem: 4000
        rt_queue: "shortq"
        rt_time: "12:00:00"
    }
    output {
        File peak_calls = "~{peaks_dir}/~{sample.sample_name}_peaks.narrowPeak"
        File macs2_xls = "~{peaks_dir}/~{sample.sample_name}_peaks.xls"
        File summits_bed = "~{peaks_dir}/~{sample.sample_name}_summits.bed"
    }
}


task misc_tasks {
    input {
        String project_dir
        String output_dir
        Int tss_slop = 2000
        String unique_tss
        String chromosome_sizes
        String whitelist

        Sample sample

        File input_bam
        File input_bai
    }

    String hub_dir = "~{project_dir}/atacseq_hub"
    String sample_dir = "~{output_dir}/~{sample.sample_name}"
    String slopped_tss = "~{output_dir}/~{sample.sample_name}/slopped_tss.bed"
    String tss_hist = "~{output_dir}/~{sample.sample_name}/~{sample.sample_name}.tss_histogram.csv"
    Int noise_upper = ( tss_slop * 2 ) - 100
    command <<<
        [ ! -d "~{hub_dir}" ] && mkdir -p ~{hub_dir};

        bamCoverage --bam ~{input_bam} \
            -p max --binSize 10  --normalizeUsing RPGC \
            --effectiveGenomeSize ~{sample.genome_size} --extendReads 175 \
            -o "~{hub_dir}/~{sample.sample_name}.bigWig" > "~{hub_dir}/~{sample.sample_name}.bigWig.log" 2>&1;

        echo "base,count" > ~{tss_hist};
        bedtools slop -b ~{tss_slop} -i ~{unique_tss} -g ~{chromosome_sizes} | bedtools intersect -a - -b ~{whitelist} -u > ~{slopped_tss};
        bedtools bamtobed -i ~{input_bam} | \
            bedtools shift -p -0.5 -m 0.5 -g ~{chromosome_sizes} -i - | \
            bedtools sort -g ~{chromosome_sizes} -i - | \
            bedtools coverage -a ~{slopped_tss} -b - -d -sorted | \
            awk '{counts[$7] += $8;} END { for(pos in counts){if(pos<100 || pos>~{noise_upper}){sum+=counts[pos]}}; for(pos in counts){print pos-~{tss_slop}","(counts[pos]/sum)*100}}' | \
            sort -t "," -k1,1n >> ~{tss_hist} ;
        rm ~{slopped_tss};
    >>>

    runtime {
        rt_cpus: 2
        rt_mem: 4000
        rt_queue: "shortq"
        rt_time: "12:00:00"
    }

    output {
        File bigWig = "~{hub_dir}/~{sample.sample_name}.bigWig"
    }
}

workflow atacseq {
    input{
        String project_path
        String project_name
        String genome
        String regulatory_regions
        String blacklisted_regions
        String whitelisted_regions
        String chromosome_sizes
        String unique_tss
        Array[String] sample_list
        String bowtie2_index
        String adapter_fasta
        String mitochondria_name
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
                adapter_fasta = adapter_fasta,
                whitelist = whitelisted_regions,
                chrM = mitochondria_name
        }

        call macs2_peak_call {
            input:
                output_dir = output_dir,
                sample = sample,
                input_bam = bowtie2_align.filtered_bam,
                input_bai = bowtie2_align.filtered_bai,
                regulatory_regions = regulatory_regions
        }

        call misc_tasks {
            input:
                project_dir = project_path,
                output_dir = output_dir,
                sample = sample,
                input_bam = bowtie2_align.filtered_bam,
                input_bai = bowtie2_align.filtered_bai,
                unique_tss = unique_tss,
                chromosome_sizes = chromosome_sizes,
                whitelist = whitelisted_regions
        }
    }
}