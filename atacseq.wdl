version 1.0

task process_atacseq_sample {
    input{
        String output_dir
        String sample_name
        String sample_dir = "~{output_dir}/~{sample_name}"
        String sample_raws
        String bowtie2_ref

        Int threads = 4
    }
    command {
        mkdir -p ~{output_dir};
        mkdir -p ~{sample_dir};
        mkdir -p "~{sample_dir}/tmp";

        for i in "~{sample_raws}"; do samtools fastq $i 2> "~{sample_dir}/samtools_fastq.log" ; done | \
            fastp --stdin --interleaved_in --stdout 2> "~{sample_dir}/~{sample_name}_fastp.log" | \
            bowtie2 --very-sensitive --no-discordant -p ~{threads} --maxins 2000 -x ~{bowtie2_ref} --met-file "~{sample_dir}/~{sample_name}_bowtie2.met" --interleaved - 2> "~{sample_dir}/~{sample_name}_bowtie.log" | \
            samtools fixmate -O SAM - - 2> "~{sample_dir}/~{sample_name}_fixmate.log" | \
            samblaster 2> "~{sample_dir}/~{sample_name}_samblaster.log" | \
            samtools sort -o "~{sample_dir}/~{sample_name}.bam" - 2> "~{sample_dir}/~{sample_name}_samtools_sort.log";
            samtools index "~{sample_dir}/~{sample_name}.bam" 2> "~{sample_dir}/~{sample_name}_samtools_index.log";
    }
    output {
        File output_bam = "~{sample_dir}/~{sample_name}.bam"
        File output_bai = "~{sample_dir}/~{sample_name}.bam.bai"
    }
}

workflow atacseq {
    input{
        String output_dir
        Array[Map[String, String]] sample_list
        Map[String, String] tools
        String bowtie2_ref
    }

    scatter(sample in sample_list) {
        call process_atacseq_sample {
            input:
                output_dir = "~{output_dir}",
                sample_name = "~{sample[sample_name]}",
                sample_raws = "~{sample[sample_raws]}",
                bowtie2_ref = "~{bowtie2_ref}"
        }
    }

}