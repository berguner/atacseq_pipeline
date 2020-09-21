# ATAC-seq Data Processing Pipeline

This pipeline is developed for processing raw/unaligned ATAC-seq data to generate peak calls, a QC report and motif 
discovery results. The pipeline can process both single-end and paired-end reads in unaligned BAM format.

## 1 - Installation

1. Clone this repository.
2. Install the software listed in the `required_software.txt`.
3. Install atacseq plugin for MultiQC: `python3 /path/to/atacseq_pipeline/multiqc_atacseq/setup.py install`
4. Collect all the reference data listed in the `atacseq_config.yaml` file and update the file with correct paths. This 
file contains the pipeline configurations common among different projects.

## 2 - Running the pipeline

### 2.1 - Setting up new analysis/project

For each project you need to prepare a project config (`.yaml`) file and a sample annotations sheet (`.csv`).
You can set the project name, output path, genome version etc. in the project config file.
Sample annotations sheet should contain one line for each input raw BAM file, multiple files with the same sample_name
will be merged during the alignment step. Make sure that the `data_source` template of each file is configured in the 
config file under the `data_sources:` attribute. These templates will be used to generate the full path of the input files.

### 2.2 - Running the pipeline
1. Generate the input files for the `cromwell` WDL engine:
```
python3 /path/to/atacseq_pipeline/configurator.py \
    -p /path/to/atacseq_pipeline/atacseq_config.yaml \
    -c /path/to/my_project_config.yaml
``` 
2. Go to the project output path and run `cromwell` workflow manager:
```
cd /my/project/path
java -Dconfig.file=/path/to/atacseq_pipeline/backends/slurm.conf -Xmx2g \
    -jar /path/to/cromwell-52.jar run /path/to/atacseq_pipeline/atacseq.wdl \
    --inputs config_files/BSA_0000_test_atac.inputs.json
```
3. After the pipeline has finished run the MultiQC to generate the QC report:
```
multiqc -c /path/to/my_project_config.yaml /my/project/path
```