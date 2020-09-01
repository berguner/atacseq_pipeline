#!/usr/bin/env python

"""
MultiQC module to parse variant calling pipeline stats
"""

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import csv

from multiqc import config
from multiqc.plots import bargraph, linegraph, table, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    """
    variant_calling module class
    """

    def __init__(self):
        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_bsf_reports', True):
            return None
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Variant Calling Pipeline', anchor='variant_calling',
                                            info=" processes genomic data to detect variations.")
        log.info('Initialized variant_calling module')

        # Load the sample annotation sheet
        sample_sas_path = config.sample_annotation
        sample_sas = csv.DictReader(open(sample_sas_path, 'r'))
        self.sample_sas_dict = {}
        for k in sample_sas:
            self.sample_sas_dict[k['sample_name']] = k

        # Update general stats table
        self.update_general_stats()

        # Add download links table
        self.add_download_table()

    def update_general_stats(self):
        data = {}
        for sample_name in self.sample_sas_dict:
            data[sample_name] = {}
            data[sample_name]['Tissue'] = self.sample_sas_dict[sample_name]['tissue']

        headers = OrderedDict()
        headers['Tissue'] = {
            'description': 'Sample tissue',
            'title': 'Tissue',
            'scale': False
        }

        self.general_stats_addcols(data=data, headers=headers)

    def add_download_table(self):
        # Create a table with download links to various files
        results_url = config.project_url
        results_path = os.path.join(config.project_path, config.genome)

        # Configuration for the MultiQC table
        table_config = {
            'namespace': 'Download links',  # Name for grouping. Prepends desc and is in Config Columns modal
            'id': 'download_links',  # ID used for the table
            'table_title': 'Download variant calling data',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_download_links_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'col1_header': 'Sample Name',  # The header used for the first column
            'no_beeswarm': True,
            'scale': False
        }

        # Configuration for the header row
        headers = OrderedDict()
        headers['BAM'] = {
            'title': 'Aligned BAM',
            'description': 'BWA-mem alignment results in BAM format',
            'scale': False
        }
        headers['VCF'] = {
            'title': 'Germline VCF',
            'description': 'VCF files containing the germline variants detected by GATK HaplotypeCaller',
            'scale': False
        }
        headers['VEP_VCF'] = {
            'title': 'Annotated Germline VCF',
            'description': 'VCF files containing the germline variants annotated by Ensembl VEP',
            'scale': False
        }
        headers['mutect2'] = {
            'title': 'Somatic VCF',
            'description': 'VCF files containing the somatic mutations detected by GATK Mutect2',
            'scale': False
        }
        headers['VEP_mutect2'] = {
            'title': 'Annotated Somatic VCF',
            'description': 'VCF files containing the somatic mutations annotated by Ensembl VEP',
            'scale': False
        }
        headers['cnvkit'] = {
            'title': 'CNV',
            'description': 'Copy number analysis results by CNVkit',
            'scale': False
        }
        data = OrderedDict()
        for sample_name in self.sample_sas_dict:
            data[sample_name] = {}
            bam_path = os.path.join(results_path, 'bam', '{}.bam'.format(sample_name))
            bai_path = os.path.join(results_path, 'bam', '{}.bam.bai'.format(sample_name))
            if os.path.exists(bam_path) and os.path.exists(bai_path):
                sample_bam_url = '{}/{}/bam/{}.bam'.format(results_url, config.genome, sample_name)
                sample_bai_url = '{}/{}/bam/{}.bam.bai'.format(results_url, config.genome, sample_name)
                data[sample_name]['BAM'] = '<a href=\"{}\">{} bam</a><br><a href=\"{}\">{} bai</a>'.format(
                    sample_bam_url, sample_name, sample_bai_url, sample_name)

            vcf_path = os.path.join(results_path, 'vcf', '{}.vcf.gz'.format(sample_name))
            vcf_tbi_path = os.path.join(results_path, 'vcf', '{}.vcf.gz.tbi'.format(sample_name))
            if os.path.exists(vcf_path) and os.path.exists(vcf_tbi_path):
                sample_vcf_url = '{}/{}/vcf/{}.vcf.gz'.format(results_url, config.genome, sample_name)
                sample_vcf_tbi_url = '{}/{}/vcf/{}.vcf.gz.tbi'.format(results_url, config.genome, sample_name)
                data[sample_name]['VCF'] = '<a href=\"{}\">{} vcf</a><br><a href=\"{}\">{} tbi</a>'.format(sample_vcf_url, sample_name,
                                                                                           sample_vcf_tbi_url, sample_name)

            vep_vcf_path = os.path.join(results_path, 'vcf', '{}.vep.vcf.gz'.format(sample_name))
            vep_vcf_tbi_path = os.path.join(results_path, 'vcf', '{}.vep.vcf.gz.tbi'.format(sample_name))
            if os.path.exists(vep_vcf_path) and os.path.exists(vep_vcf_tbi_path):
                sample_vep_vcf_url = '{}/{}/vcf/{}.vep.vcf.gz'.format(results_url, config.genome, sample_name)
                sample_vep_vcf_tbi_url = '{}/{}/vcf/{}.vep.vcf.gz.tbi'.format(results_url, config.genome, sample_name)
                data[sample_name]['VEP_VCF'] = '<a href=\"{}\">{} vep</a><br><a href=\"{}\">{} tbi</a>'.format(
                    sample_vep_vcf_url, sample_name,
                    sample_vep_vcf_tbi_url, sample_name)

            mutect2_vcf_path = os.path.join(results_path, 'mutect2', '{}.vcf.gz'.format(sample_name))
            mutect2_vcf_tbi_path = os.path.join(results_path, 'mutect2', '{}.vcf.gz.tbi'.format(sample_name))
            if os.path.exists(mutect2_vcf_path) and os.path.exists(mutect2_vcf_tbi_path):
                sample_mutect2_vcf_url = '{}/{}/mutect2/{}.vcf.gz'.format(results_url, config.genome, sample_name)
                sample_mutect2_vcf_tbi_url = '{}/{}/mutect2/{}.vcf.gz.tbi'.format(results_url, config.genome, sample_name)
                data[sample_name]['mutect2'] = '<a href=\"{}\">{} vcf</a><br><a href=\"{}\">{} tbi</a>'.format(
                    sample_mutect2_vcf_url, sample_name,
                    sample_mutect2_vcf_tbi_url, sample_name)

            mutect2_vep_vcf_path = os.path.join(results_path, 'mutect2', '{}.vep.vcf.gz'.format(sample_name))
            mutect2_vep_vcf_tbi_path = os.path.join(results_path, 'mutect2', '{}.vep.vcf.gz.tbi'.format(sample_name))
            if os.path.exists(mutect2_vep_vcf_path) and os.path.exists(mutect2_vep_vcf_tbi_path):
                sample_mutect2_vep_vcf_url = '{}/{}/mutect2/{}.vep.vcf.gz'.format(results_url, config.genome, sample_name)
                sample_mutect2_vep_vcf_tbi_url = '{}/{}/mutect2/{}.vep.vcf.gz.tbi'.format(results_url, config.genome,
                                                                                  sample_name)
                data[sample_name]['VEP_mutect2'] = '<a href=\"{}\">{} vep</a><br><a href=\"{}\">{} tbi</a>'.format(
                    sample_mutect2_vep_vcf_url, sample_name,
                    sample_mutect2_vep_vcf_tbi_url, sample_name)

            cnvkit_path = os.path.join(results_path, 'cnvkit', 'results', '{}.cnvkit.tar.gz'.format(sample_name))
            if os.path.exists(cnvkit_path):
                sample_cnvkit_url = '{}/{}/cnvkit/results/{}.cnvkit.tar.gz'.format(results_url, config.genome, sample_name)
                data[sample_name]['cnvkit'] = '<a href=\"{}\">{} CNV</a>'.format(
                    sample_cnvkit_url, sample_name)
        section_description = 'Download links for the variant calling results'

        cohort_vcf_path = os.path.join(results_path, 'vcf', '{}_cohort.vcf.gz'.format(config.project_name))
        cohort_vcf_url = '{}/{}/vcf/{}_cohort.vcf.gz'.format(results_url, config.genome, config.project_name)
        cohort_vcf_tbi_url = '{}/{}/vcf/{}_cohort.vcf.gz.tbi'.format(results_url, config.genome, config.project_name)
        cohort_vcf_vep_url = '{}/{}/vcf/{}_cohort.vep.vcf.gz'.format(results_url, config.genome, config.project_name)
        cohort_vcf_vep_tbi_url = '{}/{}/vcf/{}_cohort.vep.vcf.gz.tbi'.format(results_url, config.genome, config.project_name)
        if os.path.exists(cohort_vcf_path):
            section_description += '<br><a href=\"{}\" target=\"_blank\">' \
                                   'Click here to download the cohort (multi sample) germline VCF ' \
                                   '<span class="glyphicon glyphicon-new-window"></span></a>'.format(cohort_vcf_url)
            section_description += '<a href=\"{}\" target=\"_blank\">' \
                                   ' VCF index ' \
                                   '<span class="glyphicon glyphicon-new-window"></span></a>'.format(cohort_vcf_tbi_url)
            section_description += '<br><a href=\"{}\" target=\"_blank\">' \
                                   'Click here to download the cohort (multi sample) VEP annotated germline VCF ' \
                                   '<span class="glyphicon glyphicon-new-window"></span></a>'.format(cohort_vcf_vep_url)
            section_description += '<a href=\"{}\" target=\"_blank\">' \
                                   ' VCF index ' \
                                   '<span class="glyphicon glyphicon-new-window"></span></a>'.format(cohort_vcf_vep_tbi_url)

        # Finally add a MultiQC section together with the URL table
        self.add_section(
            name='Download Links',
            anchor='variant_calling_download',
            description=section_description,
            helptext='You can click on the table elements to download the files',
            plot=table.plot(data, headers, table_config)
        )