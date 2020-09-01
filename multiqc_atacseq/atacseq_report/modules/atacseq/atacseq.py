#!/usr/bin/env python

"""
MultiQC module to parse ATAC-seq pipeline stats
"""

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
import os
import csv
import numpy as np

from multiqc import config
from multiqc.plots import linegraph, table, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):
    """
    atacseq module class
    """

    def __init__(self):
        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_atacseq_report', True):
            return None

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='ATAC-seq Pipeline', anchor='atacseq',
                                            href='https://github.com/berguner/atacseq_pipeline',
                                            info="The ATAC-seq pipeline processes ATAC-seq data.")
        log.info('Initialized atacseq module')
        # Parse ATAC-seq stats for each sample
        self.atacseq_data = dict()
        for f in self.find_log_files(sp_key='atacseq'):
            self.atacseq_data[f['s_name']] = self.parse_atacseq_stats(f['f'])
        log.info('Found stats file for {} ATAC-seq samples'.format(len(self.atacseq_data)))

        # Raise the not found warning
        if len(self.atacseq_data) == 0:
            raise UserWarning

        # Remove ignored samples if there is any
        self.atacseq_data = self.ignore_samples(self.atacseq_data)

        # Parse TSS for each sample
        self.atacseq_tss_data = dict()
        for f in self.find_log_files(sp_key='atacseq/tss'):
            self.atacseq_tss_data[f['s_name']] = self.parse_atacseq_tss(f['f'])
        log.info('Found TSS file for {} ATAC-seq samples'.format(len(self.atacseq_tss_data)))
        # Remove ignored samples if there is any
        self.atacseq_tss_data = self.ignore_samples(self.atacseq_tss_data)

        # Load the sample annotation sheet
        sample_sas_path = config.sample_annotation
        sample_sas = csv.DictReader(open(sample_sas_path, 'r'))
        self.sample_sas_dict = {}
        for k in sample_sas:
            self.sample_sas_dict[k['sample_name']] = k

        # Check if there are any paired end sample in the current project
        self.pairedSampleExists = False
        self.organism = ''
        for sample in self.sample_sas_dict:
            if self.sample_sas_dict[sample]['read_type'] == 'paired':
                self.pairedSampleExists = True
            self.organism = self.sample_sas_dict[sample]['organism']

        # # Set the color palette for PCA plot
        # self.attribute_colors = {}
        # self.color_attribute = config.pca_color_attribute
        # for s in self.sample_sas_dict:
        #     if self.sample_sas_dict[s][self.color_attribute] not in self.attribute_colors:
        #         self.attribute_colors[self.sample_sas_dict[s][self.color_attribute]] = '#%02x%02x%02x' % tuple(
        #             np.random.choice(range(256), size=3))
        # Get the genome version
        self.genome_version = config.genome

        # Load global statistics
        self.global_stats = pd.read_hdf(config.atacseq_global_hdf)

        # Add stats to general table
        self.add_atacseq_to_general_stats()

        # Add download links table
        self.add_download_table()

        # # Add FRiP scatter plot
        # self.add_frip_plot()

        # Add TSS line graph
        self.add_tss_plot()

        # # Add PCA plots
        # self.pca_dict = {}
        # self.add_pca_plots()

    def parse_atacseq_stats(self, f):
        data = {}
        for l in f.splitlines():
            s = l.split('\t')
            data[s[0]] = s[1]
        return data

    def parse_atacseq_tss(self, f):
        data = OrderedDict()
        count = 0
        for l in f.splitlines():
            s = l.split(',')
            if s[0] == 'base':
                continue
            count += 1
            if count % 10 == 0:
                data[int(s[0])] = float(s[1])
        return data

    def add_atacseq_to_general_stats(self):
        data = {}
        global_df = pd.DataFrame(self.global_stats, columns=['sample_name', 'NSC', 'RSC'])
        global_df['NSC_Percentile'] = global_df.NSC.rank(pct=True, na_option='keep')
        global_df['RSC_Percentile'] = global_df.RSC.rank(pct=True, na_option='keep')
        global_df.set_index('sample_name', inplace=True)
        global_df['NSC'] = global_df['NSC'].astype(float)
        global_df['RSC'] = global_df['RSC'].astype(float)
        #global_df = global_df.astype('float')
        for sample_name in self.atacseq_data:
            data[sample_name] = {}
            if hasattr(self, 'color_attribute') and self.color_attribute in self.sample_sas_dict[sample_name]:
                data[sample_name][self.color_attribute] = self.sample_sas_dict[sample_name][self.color_attribute]
            if 'NSC' in self.atacseq_data[sample_name] and self.atacseq_data[sample_name]['NSC'] != 'nan':
                try:
                    value = float(self.atacseq_data[sample_name]['NSC'])
                except ValueError as err:
                    print(err)
                    value = 'NaN'
                try:
                    pct_value = global_df.loc[(global_df['NSC'] - value).abs().idxmin()]['NSC_Percentile'] * 100
                except ValueError as err:
                    print(err)
                    pct_value = 'NaN'
                data[sample_name]['NSC'] = value
                data[sample_name]['NSC_PCT'] = pct_value
            if 'RSC' in self.atacseq_data[sample_name] and self.atacseq_data[sample_name]['RSC'] != 'nan':
                try:
                    value = float(self.atacseq_data[sample_name]['RSC'])
                except:
                    value = 'NaN'
                try:
                    pct_value = global_df.loc[(global_df['RSC'] - value).abs().idxmin()]['RSC_Percentile'] * 100
                except ValueError as err:
                    print(err)
                    pct_value = 'NaN'
                data[sample_name]['RSC'] = value
                data[sample_name]['RSC_PCT'] = pct_value
            if 'peaks' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['peaks'])
                except:
                    value = None
                data[sample_name]['peaks'] = value
            if 'filtered_peaks' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['filtered_peaks'])
                except:
                    value = None
                data[sample_name]['filtered_peaks'] = value
            if 'frip' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['frip'])
                except:
                    value = None
                data[sample_name]['frip'] = value
            if 'regulatory_frip' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['regulatory_frip'])
                except:
                    value = None
                data[sample_name]['regulatory_fraction'] = value
        headers = OrderedDict()
        if hasattr(self, "color_attribute"):
            headers[self.color_attribute] = {
                'description': self.color_attribute,
                'title': self.color_attribute,
                'scale': False
            }
        headers['peaks'] = {
            'description': 'Number of detected peaks',
            'title': 'Peaks',
            'scale': 'Greens',
            'format': '{:,.0f}'
        }
        headers['filtered_peaks'] = {
            'description': 'Number of peaks remaining after filtering',
            'title': 'Filtered\nPeaks',
            'scale': 'Greens',
            'format': '{:,.0f}'
        }
        headers['NSC'] = {
            'description': 'Normalized Strand Cross-correlation Coefficient',
            'title': 'NSC',
            'scale': 'Reds',
            'min': 0.0,
            'max': 2.0,
            'format': '{:.2f}'
        }
        headers['NSC_PCT'] = {
            'description': 'NSC Percentile Among All ATAC-seq samples',
            'title': 'NSC_PCT',
            'scale': 'Reds',
            'suffix': '%',
            'max': 100,
            'format': '{:,.0f}'
        }
        headers['RSC'] = {
            'description': 'Relative Strand Cross-correlation Coefficient',
            'title': 'RSC',
            'scale': 'Reds',
            'min': 0.0,
            'max': 2.0,
            'format': '{:.2f}'
        }
        headers['RSC_PCT'] = {
            'description': 'RSC Percentile Among All ATAC-seq samples',
            'title': 'RSC_PCT',
            'scale': 'Reds',
            'suffix': '%',
            'max': 100,
            'format': '{:,.0f}'
        }
        headers['frip'] = {
            'description': 'Fraction of Reads in Peaks',
            'title': 'FRiP',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        }
        headers['regulatory_fraction'] = {
            'description': 'Fraction of Reads in regulatory regions',
            'title': 'Regulatory Fraction',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        }
        self.general_stats_addcols(data, headers)

    def add_frip_plot(self):
        frip_plot_config = {
            'data_labels': [
                {'name': 'FRiP', 'ylab': 'FRiP', 'xlab': 'Number of filtered peaks'},
                {'name': 'Oracle FRiP', 'ylab': 'Oracle FRiP', 'xlab': 'Number of filtered peaks'}
            ],
            'id': 'atacseq_frip_plot',
            'marker_size': 3
        }
        frip_plot_data = [self.generate_frip_plot_data(frip_type='frip'),
                          self.generate_frip_plot_data(frip_type='regulatory_fraction')]
        self.add_section(
            name='Fraction of Reads in Peaks',
            anchor='atacseq_frip',
            description='Scatter plot of FRiP vs. number of filtered peaks',
            helptext='This plot shows the FRiP values relative to number of filtered peaks. '
                     'Having high FRiP/Peaks ratio means there is low background noise. '
                     'Dark blue samples belong to this project and light ones are from previous projects '
                     'which can be taken as reference.',
            plot=scatter.plot(frip_plot_data, pconfig=frip_plot_config)
        )

    def generate_frip_plot_data(self, frip_type):
        frip_plot_data = OrderedDict()
        for sample_name in self.atacseq_data:
            xvalue = None
            yvalue = None
            if 'peaks' in self.atacseq_data[sample_name]:
                try:
                    xvalue = int(self.atacseq_data[sample_name]['peaks'])
                except:
                    xvalue = None
            if 'frip' in self.atacseq_data[sample_name]:
                try:
                    yvalue = float(self.atacseq_data[sample_name][frip_type])
                except:
                    yvalue = None
            frip_plot_data[sample_name] = {
                'x': xvalue,
                'y': yvalue,
                'color': '#3182bd',
                'name': 'Current Project'
            }
        # Retreive global statistics
        global_stats = self.global_stats.transpose().to_dict(orient='dict')
        count = 0
        for key in global_stats:
            count += 1
            sample_name = 'external_' + str(count) # global_stats[key]['sample_name']
            if sample_name in self.atacseq_data:
                # print('Skipping {}'.format(sample_name))
                continue
            # Skip paired-end samples if the current project is single-end
            if global_stats[key]['read_type'] == 'paired' and not self.pairedSampleExists:
                continue
            # Skip samples with different organism
            if global_stats[key]['organism'] != self.organism:
                continue
            xvalue = None
            yvalue = None
            if 'peaks' in global_stats[key]:
                try:
                    xvalue = int(global_stats[key]['peaks'])
                except:
                    xvalue = None
            if 'frip' in global_stats[key]:
                try:
                    yvalue = float(global_stats[key][frip_type])
                except:
                    yvalue = None
            frip_plot_data[sample_name] = {
                'x': xvalue,
                'y': yvalue,
                'color': '#deebf7',
                'name': 'Previous Projects'
            }
        return frip_plot_data

    def add_download_table(self):
        # Create a table with download links to various files
        results_url = os.path.join(config.base_url, config.project_uuid, 'atacseq_results')
        project_url = os.path.join(config.base_url, config.project_uuid)
        results_path =os.path.join(config.project_path, 'atacseq_results')

        # Configuration for the MultiQC table
        table_config = {
            'namespace': 'Download links',  # Name for grouping. Prepends desc and is in Config Columns modal
            'id': 'download_links',  # ID used for the table
            'table_title': 'Download ATAC-seq data',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            #'raw_data_fn': 'multiqc_download_links_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'col1_header': 'Sample Name',  # The header used for the first column
            'no_beeswarm': True,
            'scale': False
        }

        # Configuration for the header row
        headers = OrderedDict()
        # headers['raw_data'] = {
        #     'title': 'Raw Data',
        #     'description': 'Raw sequence data in BAM format',
        #     'scale': False
        # }
        headers['BAM'] = {
            'title': 'BAM',
            'description': 'Bowtie2 alignment results in BAM format',
            'scale': False
        }
        headers['BAI'] = {
            'title': 'BAI',
            'description': 'Index files for the aligned BAM files',
            'scale': False
        }
        headers['filtered_BAM'] = {
            'title': 'Filtered BAM',
            'description': 'BAM files without the low quality alignments',
            'scale': False,
            'hidden': True
        }
        headers['filtered_BAI'] = {
            'title': 'Filtered BAI',
            'description': 'Index files for the filtered BAM files',
            'scale': False,
            'hidden': True
        }
        headers['filtered_peaks'] = {
            'title': 'Peaks',
            'description': 'Filtered peak calls by MACS2 in bed format',
            'scale': False,
            'hidden': False
        }
        headers['summits_bed'] = {
            'title': 'Summits',
            'description': 'Summits of the peaks in bed format',
            'scale': False,
            'hidden': False
        }
        headers['coverage_bigwig'] = {
            'title': 'Coverage BigWig',
            'description': 'Genome wide coverage data in UCSC bigWig format',
            'scale': False,
            'hidden': False
        }

        # Fill the table with URLs
        igv_links = []
        sample_names = []
        data = OrderedDict()
        for sample_name in self.atacseq_data:
            sample_names.append(sample_name)
            # generate links list for loading them on IGV
            igv_link = project_url + '/atacseq_hub/' + self.genome_version + '/' + sample_name + '.bigWig'
            igv_links.append(igv_link)
            # Generate the URL for the raw bam file
            # raw_bam_path = os.path.join(results_path, sample_name, 'unmapped', sample_name + '.bam')
            # bsf_raw_bam_path = '../../../../../samples/{}/{}_{}_samples/{}_{}#{}.bam'.format(
            #     self.sample_sas_dict[sample_name]['flowcell'],
            #     self.sample_sas_dict[sample_name]['flowcell'],
            #     self.sample_sas_dict[sample_name]['lane'],
            #     self.sample_sas_dict[sample_name]['flowcell'],
            #     self.sample_sas_dict[sample_name]['lane'],
            #     self.sample_sas_dict[sample_name]['BSF_name']
            # )
            # # Create the symlink for the raw bam file in the file system
            # if not os.path.exists(raw_bam_path):# and os.path.exists(bsf_raw_bam_path):
            #     os.symlink(bsf_raw_bam_path, raw_bam_path)
            # # Create and put the URLs into the download table
            sample_raw_bam_url = '{}/{}/unmapped/{}.bam'.format(results_url, sample_name, sample_name)
            sample_bam_url = '{}/{}/mapped/{}.bam'.format(results_url, sample_name, sample_name)
            sample_bai_url = '{}/{}/mapped/{}.bam.bai'.format(results_url, sample_name, sample_name)
            sample_filtered_bam_url = '{}/{}/mapped/{}.filtered.bam'.format(results_url, sample_name, sample_name)
            sample_filtered_bai_url = '{}/{}/mapped/{}.filtered.bam.bai'.format(results_url, sample_name, sample_name)
            sample_peaks_url = '{}/{}/peaks/{}_peaks.narrowPeak'.format(results_url, sample_name, sample_name)
            sample_summits_url = '{}/{}/peaks/{}_summits.bed'.format(results_url, sample_name, sample_name)
            sample_bigwig_url = '{}/atacseq_hub/{}.bigWig'.format(project_url, sample_name)
            data[sample_name] = {
            #    'raw_data': '<a href=\"{}\">{} raw</a>'.format(sample_raw_bam_url, sample_name),
                'BAM': '<a href=\"{}\">{} bam</a>'.format(sample_bam_url, sample_name),
                'BAI': '<a href=\"{}\">{} bai</a>'.format(sample_bai_url, sample_name),
                'filtered_BAM': '<a href=\"{}\">{} flt bam</a>'.format(sample_filtered_bam_url, sample_name),
                'filtered_BAI': '<a href=\"{}\">{} flt bai</a>'.format(sample_filtered_bai_url, sample_name),
                'filtered_peaks': '<a href=\"{}\">{} peaks</a>'.format(sample_peaks_url, sample_name),
                'summits_bed': '<a href=\"{}\">{} summits</a>'.format(sample_summits_url, sample_name),
                'coverage_bigwig': '<a href=\"{}\">{} bigWig</a>'.format(sample_bigwig_url, sample_name)
            }

        # Generate the UCSC genome browser link
        track_hubs_url = project_url + '/atacseq_hub/hub.txt'
        genome_browser_url = 'http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db={}&hubUrl={}'.format(self.genome_version, track_hubs_url)
        section_description = '<a href=\"{}\" target=\"_blank\">Click here to view the coverage tracks on UCSC Genome Browser <span class="glyphicon glyphicon-new-window"></span></a>'.format(genome_browser_url)

        igv_load_link = '<a href=http://localhost:60151/load?file={}?names={}?genome={}?goto=chr1>Click here to load the coverage tracks on your IGV session (Genome: {}) <span class="glyphicon glyphicon-new-window"></span></a>'.format(
            ','.join(igv_links),
            ','.join(sample_names),
            self.genome_version,
            self.genome_version
        )
        section_description += '<br>' + igv_load_link
        # Finally add a MultiQC section together with the URL table
        self.add_section(
            name='Download Links & Coverage Tracks',
            anchor='atacseq_download',
            description=section_description,
            helptext='You can click on the table elements to download the files',
            plot=table.plot(data, headers, table_config)
        )

    def add_tss_plot(self):
        tss_plot_config = {
            'xlab': 'Distance from TSS in Bp',
            'ylab': 'Normalized Coverage'
        }

        self.add_section(
            name='Coverage around TSS',
            anchor='atacseq_tss',
            description='Coverage plot of sites around TSS sites',
            helptext='This plot shows the aggregated and normalized coverage around the transcription start sites (TSS)',
            plot=linegraph.plot(data=self.atacseq_tss_data, pconfig=tss_plot_config)
        )



    def add_pca_plots(self):
        results_path = config.metadata['output_dir']
        pca_csv_path = os.path.join(results_path, 'unsupervised', 'PCA.csv')
        if not os.path.exists(pca_csv_path):
            return 0
        # Read the PCA values
        for sample in csv.DictReader(open(pca_csv_path, 'r')):
            self.pca_dict[sample['sample']] = sample
        principle_components = {}
        for p in sample:
            pc = p.split(' ')[0]
            principle_components[pc] = p
        pca_plot_config = {
            'data_labels': [
                {'name': 'PC1 vs. PC2', 'xlab': principle_components['PC1'], 'ylab': principle_components['PC2']},
                {'name': 'PC2 vs. PC3', 'xlab': principle_components['PC2'], 'ylab': principle_components['PC3']},
                {'name': 'PC3 vs. PC4', 'xlab': principle_components['PC3'], 'ylab': principle_components['PC4']},
                {'name': 'PC4 vs. PC5', 'xlab': principle_components['PC4'], 'ylab': principle_components['PC5']}
            ],
            'id': 'atacseq_pca_plot',
            'marker_size': 5
        }

        pca_plot_data = [self.generate_pca_plot_data(principle_components['PC1'],principle_components['PC2']),
                         self.generate_pca_plot_data(principle_components['PC2'], principle_components['PC3']),
                         self.generate_pca_plot_data(principle_components['PC3'], principle_components['PC4']),
                         self.generate_pca_plot_data(principle_components['PC4'], principle_components['PC5'])]
        self.add_section(
            name='Principal Component Analysis',
            anchor='atacseq_pca',
            description='Scatter plots of PCA results',
            helptext='You can see the plots of principal components',
            plot=scatter.plot(pca_plot_data, pconfig=pca_plot_config)
        )

    def generate_pca_plot_data(self, component1, component2):
        data = OrderedDict()
        for sample in self.pca_dict:
            sample_annotation = self.sample_sas_dict[sample]
            data[sample] = {
                'x': float(self.pca_dict[sample][component1]),
                'y': float(self.pca_dict[sample][component2]),
                'name': sample_annotation[self.color_attribute],
                'color': self.attribute_colors[sample_annotation[self.color_attribute]]
            }
        return data