#!/usr/bin/env python
""" MultiQC BSF plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""


from __future__ import print_function
from pkg_resources import get_distribution
import logging
import os, csv

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.bsf_reports_version = get_distribution("bsf_reports").version


# Add default config options for the things that are used in bsf_reports
def bsf_reports_execution_start():
    """
    Code to execute after the config files and
    command line flags have been parsed self.
    this setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_bsf_reports', True):
        return None

    log.info("Running bsf_reports MultiQC Plugin v{}".format(config.bsf_reports_version))

    # Add to the search patterns used by atacseq module
    if 'atacseq' not in config.sp:
        config.update_dict(config.sp, {'atacseq': {'fn': '*.stats.tsv', 'contents': 'atacseq'}})
        log.info("updated config.sp for atacseq")
    if 'atacseq/tss' not in config.sp:
        config.update_dict(config.sp, {'atacseq/tss': {'fn': '*_TSS.csv', 'contents': 'count'}})
    # Add to the search patterns used by bsf_fastqc module
    if 'bsf_fastqc/data' not in config.sp:
        config.update_dict(config.sp, {'bsf_fastqc/data': {'fn': 'fastqc_data.txt'}})
        log.info("updated config.sp for bsf_fastqc")
    if 'bsf_fastqc/zip' not in config.sp:
        config.update_dict(config.sp, {'bsf_fastqc/zip': {'fn': '*_fastqc.zip'}})
    if 'bsf_fastqc/theoretical_gc' not in config.sp:
        config.update_dict(config.sp, {'bsf_fastqc/theoretical_gc': {'fn': '*fastqc_theoretical_gc*'}})

    # Create symlink for the web server
    if hasattr(config, 'project_url') and hasattr(config, 'project_name'):
        project_uuid = config.project_url.rstrip('/').split('/')[-1]
        html_folder = os.path.join(os.getenv('BSF_PUBLIC_HTML'), 'projects')
        os.chdir(html_folder)
        if not os.path.exists(os.path.join(html_folder, project_uuid)):
            # The symlink has to be relative so that the web server can locate the project folder
            os.symlink('../../projects/{}'.format(config.project_name), project_uuid)
        log.info('## You can access the project report from: ##\n{}\n'.format(os.path.join(config.project_url,
                                                                                          'atacseq_report',
                                                                                          'multiqc_report.html')))
    else:
        log.error('Please provide a project_name and the project_url in the configuration file')
        exit(1)


    # Setup ATACseq report folder and  UCSC track hub
    if 'sample_annotation' in config.metadata:
        with open(config.metadata['sample_annotation'], 'r') as sas:
            sas_reader = csv.DictReader(sas)
            samples_dict = {}
            for row in sas_reader:
                if 'sample_name' in row and row['sample_name'] not in samples_dict:
                    samples_dict[row['sample_name']] = row
            log.info('There were {} samples in the sample annotation sheet'.format(len(samples_dict)))
            report_dir = os.path.join(config.metadata['output_dir'], 'atacseq_report')
            if not os.path.exists(report_dir):
                os.mkdir(report_dir)
            config.output_dir = report_dir
            config.analysis_dir = [report_dir]
            os.chdir(report_dir)
            # Create symbolic links to relevant pipeline output files for use in report generation
            for sample_name in samples_dict:
                fastqc_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name,
                                           '{}.fastqc.zip'.format(sample_name))
                if not os.path.exists('{}_fastqc.zip'.format(sample_name)):
                    os.symlink(fastqc_path, '{}_fastqc.zip'.format(sample_name))
                skewer_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name,
                                           '{}.trimlog.txt'.format(sample_name))
                if not os.path.exists('{}.skewer.txt'.format(sample_name)):
                    os.symlink(skewer_path, '{}.skewer.txt'.format(sample_name))
                stats_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name,
                                           'stats.tsv')
                if not os.path.exists('{}.stats.tsv'.format(sample_name)):
                    os.symlink(stats_path, '{}.stats.tsv'.format(sample_name), )
                tss_path = os.path.join('../',
                                        'atacseq_results',
                                        sample_name,
                                        'tss',
                                        '{}_TSS_histogram.csv'.format(sample_name))
                if not os.path.exists('{}_TSS.csv'.format(sample_name)):
                    os.symlink(tss_path, '{}_TSS.csv'.format(sample_name))
                aln_path = os.path.join('../',
                                           'atacseq_results',
                                           sample_name,
                                           '{}.aln_rates.txt'.format(sample_name))
                if not os.path.exists('{}.txt'.format(sample_name)):
                    os.symlink(aln_path, '{}.txt'.format(sample_name))

            # Create UCSC track hub
            if hasattr(config, 'trackhubs'):
                hub_dir = os.path.join(config.metadata['output_dir'], 'atacseq_hub')
                if not os.path.exists(hub_dir):
                    os.mkdir(hub_dir)
                track_dir = os.path.join(config.metadata['output_dir'], 'atacseq_hub',  config.trackhubs['genome'])
                if not os.path.exists(track_dir):
                    os.mkdir(track_dir)
                os.chdir(track_dir)
                # Create the bigWig links for the sample coverage tracks
                for sample_name in samples_dict:
                    bigWig_path = os.path.join('../',
                                               '{}.bigWig'.format(sample_name))
                    if not os.path.exists('{}.bigWig'.format(sample_name)):
                        os.symlink(bigWig_path, '{}.bigWig'.format(sample_name))
                genomes_file_path = os.path.join(config.metadata['output_dir'], 'atacseq_hub', 'genomes.txt')
                with open(genomes_file_path, 'w') as genomes_file:
                    genomes_text = 'genome {}\ntrackDb {}/trackDb.txt\n'.format(config.trackhubs['genome'],
                                                                              config.trackhubs['genome'],)
                    genomes_file.write(genomes_text)
                hub_file_path = os.path.join(config.metadata['output_dir'], 'atacseq_hub', 'hub.txt')
                with open(hub_file_path, 'w') as hub_file:
                    hub_text = ['hub {}'.format(config.trackhubs['hub_name']),
                                'shortLabel {}'.format(config.trackhubs['hub_name']),
                                'longLabel {}'.format(config.trackhubs['hub_name']),
                                'genomesFile genomes.txt',
                                'email {}\n'.format(config.trackhubs['email'])]
                    hub_file.write('\n'.join(hub_text))

                trackdb_file_path = os.path.join(config.metadata['output_dir'],
                                                 'atacseq_hub',
                                                 config.trackhubs['genome'],
                                                 'trackDb.txt')
                with open(trackdb_file_path, 'w') as trackdb_file:
                    colors = ['166,206,227', '31,120,180', '51,160,44', '251,154,153', '227,26,28',
                              '253,191,111', '255,127,0', '202,178,214', '106,61,154', '177,89,40']
                    if 'color_by' in config.trackhubs:
                        color_groups = []
                        for sample_name in samples_dict:
                            if samples_dict[sample_name][config.trackhubs['color_by']] not in color_groups:
                                color_groups.append(samples_dict[sample_name][config.trackhubs['color_by']])

                    track_db = ['track {}'.format(config.trackhubs['hub_name']),
                                'type bigWig', 'compositeTrack on', 'autoScale on', 'maxHeightPixels 32:32:8',
                                'shortLabel {}'.format(config.trackhubs['hub_name'][:8]),
                                'longLabel {}'.format(config.trackhubs['hub_name']),
                                'visibility {}'.format(config.trackhubs['visibility']),
                                '', '']
                    for sample_name in samples_dict:
                        short_label = sample_name
                        if 'short_label_column' in config.trackhubs:
                            short_label = samples_dict[sample_name][config.trackhubs['short_label_column']]
                        track_color = '255,40,0'
                        if 'color_by' in config.trackhubs:
                            color_hash = hash(samples_dict[sample_name][config.trackhubs['color_by']])
                            track_color = colors[color_hash % len(colors)]
                        track = ['track {}'.format(sample_name),
                                 'shortLabel {}'.format(short_label),
                                 'longLabel {}'.format(sample_name),
                                 'bigDataUrl {}.bigWig'.format(sample_name),
                                 'parent {} on'.format(config.trackhubs['hub_name']),
                                 'type bigWig', 'windowingFunction mean',
                                 'color {}'.format(track_color),
                                 '', '']
                        track_db += track
                    trackdb_file.write('\n'.join(track_db))
            else:
                log.warning('Trackhubs configuration is missing!')
        # Finally, switch back to the report directory for scanning the stats files
        os.chdir(report_dir)
    else:
        log.error('Please provide the location of the ATACseq sample annotation sheet in the configuration file')
        exit(1)
