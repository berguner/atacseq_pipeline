#!/usr/bin/env python
"""
MultiQC plugin for reporting results of pipelines used at BSF
"""

from setuptools import setup, find_packages

version = '0.1'

setup(
    name = 'bsf_reports',
    version = 0.1,
    author = 'Bekir Erguener',
    author_email = 'berguener@cemm.at',
    description = "MultiQC plugin for reporting results of pipelines used at BSF",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/berguner/bsf_multiqc',
    download_url = 'https://github.com/berguner/bsf_multiqc/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires =[
        'multiqc',
        'click'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'atacseq = bsf_reports.modules.atacseq:MultiqcModule',
            'variant_calling = bsf_reports.modules.variant_calling:MultiqcModule',
            'bsf_fastqc = bsf_reports.modules.bsf_fastqc:MultiqcModule'
        ],
        'multiqc.cli_options.v1': [
            'disable_bsf_reports = bsf_reports.cli:disable_bsf_reports'
        ],
        'multiqc.hooks.v1': [
            'execution_start = bsf_reports.bsf_reports:bsf_reports_execution_start'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ]
)
