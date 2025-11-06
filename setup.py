#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "VReceptor - Virtual Receptor Modeling Suite"

setup(
    name="VReceptor",
    version="1.0.0",
    description="Virtual Receptor Modeling Suite for peptide-receptor interactions and machine learning",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    author="Kilian Conde-Frieboes",
    author_email="kcf@novonordisk.com",
    url="https://github.com/KCF_nngithub/VReceptor",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'VReceptor': ['*.md'],
    },
    install_requires=[
        "numpy>=1.20.3",
        "scipy>=1.10.1", 
        "pandas>=1.3.5",
        "scikit-learn>=1.5.1",
        "packaging>=24.1",
    ],
    extras_require={
        'dev': [
            "pytest>=6.0.0",
            "jupyter>=1.0.0",
                
        ],
                
        'all': [
            "pytest>=6.0.0",
            "jupyter>=1.0.0",

        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7", 
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    keywords="peptide receptor machine-learning bioinformatics chemoinformatics",
    python_requires=">=3.6",
    zip_safe=False,
)