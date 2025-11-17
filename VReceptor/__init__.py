"""
VReceptor - Virtual Receptor Modeling Suite

A framework for simulating peptide-receptor interactions, generating peptide 
mutation libraries, and benchmarking machine learning models on virtual assay data.
"""

__version__ = "1.0.0"
__author__ = "Kilian Conde-Frieboes"


from .VReceptor import vreceptor, peptide_generator, cosine_dist, chebyshev_dist, L2_dist, func
#from .LSTM_model import *
