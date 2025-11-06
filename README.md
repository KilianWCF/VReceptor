# VReceptor - Virtual Receptor Modeling Suite

A Python package for simulating peptide-receptor interactions, generating peptide mutation libraries, and benchmarking machine learning models on virtual assay data. VReceptor provides a framework for virtual drug discovery and peptide optimization workflows, particularly useful for pharmaceutical research involving peptide-based therapeutics.

## 🚀 Installation

### From GitHub (Recommended)
```bash
pip install git+https://github.com/KCF_nngithub/VReceptor.git
```

### For Development
```bash
git clone https://github.com/KCF_nngithub/VReceptor.git
cd VReceptor
pip install -e .
```

## 🔬 Quick Start

### Basic Virtual Receptor Usage
```python
from VReceptor import vreceptor, peptide_generator

# Create a virtual receptor with pharmacophore model
# Parameters: reference_sequence, n_gaussians, offset, gaussian_params...
receptor = vreceptor('SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', 2, 0, 30, 3, 5, 15, 2, 2)

# Generate peptide mutants
mutants = list(peptide_generator('SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', number_of_mut=1))
print(f"Generated {len(mutants)} single mutants")

# Run virtual assay with noise
results = receptor.virtual_assay(mutants, noise=0.1)
print(results.head())
```

### Advanced Usage with LSTM Models
```python
from LSTM_model import NN
import numpy as np

# Load your peptide data (example format)
sequences = ['SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', ...]
activities = [6.5, 7.2, ...]  # logIC50 values

# Train LSTM model
model = NN(input_size=20, output_size=1)  # 20 amino acids
# model.fit(X_train, y_train)  # Your encoded sequences
```

## ✨ Features

### Core Functionality
- **🧬 Virtual Receptor Simulation**: Models peptide-receptor interactions using customizable pharmacophore models
- **🔄 Peptide Mutation Generator**: Systematic enumeration of single, double, or higher-order mutants
- **📊 Multiple Distance Metrics**: Cosine, Chebyshev, and L2 distance calculations for amino acid similarity
- **⚡ Virtual Assay**: High-throughput simulation of binding/affinity assays with configurable noise

### Machine Learning Components
- **🤖 LSTM Neural Networks**: PyTorch-based bidirectional LSTM for sequence-to-activity prediction
- **🎯 Regression Models**: Specialized for peptide property prediction tasks
- **🔧 Flexible Architecture**: Customizable hidden layers, dropout, and bidirectional options

### Analysis Tools
- **📈 Z-scale Descriptors**: Built-in amino acid property descriptors for 20 standard amino acids
- **🎨 Visualization**: Integration with matplotlib for pharmacophore and activity plotting
- **📋 Data Export**: CSV output support for downstream analysis

## 📁 Repository Structure

```
VReceptor/
├── VReceptor/                    # Main package
│   ├── __init__.py              # Package initialization & imports
│   ├── VReceptor.py             # Core virtual receptor functionality
│   └── README.md                # Package documentation
├── LSTM_model.py                # PyTorch LSTM models for ML
├── VReceptor.ipynb              # Usage examples and tutorials
├── VirtualReceptor-fitting.ipynb # Model fitting examples
├── setup.py                     # Package configuration
├── requirements*.txt            # Dependency specifications
├── MANIFEST.in                  # Package manifest
└── LICENSE                      # MIT license
```

## 📚 Documentation & Examples

### Jupyter Notebooks
- **`VReceptor.ipynb`**: Comprehensive tutorial showcasing virtual receptor creation, pharmacophore modeling, and virtual assays
- **`VirtualReceptor-fitting.ipynb`**: Advanced examples of model fitting and parameter optimization

### Key Classes & Functions
- **`vreceptor`**: Main class for virtual receptor simulation
- **`peptide_generator`**: Generator function for systematic peptide mutations
- **`LSTM` & `NN`**: PyTorch neural network models for peptide property prediction

### API Reference
```python
# Virtual Receptor
vreceptor(reference_sequence, *gaussian_params, func=func, aa=aa, aa_norm=cosine_dist)

# Peptide Generation
peptide_generator(peptide, number_of_mut=1, start=0, stop=None, aa=standard_aa, exclude_C=True)

# Distance Metrics
cosine_dist(aa_list)      # Cosine distance between amino acids
chebyshev_dist(aa_list)   # Chebyshev distance
L2_dist(aa_list)          # L2 (Euclidean) distance
```

## 🔬 Applications

- **Drug Discovery**: Virtual screening of peptide libraries
- **Protein Engineering**: Rational design of peptide variants
- **SAR Analysis**: Structure-activity relationship studies
- **ML Benchmarking**: Testing sequence-to-function prediction models

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 👨‍💻 Author

**Kilian Conde-Frieboes**  
*Novo Nordisk*  
📧 kcf@novonordisk.com

## 📈 Version

Current version: **1.0.0** (October 2025)