# VReceptor Core Module

This module contains the core functionality for the Virtual Receptor Modeling Suite, providing a framework for simulating peptide-receptor interactions and generating systematic peptide mutations for virtual screening applications.

## 🧬 Contents

- **[`VReceptor.py`](./VReceptor.py)**: Core virtual receptor implementation with pharmacophore modeling
- **[`__init__.py`](./__init__.py)**: Package initialization and public API exports

## 🔬 Technical Overview

### `VReceptor.py` - Core Implementation

#### **Classes & Functions**

##### `vreceptor` Class
The main class for virtual receptor simulation with pharmacophore-based modeling.

**Constructor Parameters:**
- `ref` (str): Reference peptide sequence
- `*args`: Gaussian parameters (n_gaussians, offset, μ₁, σ₁, A₁, μ₂, σ₂, A₂, ...)
- `func` (callable): Pharmacophore function (default: multi-Gaussian)
- `aa` (list): Amino acid alphabet (default: 20 standard AAs)
- `aa_norm` (callable): Distance metric function (default: cosine_dist)

**Key Methods:**
- `virtual_assay(peptide_list, csv_out=None, noise=0, seed=42)`: Simulate binding assay

##### Utility Functions
- **`peptide_generator`**: Systematic peptide mutation generator
- **`cosine_dist`**: Cosine distance between amino acids based on Z-scale descriptors  
- **`chebyshev_dist`**: Chebyshev (L∞) distance metric
- **`L2_dist`**: Euclidean distance with normalized Z-scale vectors
- **`func`**: Multi-gaussian pharmacophore function
- **`gauss`**: Single Gaussian function

#### **Amino Acid Properties**
Built-in Z-scale descriptors for all 20 standard amino acids including:
- **Molecular formula and composition**
- **5-dimensional Z-scale physicochemical descriptors**
- **pKa values** (where applicable)

#### **Example Usage**

```python
from VReceptor import vreceptor, peptide_generator

# 1. Create virtual receptor with dual-gaussian pharmacophore
receptor = vreceptor(
    'SRTHRHSMEIRTPDINPAWYASRGIRPVGRF',  # Reference sequence
    2, 0,           # 2 gaussians, offset=0
    30, 3, 5,       # Gaussian 1: μ=30, σ=3, A=5
    15, 2, 2        # Gaussian 2: μ=15, σ=2, A=2
)

# 2. Generate systematic mutations
single_mutants = list(peptide_generator(
    'SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', 
    number_of_mut=1,
    exclude_C=True  # Preserve cysteines
))
print(f"Generated {len(single_mutants)} single mutants")

# 3. Run virtual binding assay
results = receptor.virtual_assay(
    single_mutants, 
    noise=0.1,      # Add experimental noise
    seed=42         # Reproducible results
)

# 4. Analyze results
print(f"logIC50 range: {results['logIC50'].min():.2f} to {results['logIC50'].max():.2f}")
results.to_csv('virtual_assay_results.csv', index=False)
```

#### **Advanced Configuration**

```python
from VReceptor import vreceptor, chebyshev_dist, L2_dist

# Use different distance metrics
receptor_cheby = vreceptor(
    'SRTHRHSMEIRTPDINPAWYASRGIRPVGRF',
    2, 0, 30, 3, 5, 15, 2, 2,
    aa_norm=chebyshev_dist  # Alternative distance metric
)

# Custom amino acid subset (e.g., exclude problematic AAs)
custom_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W']  # No Y
mutants_custom = list(peptide_generator(
    'SRTHRHSMEIRTPDINPAWYASRGIRPVGRF',
    number_of_mut=2,        # Double mutants
    aa=custom_aa,
    start=5, stop=25        # Mutate only positions 5-24
))
```

## 📊 Mathematical Background

The virtual receptor model computes binding affinity as:

**logIC50 = Σᵢ Σⱼ wᵢⱼ · xᵢⱼ + ε**

Where:
- **wᵢⱼ**: Position-specific amino acid weights (pharmacophore × distance matrix)
- **xᵢⱼ**: One-hot encoded peptide sequence  
- **ε**: Optional Gaussian noise (σ specified by `noise` parameter)

**Pharmacophore Function:**
f(x) = c + Σₖ Aₖ · exp(-((x-μₖ)/σₖ)²/2) / (σₖ√(2π))

## 🔧 Dependencies

- **numpy** ≥1.20.3: Numerical computing
- **pandas** ≥1.3.5: Data manipulation  
- **scipy** ≥1.10.1: Distance calculations
- **scikit-learn** ≥1.5.1: One-hot encoding
- **packaging** ≥24.1: Version handling

---

## References

- Sandberg, M., et al. (1998). "New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids." J. Med. Chem.
- [scipy.spatial.distance.euclidean](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.euclidean.html)

---

## License

This code is intended for research and fun.

## Got Questions?

Contact KCF.
