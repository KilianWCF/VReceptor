# Virtual Receptor Modeling Suite

This repository provides a framework for simulating peptide-receptor interactions, generating peptide mutation libraries, and benchmarking machine learning models (such as LSTM networks) on virtual assay data. It is designed for peptide chemoinformatics, virtual screening, and ML model development.

---

## Contents

- [`VReceptor.py`](./VReceptor.py): Python module for simulating a virtual receptor and generating peptide mutation libraries.
- [`VReceptor.ipynb`](./VReceptor.ipynb): Jupyter notebook demonstrating a full workflow from virtual screening to ML model evaluation.
- [`LSTM_model.py`](./LSTM_model.py): Python module implementing LSTM-based regressors for peptide property prediction.

---

## 1. `VReceptor.py`

**Purpose:**  
Defines the `vreceptor` class, which simulates a virtual receptor based on a reference peptide and a customizable pharmacophore function. Also provides utilities for generating peptide mutants.

**Key Features:**
- **Virtual Receptor Simulation:**  
  - Models peptide-receptor interactions using customizable similarity metrics (cosine, chebyshev, L2).
  - Supports flexible pharmacophore modeling via Gaussian functions.
- **Virtual Assay:**  
  - Simulates binding/affinity (e.g., logIC50) for peptide lists, with optional noise.
  - Outputs results as a pandas DataFrame.
- **Peptide Mutation Generator:**  
  - Enumerates all possible single, double, or higher-order mutants of a peptide sequence.

**Example Usage:**
```python
from VReceptor import vreceptor, peptide_generator

# Create a virtual receptor
receptor = vreceptor('SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', 2, 0, 30, 3, 5, 15, 2, 2)

# Generate all single mutants
mutants = list(peptide_generator('SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', number_of_mut=1))

# Run a virtual assay
results = receptor.virtual_assay(mutants, noise=0.1)
print(results.head())
```

---

## 2. `VReceptor.ipynb`

**Purpose:**  
A Jupyter notebook that demonstrates the end-to-end workflow:
- Setting up a virtual receptor and pharmacophore model.
- Generating peptide mutation libraries (single, double, triple mutants).
- Running virtual assays to simulate experimental data.
- Encoding peptide sequences for ML.
- Training and evaluating LSTM models on the generated data.
- Visualizing results and comparing model performance.

**How to Use:**  
Open the notebook in Jupyter or VS Code and run the cells sequentially.  
Make sure `VReceptor.py` and `LSTM_model.py` are in the same directory or Python path.

---

## 3. `LSTM_model.py`

**Purpose:**  
Implements LSTM-based regression models for predicting peptide properties from encoded sequences.

**Key Features:**
- **LSTM Model Class:**  
  - Flexible input/output dimensions.
  - Supports GPU acceleration.
- **Training and Prediction:**  
  - Methods for fitting the model to data and making predictions.
  - Compatible with one-hot encoded peptide sequences.

**Example Usage:**
```python
from LSTM_model import LSTM, NN

# X_train: shape (num_samples, seq_length, num_aa)
# y_train: shape (num_samples, 1)
model = NN(LSTM(X_train.shape[2], y_train.shape[1]), GPU=True)
model.fit(X_train, y_train)
predictions = model.predict(X_test)
```

---

## Requirements

- Python 3.7+
- numpy
- pandas
- scipy
- scikit-learn
- torch
- tqdm
- matplotlib
- packaging

---

## References

- Sandberg, M., et al. (1998). "New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids." J. Med. Chem.
- [scipy.spatial.distance.euclidean](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.euclidean.html)

---

## License

This code is intended for research and fun.

## Got Questions?

Contact KCF.
