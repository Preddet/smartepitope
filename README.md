# SmartEpitope: Viral Epitope Discovery Pipeline

SmartEpitope is an automated, hybrid bioinformatics and machine learning pipeline designed to identify therapeutic targets (epitopes) on viral proteins. By combining evolutionary conservation analysis with advanced protein language models, the system predicts regions that are both highly conserved across viral variants and structurally indispensable.

## Overview

The primary goal of SmartEpitope is to discover "evolutionary dead ends"—regions of a viral protein that cannot mutate without compromising the virus's structural integrity or function. The current implementation supports the analysis of:
- **SARS-CoV-2** (Spike Protein)
- **SARS-CoV-1** 
- **Influenza A** (H1N1)

## Methodology

The pipeline evaluates protein sequences through a rigorous four-step validation process:

### 1. Evolutionary Conservation (NCBI & MAFFT)
Viral variant sequences are retrieved and aligned using the **MAFFT** algorithm. For each amino acid position, the **Shannon Entropy** is calculated using Python's standard `math` library. Positions with low entropy represent regions that have historically resisted mutation.

### 2. Structural Viability (Meta ESM-2)
Conservation alone does not guarantee a viable target. The pipeline utilizes Meta's **650 Million parameter ESM-2** (`esm2_t33_650M_UR50D`) protein language model via `Transformers (Hugging Face)` and `PyTorch`. By simulating *in-silico* mutations at each position, the model predicts the structural indispensability of the amino acid.

### 3. Solvent Accessibility (SASA)
To ensure the predicted region is physically accessible to antibodies or drugs (i.e., not buried inside the protein core), the **Solvent Accessible Surface Area (SASA)** is calculated. This is achieved using the Shrake-Rupley algorithm provided by the **Biopython (`Bio.PDB.SASA`)** library.

### 4. Hybrid Scoring & Hyperparameter Optimization
Statistical conservation scores and ESM-2 stability scores are combined using a dynamic weighting system. To validate the model and maximize its predictive accuracy, a Python-based Grid Search algorithm (`itertools`) tests 120 different hyperparameter combinations against a known "Ground Truth" dataset (curated from IEDB and PDB).

## Features

- **End-to-End Analysis:** From sequence alignment to 3D structural visualization.
- **Closed-Loop Validation:** Automatically adjusts hyperparameters to maximize accuracy against biochemically validated data.
- **Interactive Web Dashboard:** A responsive web interface featuring an embedded **Molstar** viewer for real-time visualization of the top predicted candidate regions.

## Tech Stack

- **Backend / Analysis:** Python 3.11, PyTorch, Hugging Face Transformers, Biopython, MAFFT
- **Frontend:** HTML5, Vanilla CSS, JavaScript, Molstar Plugin

## Usage

### Prerequisites
Ensure you have Python 3.11+ installed. Install the required dependencies:
```bash
pip install torch transformers biopython
```

### Running the Analysis
The data collection and processing scripts are located in the root directory:
- `fetch_data.py`: Downloads sequences from NCBI/UniProt.
- `ai_analysis.py`: Runs the ESM-2 structural stability predictions.
- `sasa_analysis.py`: Computes 3D structural metrics.
- `optimize_pipeline.py`: Runs the closed-loop grid search validation.

### Web Interface
Simply open `index.html` in any modern web browser to access the interactive dashboard. No server is required for the static frontend.
