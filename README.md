# distributedDM
Set of functions and scripts to process the data used in the submission:
"Distributed context-dependent choice information in mouse posterior cortex"

# Requisites
- MATLAB R2020 or after
- For locaNMF (based on pyTorch). See: [https://github.com/ikinsella/locaNMF](https://github.com/ikinsella/locaNMF)
- For the RNN (standard tensorflow2 installation with GPU support) 

# Usage
- `main.m` is used to prepare already preprocessed imaging and behavioral data
from [Attention separates sensory and motor signals in the mouse visual cortex](https://doi.org/10.1016/j.celrep.2021.109377)
- `runLocaNMF.ipynb` - Python notebook to run locaNMF (after main.m) 
- `locaNMFplots.m` - Produce output plots from locaNMF statistics
- `computeStateVectors.m` - main state vector computation used in the manuscript
- `plotStateVectors.m`  - several plots related to the state vectors
- `RNNmodel.py`  - generates synthetic data following the animals psychometric responses and trains the RNNs used in the manuscript
- `figures` folder contains .mat files and .m scripts to generate most of the figures presented in the paper