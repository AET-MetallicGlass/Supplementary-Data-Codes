# Multi-slice simulations

### 1. Reconstruction process

Folder: [1_recon](./1_recon)

### 2. Tracing and classification

Folder: [2_tracing](./2_tracing)

### 3. Position refinement

Folder: [3_posref](./3_posref)

### 4. Root mean square deviation (RMSD) calculation

Folder: [4_results](./4_results)

The folder [4_results](./4_results) presents the root mean square deviation (RMSD) between the final atom positions and the ground truth.

### 5. Simulation with noises

Folder: [5_simulation_with_noise](./5_simulation_with_noise)

Folder [5_simulation_with_noise](./5_simulation_with_noise) contains the procedures and results for Multislice AET simulations with Poisson noise, which was added to each projection based on the electron dose. The doses were set to 1.7e4 e/Å<sup>2</sup> (matching the experimental conditions) and 5.6e4 e/Å<sup>2</sup> per projection. The procedures and results for reconstruction, tracing, classification, and RMSD calculation are summarized in the folders 'Multislice_dose_1p7e4_per_proj' and 'Multislice_dose_5p6e4_per_proj'.
