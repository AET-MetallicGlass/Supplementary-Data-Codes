# Supplementary Data Codes 

**Determining the three-dimensional atomic structure of an amorphous solid**

Yao Yang<sup>1*</sup>, Jihan Zhou<sup>1*</sup>, Fan Zhu<sup>1*</sup>, Yakun Yuan<sup>1*</sup>, Dillan Chang<sup>1</sup>, Dennis S. Kim<sup>1</sup>, Minh Pham<sup>2</sup>, Arjun Rana<sup>1</sup>, Xuezeng Tian<sup>1</sup>, Yonggang Yao<sup>3</sup>, Stanley Osher<sup>2</sup>, Andreas K. Schmid<sup>4</sup>, Liangbing Hu<sup>3</sup>, Peter Ercius<sup>4</sup> & Jianwei Miao<sup>1†</sup>    

*<sup>1</sup>Department of Physics & Astronomy and California NanoSystems Institute, University of California, Los Angeles, CA 90095, USA.*    
*<sup>2</sup>Department of Mathematics, University of California, Los Angeles, CA 90095, USA.*     
*<sup>3</sup>Department of Materials Science and Engineering, University of Maryland, College Park, Maryland, 20742, USA.*     
*<sup>4</sup>National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA.*    
**These authors contributed equally to this work.*    
*†Correspondence and requests for materials should be addressed to J.M. (miao@physics.ucla.edu).*     
 
[If you use any of the data and source codes in your publications and/or presentations, we request that you cite our paper: Y. Yang, J. Zhou, F. Zhu, Y. Yuan, D. Chang, D. S. Kim, M. Pham, A. Rana, X. Tian, Y. Yao, S. Osher, A. K. Schmid, L. Hu, P. Ercius and J. Miao, "Determining the three-dimensional atomic structure of an amorphous solid", _Nature_ **592**, 60–64 (2021).](https://www.nature.com/articles/s41586-021-03354-0)

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)

# Overview

Amorphous solids such as glass are ubiquitous in our daily life and have found broad applications ranging from window glass and solar cells to telecommunications and transformer cores. However, due to the lack of long-range order, the three dimensional (3D) atomic structure of amorphous solids have thus far defied any direct experimental determination. Here, using a glass-forming alloy
as a proof-of-principle, we advanced atomic electron tomography to determine the 3D atomic positions and chemical species in an amorphous solid with a precision of 21 picometer. We quantified the short-range order (SRO) and medium-range order (MRO) of the 3D atomic arrangement and measured the size, shape, volume, and structural distortion of the MROs. The experimental data and source codes for the 3D image reconstruction and post analysis are provided here.

# System Requirements

## Hardware Requirements

We recommend a computer with 16G DRAM, standard i7 4-core CPU, and a GPU to run most data analysis source codes. But for the 3D reconstruction of the experimental data with RESIRE, atomic tracing and the determination of the MROs, we recommend a computer with large memory (256G DRAM, 16-core CPU and 1 GPU).

## Software Requirements

### OS Requirements

This package has been tested on the following Operating System:

Linux: CentOS 6 2.6.32    
Windows: Windows 10 18368.778    
Mac OSX: We have not tested it on a Mac yet, but it should in principle work.     

### Matlab Version Requirements

This package has been tested with `Matlab` R2019b. All the codes have to run in their own folders. We recommend the use of `Matlab` version R2018a or higher to test the data and source codes.

# Repositary Contents

### 1. Experiment Data

Folder: [Measured_data](./1_Measured_data)

This folder contains 55 experimental projections after denoising and alignment as well as their corresponding angles.

### 2. The REal Space Iterative REconstruction (RESIRE) Package

Folder: [RESIRE_package](./2_RESIRE_package)

Run the sample code Main_RESIRE_sample.m to get the 3D reconstruction of a smaller test object. Run the main code `Main_RESIRE_MG.m` to obtain the 3D reconstruction of the multi-component glass-forming sample.

### 3. Reconstructed 3D Volume

Folder: [Final_reconstruction_volume](./3_Final_reconstruction_volume)

This folder contains the 3D volume of the glass-forming nanoparticle reconstructed from `Main_RESIRE_MG.m`.

### 4. Atom Tracing and Classification

Folder: [Tracing_and_classification](./4_Tracing_and_classification)

Run the code `Main_polynomial_tracing.m` to trace the initial atomic positions from the reconstructed 3D volume. After the manual checking of the 3D atomic positions, run the code Main_classification.m to classify the eight elements in the sample into three different types: Co and Ni as type 1, Ru, Rh, Pd and Ag as type 2, and Ir and Pt as type 3.

### 5. Atomic Position Refinement

Folder: [Position_refinement](./5_Position_refinement)

Run the code `Main_position_refinement.m` to refine the 3D atomic coordinates in the glass-forming nanoparticle.

### 6. Experimental Atomic Model

Folder: [Final_coordinates](./6_Final_coordinates)

The final 3D atomic model and chemical species (i.e. type 1, 2 and 3) of the glass-forming nanoparticle.

### 7. Post Data Analysis —— Short Range Order

Folder: [Data_analysis_sro](./7_Data_analysis_sro)

Run the code `Main_1_pdf_and_boo_calculation_all_atoms.m` to calculate the radial distribution function and the bond orientation order parameter for all the atoms in the glass-forming nanoparticle; Run the code `Main_2_pdf_calculation_amorphous_region.m` to compute the radial distribution function and pair distribution function for all the amorphous atoms in the sample; Run the code `Main_3_voronoi_calculation_amorphous_region.m` to determine the Voronoi indices for all the atoms in the sample.

### 8. Post Data Analysis —— Medium Range Order

Folder: [Data_analysis_mro](./8_Data_analysis_mro)

Run the code `Main_1_potential_mro.m` to identify the possible MROs based on the breadth first search algorithm; Run the code `Main_2_final_mro.m` to determine the final MROs in the glass-forming nanoparticle.

### 9. Supplementary Figures

Folder: [Supplementary_Figures](./9_Supplementary_figures)

This file contains four Supplementary Figures.
