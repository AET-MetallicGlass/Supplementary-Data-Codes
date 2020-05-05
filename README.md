# MetallicGlass-Package

Determining the three-dimensional atomic structure of an amorphous solid

Coherent Imaging Group, UCLA, 

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repo Contents](#repo-contents)

# Overview

Amorphous solids such as glass are ubiquitous in our daily life and have found broad applications ranging from window glass and solar cells to telecommunications and transformer cores. However, due to the lack of long-range order, the three dimensional (3D) atomic structure of amorphous solids have thus far defied any direct experimental determination. Here, using a multi-component metallic glass
as a model, we advance atomic electron tomography technique with newly developed RESIRE (REal Space Iterative Algorithm) package to determine its 3D atomicpositions and chemical species with a precision of 21 picometer. We quantify the short-range order (SRO) and medium-range order (MRO) of the 3D atomicarrangement, as well as the size, shape,
volume, and structural distortion of these MROs with unprecedented detail.

# System Requirements

## Hardware Requirements

Most of the AET processing codes require only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with 16G DRAM, standard i7 4-core CPU, and a GPU, which could support running of `RESIRE` package on a projection set with less than 100 pixel size.
Users could check the code `Main_RESIRE_sample.m` in folder `2_RESIRE_package`.

When the matrix size is larger than 100^3 (such as 320^3 for the experimental data here). It is better to use super computer with larger memory (Testing environment: 256G DRAM, 16-core CPU, 1 GPU).

## Software Requirements

### OS Requirements

The package development version is tested on Linux operating systems. The developmental version of the package has been tested on the following systems:

Linux: CentOS 6 2.6.32 
Mac OSX: 
Windows: Windows 10 18368.778 

### Matlab Version Requirements

The package is tested with `Matlab` R2019b. For correctly using of this package, we suggest `Matlab` version R2018a or higher.

# Repo Contents

## 1. Input Experiment Data

Folder:[1_Measured_data](./1_Measured_data)

Denoised and aligned experimental projections by applying the Block-Matching and 3D filtering (BM3D) algorithm for the Metallic Glass sample. Please visit [BM3D package](http://www.cs.tut.fi/~foi/GCF-BM3D/) for more details.

## 2. RESIRE Package

Folder:[2_RESIRE_package](./2_RESIRE_package)

The Real Space Iterative Reconstruction (RESIRE) algorithm package. Run the sample code `Main_RESIRE_sample.m` to see the reconstruction from smaller size sample. Run the main code `Main_RESIRE_MG.m` to get the reconstruction of Metallic Glass nanoparticle.

## 3. Output Experiment Reconstruction Volume

Folder:[3_Final_reconstruction_volume](./3_Final_reconstruction_volume)

The Final reconstruction volume of the Metallic Glass nanoparticle achieved from Main_RESIRE_MG.m.

## 4. Tracing and Classification

Folder:[4_Tracing_and_classification](./4_Tracing_and_classification)

Run the main code `Main_polynomial_tracing.m` to get the initial tracing results from the reconstruction volume. After manually check the peak position, run the main code `Main_classification.m` to get the atomic species results of the Metallic Glass nanoparticle.

## 5. Position Refinement

Folder:[5_Position_refinement](./5_Position_refinement)

Run the main code `Main_position_refinement.m` to get the finalized atomic position of Metallic Glass nanoparticle.

## 6. Output Experimental Atomic Model

Folder:[6_Final_coordinates](./6_Final_coordinates)

The Final atomic model and species of the Metallic Glass nanoparticle.

## 7. Post Data Analysis —— Short Range Order

Folder:[7_Data_analysis_sro](./7_Data_analysis_sro)

Run the code `Main_1_rdf_and_boo_calculation_all_atoms.m` to get the Radial Distribution Function (RDF) and Bond Orientation Order (BOO) for all atoms in the Metallic Glass nanoparticle; Run the code `Main_2_rdf_calculation_amorphous_region.m` to get the Radial Distribution Function (RDF) and Pair Distribution Function (PDF) for amorphous atoms in the Metallic Glass nanoparticle; Run the code `Main_3_voronoi_calculation_amorphous_region.m` to get the Voronoi index for all atoms in the Metallic Glass nanoparticle.

## 8. Post Data Analysis —— Medium Range Order

Folder:[8_Data_analysis_mro](./8_Data_analysis_mro)

Run the code `Main_1_potential_mro.m` to calculate the potential Medium Range Order (MRO) networks based on breadth first search algorithm; Run the code `Main_2_final_mro.m` to get the final MRO networks to fill in the whole nanoparticle.
