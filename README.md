# Flex Meta-Storms

![Version](https://img.shields.io/badge/Version-1.0%20-brightgreen)
![Release date](https://img.shields.io/badge/Released%20date-Sept.%2021%2C%202022-brightgreen)



# Contents

- [Introduction](#introduction)
- [System Requirement and dependency](#system-requirement-and-dependency)
- [Installation guide](#installation-guide)
- [Usage](#usage)
- [Example dataset](#example-dataset)
- [Tools in this package](#tools-in-this-package)
- [Supplementary](#supplementary)

  

# Introduction

Flex Meta-Storms (FMS) algorithm implements the “local alignment” concept of microbiome beta diversity for the first time. With microbes of interest (e.g. biomarkers), FMS produces a normalized phylogenetic distance between microbiome pairs. Therefore, FMS can mine potential relationships generated by a small subset of microbes in community samples, elucidating beta diversity with greater sensitivity and flexibility, thereby contributing to a deeper understanding of microbe-host interactions effect. FMS is developed as a standalone package.

# System requirement and dependency

## Hardware requirements

Flex Meta-Storms only requires a standard computer with sufficient RAM to support the operations defined by a user. For typical users, this would be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

  RAM: 8+ GB  
  CPU: 4+ cores

## Software requirements

OpenMP library is the C/C++ parallel computing library. Most Linux releases have OpenMP already been installed in the system. In Mac OS X, to install the compiler that supports OpenMP, we recommend using the Homebrew package manager:
```
brew install gcc
```

## Rscript environment (optional for biomarker selection))

For [biomarkers selection](#b-exact-biomarker-selection-optional) (optional), FMS requires cran-R (http://cran.r-project.org/) 3.2 or higher for the execution of “.R” scripts. Then all packages could be automatically installed and updated by the FMS installer.

# Installation guide

## Automatic installation (recommended)

At present, Flex Meta-Storms provides a fully automatic installer for easy installation.

#### **a. Download the package**

```
git clone https://github.com/qdu-bioinfo/flex-meta-storms.git	
```

#### **b. Install by installer**
```
cd flex-meta-storms
source install.sh
```

The package should take less than 5 minutes to install on a computer with the specifications recommended above.

The example dataset could be found at “example” folder. Check the “example/Readme” for details about the demo run.

## Manual installation

If the automatic installer fails, Flex Meta-Storms can still be installed manually.

#### **a. Download the package**
```
git clone https://github.com/qdu-bioinfo/flex-meta-storms.git	
```

#### **b. Configure the environment variables (the default environment variable configuration file is “~/.bashrc”)**
```
export FlexMetaStorms=Path to Flex Meta-Storms
export PATH=”$PATH:$FlexMetaStorms/bin/”
export PATH=”$PATH:$FlexMetaStorms/Rscript/”
source ~/.bashrc
```
#### **c. Compile the source code**

```
cd flex-meta-storms
make
```
#### **d. Install R packages (optional for [biomarker selection](#b-exact-biomarker-selection-optional))**

```
Rscript $FlexMetaStorms/Rscript/config.R
```

# Usage
#### **a.   Input data formats**  
 FMS requires two files to calculate the “local distances” among microbiomes:
1. Microbial abundance table (e.g. OTU table). Currently FMS supports OTUs of Greengenes (v13-8). More reference database will be released soon. For example  

```
            OTU_1   OTU_2   OTU_3   ...     OTU_M
Sample_1    0.1     0        0.1    ...     0.2
Sample_2    0.2     0.1      0      ...     0.1
...         ...     ...      ...    ...     ...
Sample_N    0       0.3      0.2    ...     0.3
```
2. Exact markers of interest  

```
#Markers
OTU_2
OTU_5
...
OTU_x
```
Biomarkers can either be manually assigned by users in the above format (e.g. parsed from LefSe), or be automatically selected by rank-sum test in the following step.


#### **b. Exact biomarker selection (optional)**  
FMS provides a biomarkers selection tool based on rank-sum test. This biomarker selection requires [microbial abundance table](#a---input-data-formats-required) and metadata in the follow format:
```
            Group
Sample_1    H
Sample_2    D
...         ...
Sample_N    D
```
```
PM_Marker_Test.R -i dataset.abd -m dataset.meta -o Marker
```
Here “dataset.abd” is the microbial abundance table, “dataset.meta” is the metadata, and “Marker” is the output directory of selected biomarkers. In “Marker” directory, the “Out.Group.sig.meanTests.xls” file is the selected biomarkers with significant differences between groups.  


#### **c. Calculate Flex Meta-Storms distance (local alignment distances)**  

With the exact markers, FMS can extract all related community members (target members) using a flexible extraction that considers the weighted taxonomic and functional relations of microbes, and then calculate the normalized phylogeny-based distance as local alignment distance.
```
FMS-comp-taxa -T dataset.abd -M Marker/Out.Group.sig.meanTests.xls -o target.dist
```
The output file “target.dist” is the pairwise matrix of FMS distances. The “-M” assigns the biomarkers, which can either be manually appointed by users (e.g. parsed from LefSe, etc.), or be selected by “[PM_Marker_Test.R](#b-exact-biomarker-selection-optional)” in the FMS package.  



#### **d. Calculate distance on exact markers (optional)**  

The FMS can also calculate the distances ONLY on exact markers (without flexible target member extraction).
```
FMS-comp-taxa -T dataset.abd -M Marker/Out.Group.sig.meanTests.xls -k -o exact.dist
```
The output file “exact.dist” is the pairwise matrix of distances on exact markers. The “-k” is the switch for distances on exact markers.  The “-M” assigns the biomarkers, which can either be manually appointed by users (e.g. parsed from LefSe, etc.), or be selected by “[PM_Marker_Test.R](#b-exact-biomarker-selection-optional)” in the FMS package. 


# Example dataset
Here, we provide a demo dataset in the "example" folder with 136 real microbiomes. In this package, "dataset.abd" is the relative abundance of OTU table, and "dataset.meta" is the metadata of samples. To run the demo, you can either automatically start:
To run the demo, you can either:
```
cd example
sh Readme
```
or type the following command to calculate the Flex Meta-Storms distance:

```
#Biomarker selection for samples
PM_Marker_Test.R -m dataset.meta -i dataset.abd -o Marker

#Calculate the Flex Meta-Storms distance of the samples
FMS-comp-taxa -T dataset.abd -M ./Marker/Out.Group.sig.meanTests.xls -o target.dist

#Calculate the Meta-Storms distance on exact markers of the samples
FMS-comp-taxa -T dataset.abd -M ./Marker/Out.Group.sig.meanTests.xls -k -o exact.dist
```
### Output description
We also provide the example [output](https://github.com/qdu-bioinfo/flex-meta-storms/tree/master/example/example_output)  
target.dist: Regular pairwise Flex Meta-Storms distance of 136 samples  
target.dist.target_marker: Relative abundance of target members in 136 samples  
exact.dist: Regular pairwise Meta-Storms distance on exact markers of 136 samples  
exact.dist.exact_marker: Relative abundance of exact markers in 136 samples

# Tools in this package
#### **a. PM_Marker_Test.R**

Biomarker selection from a microbial feature table. Run:
```
PM_Marker_Test.R -h
```
for detailed parameters.


#### **b. FMS-comp-taxa**

Calculate Flex Meta-Storms distance between samples. Run:
```
FMS-comp-taxa -h 
```
for detailed parameters.



# Supplementary

[Artificial Dataset](http://bioinfo.single-cell.cn/Released_Software/flex-meta-storms/data/Artificial_dataset.tar.gz) contains 100 simulated microbiomes from Greengenes (v13-8) OTUs.  

[Real Dataset 1](http://bioinfo.single-cell.cn/Released_Software/flex-meta-storms/data/Real_dataset_I.tar.gz) contains 88 ASD 16S amplicon samples processed by Parallel-Meta Suite.

[Real Dataset 2](http://bioinfo.single-cell.cn/Released_Software/flex-meta-storms/data/Real_dataset_II.tar.gz) Contains 104 CRC 16S amplicon samples processed by Parallel-Meta Suite.
