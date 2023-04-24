# Analysis and visualize code of *P.sibiricum* research.

[![R version](https://img.shields.io/badge/R-v4.2.2-salmon)](https://www.r-project.org) [![license](https://img.shields.io/badge/license-MIT-cyan)](https://opensource.org/licenses/MIT) [![Myblog](https://img.shields.io/badge/Blog-ShanwLearnBioinfo-purple)](https://shawnwx2019.github.io/) [![DOI](https://img.shields.io/badge/DOI-10.3389/fpls.2023.1181861-blue)](https://www.frontiersin.org/articles/10.3389/fpls.2023.1181861/full) [![Accept](https://img.shields.io/badge/Accept%20date-18%20Apr%202023-orange)](https://www.frontiersin.org/articles/10.3389/fpls.2023.1181861/full)

# Citation

Ou X<sup>#\*</sup>., Wang X<sup>#</sup>., Zhao B., Zhao Y., Liu HQ., Chang Y., Wang Z., Yang W., Zhang X<sup>\*</sup>., Yu K<sup>\*</sup>. Metabolome and transcriptome signatures shed light on the anti-obesity effect of *Polygonatum sibiricum*. 2023, *Frontiers in Plant Science.*, 14:1273.

# Metabolomics data analysis

Metabolomics data was processed using R packages: [TidyMass](https://github.com/tidymass/tidymass).

From peak picking to metabolite annotation was completed using the [HPC-TidyMass pipeline](https://github.com/ShawnWx2019/HPC_tidymass/tree/main) , in which only the peak picking step was utilized.

``` shell
## The help information of Tidymass pipeline
###############################################################################################
/home/data/shawn/01.src/02.script/02.Tidymass/runTidymass.sh: option requires an argument -- h
Tidymass pipeline part1. Date transform ,peak picking and annotation
Usage:
-------------------------------
runTidymass [-i input] [-t type] [-c column]
-------------------------------
Author Shawn Wang (shawnwang2016@126.com)
-------------------------------

Date Tus Feb 09, 2023
-------------------------------

Version V.0.0.0.99 beta
-------------------------------

Description
-------------------------------

[-i]:input,   The file path of .raw data
[-t]:type,  The type of ion model, 1: NEG+POS in one file, 2: NEG and POS in differnet files
[-c]:column,  rp or hilic
###############################################################################################

bash runTidymass.sh -i raw/p_sibiricum -t 1 -c rp
```

# Content

**Step1. convert .raw data to .mzXML and .mgf**

-   .raw 2 .mzXML: [01.msconvert.sh](https://github.com/ShawnWx2019/HPC_tidymass/blob/main/src/shell/01.msconvert.sh)

-   .raw 2 .mgf: [03.ms2convert.sh](https://github.com/ShawnWx2019/HPC_tidymass/blob/main/src/shell/03.ms2convert.sh)

**Step2. Data cleaning and normalization**

-   Raw data overview: [OverviewMetabolomicsData.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/02.Metabolomics/OverviewMetabolomicsData.R)

-   Data cleaning and normalization: [DataCleaning.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/02.Metabolomics/DataCleaning.R)

-   Data annotation: [DataAnnotation.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/02.Metabolomics/DataAnnotation.R)

-   Remove redundancy and classification: [RemoveRedundancyandClassification.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/02.Metabolomics/RemoveRedundancyandClassification.R)

-   Metabolomics Mfuzz: [MetaMfuzz.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/02.Metabolomics/MetaMfuzz.R)

**Step3. Pacbio visualization**

-   Transcript annotation integration: [PacbioVisualization.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/03.PacBio/PacbioVisualization.R)

-   Picbio BUSCO analysis: [BUSCO.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/03.PacBio/BUSCO.R)

**Step4. NGS data analysis**

-   DEG analysis and visualization, (hclust, 3D pca, DEG, mFuzz): [DEG_analysis.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/04.Transcriptome/DEG_analysis.R)

-   PSP pathway: [PSP.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/04.Transcriptome/PSP.R)

-   unFA pathway: [unFA.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/04.Transcriptome/unFA.R)

-   phloretin pathway: [phloretin.R](https://github.com/ShawnWx2019/P.sibiricum_code/blob/main/04.Transcriptome/phloretin.R)
