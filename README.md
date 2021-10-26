# SparsePro for efficient genome-wide fine-mapping with summary statistics and functional annotations

Identifying causal variants from genome-wide association studies (GWASs) is challenging due to widespread linkage disequilibrium (LD). Functional annotations of the genome may help prioritize variants that are biologically relevant and thus improve fine-mapping of GWAS results.

Here, we have developed SparsePro to efficiently conduct functionally informed fine-mapping. Our method has two key features: First, by creating a sparse low-dimensional projection of the high-dimensional genotype, we enable a linear search of causal variants instead of an exponential search of causal configurations in most existing methods; Second, we adopt a probabilistic framework with a highly efficient variational expectation-maximization algorithm to integrate statistical associations and functional priors.

Full description is available at our [preprint paper](https://www.biorxiv.org/content/10.1101/2021.10.04.463133v1). 


## Table of Contents

- [SparsePro for efficient genome-wide fine-mapping with summary statistics and functional annotations](#sparsepro-for-efficient-genome-wide-fine-mapping-with-summary-statistics-and-functional-annotations)
  * [Overview](#overview)
  * [Installation](#installation)
  * [Usage](#usage)
    + [SparsePro-: statistical finemapping with summary statistics](#sparsepro---statistical-finemapping-with-summary-statistics)
    + [Enrich: testing function enrichment of annotations](#enrich--testing-function-enrichment-of-annotations)
    + [SparsePro+: annotated finemapping with summary statistics and functional annotations](#sparsepro---annotated-finemapping-with-summary-statistics-and-functional-annotations)
  * [License](#license)

## Overview 

![](doc/Fig1.png "SparsePro overview")

## Installation

## Usage

### SparsePro-: statistical finemapping with summary statistics

```python
python src/sparsepro.py \
    --ss dat/22.ss \
    --var_Y 1.0 \
    --N 283677 \ 
    --LDdir dat/LD \
    --LDlst dat/22.lst \
    --save dat/res \
    --prefix 22 \
    --tmp True \
    --K 9 
```

### Enrich: testing function enrichment of annotations

```python
python src/enrich.py \
    --save dat/res \
    --prefix 22 \
    --anno dat/22.anno \
    --pip dat/res/22.pip \
    --pthres 0.5
```

### SparsePro+: annotated finemapping with summary statistics and functional annotations 

```python
python src/sparsepro_plus.py \
    --ss dat/22.ss \
    --var_Y 1.0 \ 
    --N 283677 \
    --LDdir dat/LD \
    --LDlst dat/22.lst \
    --save dat/res \
    --prefix 22 \ 
    --K 9 \
    --anno dat/22.anno \
    --W dat/res/22.W0.5 
```

## License

This project is licensed under the MIT License.

## Authors

- Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)
- Hamed Najafabadi (hamed.najafabadi@mcgill.ca)
- Yue Li (yueli@cs.mcgill.ca)

## Citations

If you use this software, please cite:

[Wenmin Zhang, Hamed Najafabadi, Yue Li. SparsePro: an efficient genome-wide fine-mapping method integrating summary statistics and functional annotations. bioRxiv 2021.10.04.463133](https://doi.org/10.1101/2021.10.04.463133)
