# HHTAD
HHTAD: Topologically Associating Domain Detection Method based Higher-Order Reads and Hypergraph
# About HHTAD
A weighted hypergraph is constructed from higher-order reads, and random walks are performed on the hypergraph to obtain a feature representation for each genomic bin. Multi-scale convolution and a bidirectional GRU are jointly employed to extract local and long-range correlation features. The hypergraph flux insulation score, combined with the directionality score and intra-cluster sharpness, is used to quantify boundary insulation capability and filter out false-positive boundaries. Finally, hypergraph random walks and cosine similarity are integrated for iterative clustering, yielding the final set of nested TADs.
# Requirements
tensorflow =2.15.0 ，numpy =1.26.4 ，pandas= 2.3.3 ，scikit-learn=1.7.2
# Usage
## Source of the Datasets
The HiPore-C datasets analysed in HHTAD are available from the Gene Expression Omnibus (GEO) repository under accession number GSE202539 (link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202539). The Hi-C datasets can be accessed in the Zenodo repository at https://zenodo.org/records/10822184. ChIP-seq data for the human GM12878 and K562 cell lines were downloaded from the ENCODE portal (link: https://www.encodeproject.org/).

| Cell lines | Datasets | Brief description | Data source |
| --- | --- | --- | --- |
| GM12878 | HiPore-C | G_FC5 | GSE202539 |
| GM12878 | Hi-C | GM12878_FC5.mcool | https://zenodo.org/records/10822184/files/GM12878_FC5.mcool?download=1 |
| K562 | HiPore-C | K_FC2 | GSE202539 |
| K562 | Hi-C | K562_FC2.mcool | https://zenodo.org/records/10822184/files/K562_FC2.mcool?download=1 |

| Species | Cell lines | Name | Download information |
| --- | --- | --- | --- |
| Human | GM12878 | CTCF | ENCFF796WRU |
| Human | GM12878 | SMC3 | ENCFF572RPI |
| Human | GM12878 | RAD21 | ENCFF004RJL |
| Human | GM12878 | H3K36me3 | ENCFF432EMI |
| Human | GM12878 | H3K4me3 | ENCFF587DVA |
| Human | GM12878 | H3K9me3 | ENCFF874UEV |
| Human | GM12878 | Polymerase II | ENCFF912DZY |
| Human | GM12878 | SINEs | http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/msk.txt.gz |
| Human | K562 | CTCF | ENCFF973WXW |
| Human | K562 | SMC3 | ENCFF151QIM |
| Human | K562 | RAD21 | ENCFF265WWW |
| Human | K562 | H3K36me3 | ENCFF119BYM |
| Human | K562 | H3K4me3 | ENCFF112BHN |
| Human | K562 | H3K9me3 | ENCFF371GMJ |

## First --Model training
Model training

```python
python main.ipynb
```

## Second --Model prediction
#predicted
```python
load_model.ipynb
```
## Third --Filtering false-positive borders and assembling merged TADs
1. hypergraph flux coefficient calculation
```python
python flux.py
```

2. Filtering false-positive borders based hypergraph flux coefficient
```python
python fliter_borders.py
```
3. Nested TADs are obtained

```python
python merge_TAD.py
```
