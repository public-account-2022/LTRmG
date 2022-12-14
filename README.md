# Lineage Tree Reconstruction with macsGESTALT

## Background

Single-cell lineage tracing technology enables an unprecedented exploration of lineage history and transcriptional states. By tracking clonal relationships at the single-cell level, this method has led to a new level of understanding of several biological processes, such as clonal evolution (Wagner *et al.*, 2020).

macsGASTALT (multiplexed, activatable, clonal and subclonal GESTALT) is an inducible CRISPR-Cas9-based single-cell lineage tracing method. Each cell is cotransfected with three genetic constructs, dox-inducible Cas9, the gRNA array and multiplexed barcodes. The gRNA array consists of 5 differ-ent gRNA separated by tRNA sequences, and each barcode consists of a static 10bp region and 250bp evolving region that contains 5 CRISPR-Cas9 target sites. After the doxycycline induction, Cas9 expression is triggered, and the individual gRNAs that are randomly released from the gRNA array can guide the Cas9 to edit the evolving regions. As a result, transfected cells will contain a random number of integrated barcodes, and these cells are used for in vivo or in vitro studies (Simeonov *et al.*, 2021).

In this method, the static barcodes in each cell remain unchanged throughout the development, which enables the discovery of origin clones. With the inducible CRISPR-Cas9 system, the evolving region of barcodes changed progressively (Simeonov *et al.*, 2021). During the continuous cell division, the evolving regions are edited and inherited, thus providing phylogenetic or subclonal information for bioinformatic analyses (Simeonov *et al.*, 2021). Based on static and dynamic barcode information generalized from the macsGASTALT technique and single cell sequencing, our group designed a pipeline that reconstructs the phylogenetic trees of lineage history with the identification of cell clones and subclones.

![](github-figs/Original.png)

## Analysis Workflow

The analysis workflow of this algorithm project can be divided into four steps.

- Cluster cells into clones by static barcode ID.
- Align evolving barcodes.
- Cluster cells from the same clone into subclones.
- Visualize the results.

![](github-figs/workflow.jpg)

## Install

You can directly download this package from GitHub.

## Setting up running environment

Note that we use Python 3.8.15. 

| Package    | Version |
|------------|---------|
| circlify   | 0.15.0  |
| graphiz    | 0.20.1  |
| igraph     | 0.10.2  |
| matplotlib | 3.6.2   |
| networkx   | 2.8.8   |
| numpy      | 1.23.5  |
| pandas     | 1.5.2   |
| pygraphiz  | 1.10    |
| seaborn    | 0.12.1  |

## Input

Barcode datasets sequenced from macsGESTALT technique, which should includes information such as cell barcode IDs, static barcode IDs, evolving barcode sequences (1-5).

## Parameters

Eight parameters in total can be specified when running this package in the terminal environment.

| Parameters                | Explanation                                                                                                                                                                                                                                                   |
|-------------------|-----------------------------------------------------|
| **-f/--filename**         | a file path for input, usually is a file contains cell IDs and barcode sequences.                                                                                                                                                                             |
| **-o/--output**           | a file path for output, including a .txt file and figures in .png format, stored in separate directories under the current working directory.                                                                                                                 |
| **-s/--sample**           | a table file contains two columns that describes harvest sites and samples. ’Tissue’ indicates the harvest sites from where cells are extracted, ’sample’ indicates whether cells from each site should be sampled when identifying clones, indicating by True or False. |
| **-n/--numclus** | the number of cell clones to call during clone clustering.                                                                                                                                                                  |
| **-a/--score1**           | score of match in alignment score calculation process.                                                                                                                                                                                                        |
| **-b/--score2**           | score of mismatch in alignment score calculation process.                                                                                                                                                                                                     |
| **-c/--score3**           | score of indel in alignment score calculation process.                                                                                                                                                                                                        |
| **-p/--position**           | the name (in column "sample") of tumor initial site.

## Usage

You can run the following code to test our package after downloading it. 

``` bash
python3 LTRmG.py -f example.txt \
-o ./results/subclones.txt \
-a 1 -b -2 -c -1
```

## Output

A *.txt* file contains a list describes clone and subclone information of cells, which can be used for visualization, and also some figures that visualize clones and subclones information. Stored in the different directories under the current working directory.

## Example results

We randomly select cells from the original barcode datasets provided by the paper and form testing datasets contains 20, 30 and 60 cells, respectively. Those datasets are provided in the package for test running. Below are some examples results of visualization, representing analyses results of 60 (left) or 30 (right) cells, respectively. 

<img src="github-figs/circle.jpg" width="400"/><img src="github-figs/circletree.jpg" width="400"/>


## References

Simeonov, K. P., Byrns, C. N., Clark, M. L., Norgard, R. J., Martin, B., Stanger, B. Z., Shendure, J., McKenna, A., & Lengner, C. J. (2021). Single-cell lineage tracing of metastatic cancer reveals selection of hybrid EMT states. Cancer cell, 39(8), 1150--1162.e9. <https://doi.org/10.1016/j.ccell.2021.05.005>

Wagner, D. E., & Klein, A. M. (2020). Lineage tracing meets single-cell omics: opportunities and challenges. Nature reviews. Genetics, 21(7), 410--427. <https://doi.org/10.1038/s41576-020-0223-2>

## Maintainers

Zhejiang University-University of Edinburgh Institute (ZJE)

Biomedical Informatics Class 2020

BMI3 Course: Group 4, Project 2

Supervisor: Wanlu Liu

## Contributing

Feel free to contact us and join in!

## License

[MIT](LICENSE)
