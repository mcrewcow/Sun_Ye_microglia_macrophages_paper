# SOCS3/Spp1 axis controls retinal angiogenesis through modulating neovascularization associated microglia 

![image](https://user-images.githubusercontent.com/77118598/223598115-08689a47-3b4c-498d-a02d-7416eceb2c2f.png)
<br />
The datasets are available under the following GEO number: [GSE227861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227861). <br />
For now, as the dataset is private, the token could be given upon request. <br /> 
The analysis consists of the following parts:
1. Seurat datasets preprocessing and annotation, GOBP pathway analysis - file: [Seurat](https://github.com/mcrewcow/Sun_Ye_microglia_macrophages_paper/blob/main/seurat_analysis_P13%5B1%5D.R), [GOBP](https://github.com/mcrewcow/Sun_Ye_microglia_macrophages_paper/blob/main/GOBP_pathway_analysis.R)
2. .loom file generation for RNA Velocity analysis - file: [RNA Velocity preprocessing](https://github.com/mcrewcow/Sun_Ye_microglia_macrophages_paper/blob/main/fastq_to_loom.sh)
3. RNA Velocity analysis + scFates - file: [Jupyter Notebook](https://github.com/mcrewcow/Sun_Ye_microglia_macrophages_paper/blob/main/Socs3-Spp1-no_3D.ipynb), [Python code with 3D](https://github.com/mcrewcow/Sun_Ye_microglia_macrophages_paper/blob/main/Socs3-Spp1.py)
4. CellChat analysis - file: [CellChat](https://github.com/mcrewcow/Sun_Ye_microglia_macrophages_paper/blob/main/cellchat.R)

For questions and suggestions send your correspondence to: <br />
Tianxi Wang - Contributing author - tianxi.wang@childrens.harvard.edu <br />
Guoshuai Cai - Bioinformatical analysis, p.1 - gcai@mailbox.sc.edu <br />
Emil Kriukov - Bioinformatical analysis, p.2-4 - ekriukov@meei.harvard.edu <br />
Ye Sun - Corresponding author - ye.sun@childrens.harvard.edu
