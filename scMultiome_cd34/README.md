This directory contains the notebook to reproduce the analysis on the CD34+ cells multiome analyisis from [Setty et al.](https://www.nature.com/articles/s41587-019-0068-4). 
To download the data you can just run 

```bash
mkdir data/
wget https://dp-lab-data-public.s3.amazonaws.com/SEACells-multiome/cd34_multiome_rna.h5ad -O data/cd34_multiome_rna.h5ad  # RNA
wget https://dp-lab-data-public.s3.amazonaws.com/SEACells-multiome/cd34_multiome_atac.h5ad -O data/cd34_multiome_atac.h5ad # ATAC
```

You can then follow the notebook `Analysis_multiome_10x.ipynb` for the actual analysis. I am not sure setting the torch seed is enough to reproduce exactly the analysis also on differen hardware (see [this](https://pytorch.org/docs/stable/notes/randomness.html) for more info), so we share the state dictionary of the exact models we used in the paper in this repo (more info in the noebook)
