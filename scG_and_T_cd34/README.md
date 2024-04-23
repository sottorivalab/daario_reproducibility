This directory contains the notebook to reproduce the analysis on the CD34+ cells G&T (Transcriptomics and Methylome) analyisis from [Nam et al.](https://www.nature.com/articles/s41588-022-01179-9). 
To download the data you can go on GEO at this [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158067). We assume the data is placed in a `data` directory.


You can then follow the notebook `analysis_multiome_G_and_T.ipynb` for the actual analysis. I am not sure setting the torch seed is enough to reproduce exactly the analysis also on differen hardware (see [this](https://pytorch.org/docs/stable/notes/randomness.html) for more info), so we share the state dictionary of the exact models we used in the paper in this repo (more info in the noebook)