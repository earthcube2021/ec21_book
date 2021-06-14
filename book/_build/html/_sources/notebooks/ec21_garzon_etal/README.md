# Exploring Research Data Repositories with geoextent.
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/earthcube2021/ec21_garzon_etal/HEAD)

This notebook presents [geoextent](https://github.com/o2r-project/geoextent), a Python library for reliably extracting the geospatial and temporal extents of files, directories, and repository records. 

## Showcase

To run the showcase notebook, install [JupyterLab](https://jupyter.org/) or the classic Jupyter Notebook and then start a local server as shown below.
If your IDE has support for the Jupyter format, installing `ipykernel` might be enough.
We recommend running the below commands in a virtual environment as described [here](https://jupyter-tutorial.readthedocs.io/en/latest/first-steps/install.html).
The notebook must be [trusted](https://jupyter-notebook.readthedocs.io/en/stable/security.html#notebook-security).

```bash
cd showcase
pip install -r requirements.txt
pip install -r showcase/requirements.txt
pip install -e .

jupyter trust showcase/SG_01_Exploring_Research_Data_Repositories_with_geoextent.ipynb
jupyter lab
```

Then open the local Jupyter Notebook server using the displayed link and open the notebook (`*.ipynb` files) in the `showcase/` directory.
Consult the documentation on [paired notebooks based on Jupytext](https://github.com/mwouts/jupytext/blob/master/docs/paired-notebooks.md) before editing selected notebooks.

## How to cite

> Garzón, Sebastian and Nüst, Daniel, 2021. Exploring Research Data Repositories with geoextent.

## License

`geoextent` is licensed under MIT license, see file [LICENSE](https://github.com/earthcube2021/ec21_garzon_etal/blob/master/LICENSE).

This project is developed as part of the [DFG-funded](https://o2r.info/about/#funding) research project Opening Reproducible Research (o2r, [https://o2r.info](https://o2r.info)) under project number PE 1632/17-1.

Copyright (C) 2021 - o2r project.
