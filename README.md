# Kerr Circular Equatorial orbits: Vaccum and non vacuum 
Please refer and cite:  
https://arxiv.org/abs/2410.17310

### Installing
Create a Conda env:

```
conda create -n few_Kerr_NonVac -c conda-forge gcc_linux-64 gxx_linux-64 wget gsl lapack=3.6.1 hdf5 numpy Cython scipy tqdm jupyter ipython h5py requests matplotlib python
```


Clone the repository:
```
git clone git@github.com:Hassankh92/FastEMRIWaveforms_KerrCircNonvac.git
cd FastEMRIWaveforms_KerrCircNonvac
```


```
python setup.py install
```


Run the tests: 

For trajectories
```
Kerr_vac_nonvac_test.py
```

For Waveforms test, refer to the Test_Notebook.ipynb







### Feel free to email me if you have any questions!


# few: Fast EMRI Waveforms

This package contains the highly modular framework for fast and accurate extreme mass ratio inspiral (EMRI) waveforms from [arxiv.org/2104.04582](https://arxiv.org/abs/2104.04582) and [arxiv.org/2008.06071](https://arxiv.org/abs/2008.06071). The waveforms in this package combine a variety of separately accessible modules to form EMRI waveforms on both CPUs and GPUs. Generally, the modules fall into four categories: trajectory, amplitudes, summation, and utilities. Please see the [documentation](https://bhptoolkit.org/FastEMRIWaveforms/) for further information on these modules. The code can be found on Github [here](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms). The data necessary for various modules in this package will automatically download the first time it is needed. If you would like to view the data, it can be found on [Zenodo](https://zenodo.org/record/3981654#.XzS_KRNKjlw). The current and all past code release zip files can also be found on Zenodo [here](https://zenodo.org/record/8190418). Please see the [citation](#citation) section below for information on citing FEW.

This package is a part of the [Black Hole Perturbation Toolkit](https://bhptoolkit.org/).

If you use all or any parts of this code, please cite [arxiv.org/2104.04582](https://arxiv.org/abs/2104.04582) and [arxiv.org/2008.06071](https://arxiv.org/abs/2008.06071). See the [documentation](https://bhptoolkit.org/FastEMRIWaveforms/) to properly cite specific modules.




## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

If you want to develop FEW and produce documentation, install `few` with
```
bash install.sh install_type=development
```
This will install necessary packages for building the documentation (`sphinx`, `pypandoc`, `sphinx_rtd_theme`, `nbsphinx`). The documentation source files are in `docs/source`. To compile the documentation, change to the `docs` directory and run `make html`. 

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms/tags).

Current Version: 1.5.5

## Authors/Developers

* **Michael Katz**
* Lorenzo Speri
* Christian Chapman-Bird
* Alvin J. K. Chua
* Niels Warburton
* Scott Hughes
* Hassan Khalvati

### Contibutors

* Philip Lynch
* Soichiro Isoyama
* Ryuichi Fujita
* Monica Rizzo

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details.

## Citation

Please make sure to cite FEW papers and the FEW software on [Zenodo](https://zenodo.org/record/8190418). There are other papers that require citation based on the classes used. For most classes this applies to, you can find these by checking the `citation` attribute for that class. Below is a list of citable papers that have lead to the development of FEW. 

```
@article{Chua:2020stf,
    author = "Chua, Alvin J. K. and Katz, Michael L. and Warburton, Niels and Hughes, Scott A.",
    title = "{Rapid generation of fully relativistic extreme-mass-ratio-inspiral waveform templates for LISA data analysis}",
    eprint = "2008.06071",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevLett.126.051102",
    journal = "Phys. Rev. Lett.",
    volume = "126",
    number = "5",
    pages = "051102",
    year = "2021"
}

@article{Katz:2021yft,
    author = "Katz, Michael L. and Chua, Alvin J. K. and Speri, Lorenzo and Warburton, Niels and Hughes, Scott A.",
    title = "{FastEMRIWaveforms: New tools for millihertz gravitational-wave data analysis}",
    eprint = "2104.04582",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "4",
    year = "2021"
}

@software{michael_l_katz_2023_8190418,
  author       = {Michael L. Katz and
                  Lorenzo Speri and
                  Alvin J. K. Chua and
                  Christian E. A. Chapman-Bird and
                  Niels Warburton and
                  Scott A. Hughes},
  title        = {{BlackHolePerturbationToolkit/FastEMRIWaveforms: 
                   Frequency Domain Waveform Added!}},
  month        = jul,
  year         = 2023,
  publisher    = {Zenodo},
  version      = {v1.5.1},
  doi          = {10.5281/zenodo.8190418},
  url          = {https://doi.org/10.5281/zenodo.8190418}
}

@article{Chua:2018woh,
    author = "Chua, Alvin J.K. and Galley, Chad R. and Vallisneri, Michele",
    title = "{Reduced-order modeling with artificial neurons for gravitational-wave inference}",
    eprint = "1811.05491",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.IM",
    doi = "10.1103/PhysRevLett.122.211101",
    journal = "Phys. Rev. Lett.",
    volume = "122",
    number = "21",
    pages = "211101",
    year = "2019"
}

@article{Fujita:2020zxe,
    author = "Fujita, Ryuichi and Shibata, Masaru",
    title = "{Extreme mass ratio inspirals on the equatorial plane in the adiabatic order}",
    eprint = "2008.13554",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.102.064005",
    journal = "Phys. Rev. D",
    volume = "102",
    number = "6",
    pages = "064005",
    year = "2020"
}

@article{Stein:2019buj,
    author = "Stein, Leo C. and Warburton, Niels",
    title = "{Location of the last stable orbit in Kerr spacetime}",
    eprint = "1912.07609",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.101.064007",
    journal = "Phys. Rev. D",
    volume = "101",
    number = "6",
    pages = "064007",
    year = "2020"
}

@article{Chua:2015mua,
    author = "Chua, Alvin J.K. and Gair, Jonathan R.",
    title = "{Improved analytic extreme-mass-ratio inspiral model for scoping out eLISA data analysis}",
    eprint = "1510.06245",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1088/0264-9381/32/23/232002",
    journal = "Class. Quant. Grav.",
    volume = "32",
    pages = "232002",
    year = "2015"
}

@article{Chua:2017ujo,
    author = "Chua, Alvin J.K. and Moore, Christopher J. and Gair, Jonathan R.",
    title = "{Augmented kludge waveforms for detecting extreme-mass-ratio inspirals}",
    eprint = "1705.04259",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.96.044005",
    journal = "Phys. Rev. D",
    volume = "96",
    number = "4",
    pages = "044005",
    year = "2017"
}

@article{Barack:2003fp,
    author = "Barack, Leor and Cutler, Curt",
    title = "{LISA capture sources: Approximate waveforms, signal-to-noise ratios, and parameter estimation accuracy}",
    eprint = "gr-qc/0310125",
    archivePrefix = "arXiv",
    doi = "10.1103/PhysRevD.69.082005",
    journal = "Phys. Rev. D",
    volume = "69",
    pages = "082005",
    year = "2004"
}

@article{Speri:2023jte,
    author = "Speri, Lorenzo and Katz, Michael L. and Chua, Alvin J. K. and Hughes, Scott A. and Warburton, Niels and Thompson, Jonathan E. and Chapman-Bird, Christian E. A. and Gair, Jonathan R.",
    title = "{Fast and Fourier: Extreme Mass Ratio Inspiral Waveforms in the Frequency Domain}",
    eprint = "2307.12585",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "7",
    year = "2023"
}
```

## Acknowledgments

* This research resulting in this code was supported by National Science Foundation under grant DGE-0948017 and the Chateaubriand Fellowship from the Office for Science \& Technology of the Embassy of France in the United States.
* It was also supported in part through the computational resources and staff contributions provided for the Quest/Grail high performance computing facility at Northwestern University.
