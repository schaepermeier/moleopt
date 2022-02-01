# An R Package Implementing the Multi-Objective Landscape Explorer (MOLE)

The `{moleopt}` package implements an the bi-objective MOLE algorithm. MOLE itself is implemented in C++, and an R interface for simple interoperability with [`{smoof}`](https://github.com/jakobbossek/smoof) functions and visualization tools like [`{moPLOT}`](https://github.com/kerschke/moPLOT) is provided by the package. Further, some tools for visualizing and analyzing the output of MOLE, as well as some simple benchmark functions are implemented.

For now, only the "developer" version is available from GitHub. To install it, either clone this repository or run:

```r
devtools::install_github("schaepermeier/moleopt")
```

Note that some of the documentation is still in-the-works. For now, you may find the [`scripts/experiments.R`](./scripts/experiments.R) useful as a reference on using the package.

If you find the software provided by this repository helpful, please cite it using the metadata provided in CITATION.cff or by:

```bibtex
@software{schaepermeier2022mole,
  author = {Schäpermeier, Lennart},
  month = {2},
  title = {{An R Package Implementing the Multi-Objective Landscape Explorer (MOLE)}},
  url = {https://github.com/schaepermeier/moleopt},
  year = {2022}
}
```
