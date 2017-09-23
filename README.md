## Live website

nextflu is live at [nextflu.org](https://nextflu.org).

## Introduction

nextflu is designed to perform near real-time tracking of influenza virus evolution. The current version of nextflu is focused on tracking seasonal influenza evolution in humans, looking at sequences from the hemagglutinin (HA) gene. It's divided into two components:
* _augur_, which takes a `.fasta` file of flu sequences and builds an annotated phylogeny
* _auspice_, which displays this annotated phylogeny in an interactive web-based visualization
Augur build scripts are housed at [nextstrain/augur](https://github.com/nextstrain/augur). These produce a series of JSON files that are displayed interactively on the web. Currently, the same JSONs produced by augur can be displayed with auspice v1 housed in [this repo](auspice/) or with auspice v2 housed at [nextstrain/auspice](https://github.com/nextstrain/auspice). Auspice v1 is live at [nextflu.org](https://nextflu.org) and auspice v2 is live at [nextstrain.org/flu](http://nextstrain.org/flu). Auspice v1 still provides greater functionality for influenza, but the intention is to eventually migrate to auspice v2.

## Citation

Please cite nextflu as:

[Neher RA, Bedford T. 2015. nextflu: real-time tracking of seasonal influenza virus evolution in humans. Bioinformatics DOI: 10.1093/bioinformatics/btv381.](http://dx.doi.org/10.1093/bioinformatics/btv381)

## License and copyright

Copyright 2014-2017 Trevor Bedford and Richard Neher.

Source code to nextflu is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextflu is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
