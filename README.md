## Current status

_All Nextflu development has been folded into the Nextstrain project. Python build pipine is available at [nextstrain/augur](https://github.com/nextstrain/augur) (which supercedes build pipeline in `blab/nextflu/augur`) and JavaScript visualization is available at [nextstrain/auspice](https://github.com/nextstrain/auspice) (which supercedes visualization in `blab/nextflu/auspice`)._

_Live display of seasonal flu evolution is available at [nextstrain.org/flu](http://nextstrain.org/flu)._

## Introduction

nextflu is designed to perform near real-time tracking of influenza virus evolution. The current version of nextflu is focused on tracking seasonal influenza evolution in humans, looking at sequences from the hemagglutinin (HA) gene. It's divided into two components:
* _augur_, which takes a `.fasta` file of flu sequences and builds an annotated phylogeny
* _auspice_, which displays this annotated phylogeny in an interactive web-based visualization

Augur build scripts are housed at [nextstrain/augur](https://github.com/nextstrain/augur) and auspice visualization tool is housed at [nextstrain/auspice](https://github.com/nextstrain/auspice). Legacy v1 versions of both `augur` and `auspice` are maintained in this repository for purposes of posterity.

## Citation

Please cite nextflu as:

[Neher RA, Bedford T. 2015. nextflu: real-time tracking of seasonal influenza virus evolution in humans. Bioinformatics DOI: 10.1093/bioinformatics/btv381.](http://dx.doi.org/10.1093/bioinformatics/btv381)

## License and copyright

Copyright 2014-2017 Trevor Bedford and Richard Neher.

Source code to nextflu is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextflu is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
