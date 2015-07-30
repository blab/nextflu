## Live website

nextflu is live at [nextflu.org](http://nextflu.org).

## Introduction

nextflu is designed to perform near real-time tracking of influenza virus evolution. It's divided into two components: [augur](augur/), which takes a `.fasta` file of flu sequences and builds an annotated phylogeny, and [auspice](auspice/), which displays this annotated phylogeny in an interactive web-based visualization.

The current version of nextflu is focused on tracking seasonal influenza H3N2 evolution in humans, looking at sequences from the hemagglutinin (HA) gene. Future versions may extend this analysis to other genes in H3N2 or other influenza subtypes. We would also like to implement formal predictive models to make nextflu a platform for forecasting evolution in addition to up-to-date tracking.

## Citation

Please cite nextflu as:

[Neher RA, Bedford T. 2015. nextflu: real-time tracking of seasonal influenza virus evolution in humans. Bioinformatics DOI: 10.1093/bioinformatics/btv381.](http://dx.doi.org/10.1093/bioinformatics/btv381)

## License and copyright

Copyright 2014-2015 Trevor Bedford and Richard Neher.

Source code to nextflu is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextflu is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
