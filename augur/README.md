## Augur

_Note: As of Sep 2017, these processing scripts are deprecated in favor of [nextstrain/augur](https://github.com/nextstrain/augur). Current Nextflu builds run off the flu build detailed [here](https://github.com/nextstrain/augur/tree/master/flu). The code in this directory is kept in place for archival reasons._

Augur is the processing pipeline to track flu evolution.  It currently

* imports public sequence data
* subsamples, cleans and aligns sequences
* builds a phylogenetic tree from this data
* reports statistics about mutations and branching patterns of the tree
* infers mutation frequency trajectories through time
* infers antigenic phenotypes from titer data

## Pipeline

The entire pipeline is run with [`process.py`](src/process.py).

### Sequence download, cleaning and alignment

#### Download

Virus sequence data is manually downloaded from the [GISAID EpiFlu database](http://gisaid.org). Data from GISAID may not be disclosed outside the GISAID community. We are mindful of this and raw GISAID data has not been released publicly as part of this project. The current pipeline is designed to work specifically for HA from influenza H3N2. Save GISAID sequences as `data/gisaid_epiflu_sequence.fasta`.

#### [Filter](src/virus_filter.py)

Keeps viruses with fully specified dates, cell passage and only one sequence per strain name. Subsamples to 50 (by default) sequences per month for the last 3 (by default) years before present. Appends geographic metadata. Subsampling prefers longer sequences over shorter sequences and prefer more geographic diversity over less geographic diversity.

#### [Align](src/virus_align.py)

Aligns sequences with [mafft](http://mafft.cbrc.jp/alignment/software/).  Testing showed a much lower memory footprint than [muscle](http://www.drive5.com/muscle/).

#### [Clean](src/virus_clean.py)

Clean up alignment so that reference frame is kept intact. Remove sequences that don't conform to a rough molecular clock and remove known reassortant sequences and other outliers.

### Tree processing

#### [Infer](src/tree_infer.py)

Uses [FastTree](http://meta.microbesonline.org/fasttree/) to get a starting tree, and then refines this tree with [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/).

#### [Refine](src/tree_refine.py)

Reroot the tree based on outgroup strain, collapse nodes with zero-length branches, ladderize the tree and collect strain metadata.

#### [Frequency estimation](src/bernoulli_frequency.py)

Estimate genotype and clade frequency trajectories using a Bernoulli observation model combined with a genetic drift model of process noise.

#### [Streamline](src/streamline.py)

Prep and remove cruft from data files for [auspice](../auspice/) visualization.
