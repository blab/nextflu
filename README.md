## Augur

Augur is Python package to forecast flu evolution.  It will

* import public sequence data
* build a phylogenetic tree from this data
* estimate clade fitnesses
* estimate clade frequencies
* project frequencies foreword

It is intended to be run in an always-on fashion, recomputing predictions daily and pushing predictions to a (static) website.

At the moment, I'm experimenting with [RethinkDB](http://www.rethinkdb.com/) as the database backend.  It's really pleasant to work with.

## Data download

Using [Selenium](https://github.com/SeleniumHQ/selenium) and Python bindings to automate web crawling. [GISAID](http://platform.gisaid.org/epi3/) requires login access.  User credentials are stored in the ENV as `GISAID_USER` and `GISAID_PASS`.

