## Auspice

_Note: As of Sep 2017, this JavaScript app is deprecated in favor of [nextstrain/auspice](https://github.com/nextstrain/auspice). The code in this directory is kept in place for archival reasons._

Website for visualizing and interacting with phylogenies produced augur pipeline. Interactive visualization is done with [d3.js](http://d3js.org/).

### Build and compile

1. Build JSONs via augur. These will be exported to `augur/flu/auspice/`. Move JSON files to `auspice/data/`.
2. Create index files by running `python provision_directories.py` from `auspice/`.
3. Compile site with [Jekyll](http://jekyllrb.com/) by running `jekyll build` from within the `auspice/` directory. This creates a `_site/` directory containing compiled resources.

### Deployment

Website is hosted on Amazon S3. Deploy with `s3_website push` from within the `auspice/` directory. This pushes the `_site/` directory to the specified S3 bucket. S3 credentials are stored in ENV as `S3_BUCKET`, `S3_KEY` and `S3_SECRET`. Environment variables can be updated locally without fear of committing private information.

### Local development

Develop locally by running `jekyll serve` and then going to `http://localhost:4000/` in the browser. Jekyll will recompile as local files are modified.
