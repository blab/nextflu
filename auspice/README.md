## Auspice

Website for visualizing and interacting with phylogenies produced by the [nextstrain/augur](https://github.com/nextstrain/augur) pipeline. Interactive visualization is done with [d3.js](http://d3js.org/). This is v1 of auspice and is in maintenance mode. Active development is continuing on v2 of auspice at [nextstrain/auspice](https://github.com/nextstrain/auspice).

### Build and compile

Build JSONs via augur. These will be exported to `augur/flu/auspice/`. Move JSON files to `auspice/data/`. Compile site with [Jekyll](http://jekyllrb.com/) by running `jekyll build` from within the `auspice/` directory. This creates a `_site/` directory containing compiled resources.

### Deployment

Website is hosted on Amazon S3. Deploy with `s3_website push` from within the `auspice/` directory. This pushes the `_site/` directory to the specified S3 bucket. S3 credentials are stored in ENV as `S3_BUCKET`, `S3_KEY` and `S3_SECRET`. Environment variables can be updated locally without fear of committing private information.

### Local development

Develop locally by running `jekyll serve` and then going to `http://localhost:4000/` in the browser. Jekyll will recompile as local files are modified.
