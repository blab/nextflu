## Auspice

Website for visualizing and interacting with phylogenies produced by the [augur](../augur/) pipeline.

### Compiling

Compile site with [Jekyll](http://jekyllrb.com/) by running `jekyll build` from within the `auspice/` directory. This creates a `_site/` directory containing compiled resources.

### Deployment

Website is hosted on Amazon S3. Deploy with `s3_website push` from within the `auspice/` directory. This pushes the `_site/` directory to the specified S3 bucket. S3 credentials are stored in ENV as `S3_BUCKET`, `S3_KEY` and `S3_SECRET`. Environment variables can be updated locally without fear of committing private information.

### Local development

Develop locally by running `jekyll serve` and then going to `http://localhost:4000/` in the browser. Jekyll will recompile as local files are modified.
