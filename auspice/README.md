## Auspice

Website for visualizing and interacting with phylogenies produced by the [augur](../augur/) pipeline.

### Deployment

Website is hosted on Amazon S3. Deploy with:

```
s3_website push --site=.
```

from within this directory. S3 credentials are stored in the PATH as `S3_BUCKET`, `S3_KEY` and `S3_SECRET`. Path can be updated locally without fear of committing private information.
