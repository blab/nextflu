# Docker

Previously, we had been using [Docker](https://www.docker.com/) to provide portability of the build process. Currently Nextflu is built with [augur](https://github.com/nextstrain/augur/tree/master/flu) and portability is provided through [janus](https://github.com/nextstrain/janus). Docker functionality is deprecated.

# Command notes

Remove all containers

    docker rm `docker ps --no-trunc -aq`

Remove all images (careful)

    docker rmi $(docker images -q)

# Build

## Build vdb

From within the vdb/ directory

Build vdb docker image

    docker build -t trvrb/vdb:latest .

## Build augur

From within the nextflu/augur/ directory

Build augur docker image

    docker build -t trvrb/augur:latest .

Push image to hub

    docker push trvrb/augur:latest

## Build auspice

From within the nextflu/auspice/ directory

Build auspice docker image

    docker build -t trvrb/auspice:latest .

Push image to hub

    docker push trvrb/auspice:latest

# Run

## Setup

Create a named data volume for augur-data

    docker create --name augur-data -v /Users/trvrb/Dropbox/current-projects/nextstrain/augur/data:/nextflu/augur/data trvrb/augur

Create a named data volume for auspice-data

    docker create --name auspice-data -v /Users/trvrb/Dropbox/current-projects/nextstrain/auspice/data:/nextflu/auspice/data trvrb/auspice

## Augur

From within the nextflu/augur/ directory, run a shell from within the container

    docker run -t -i --volumes-from augur-data --volumes-from auspice-data trvrb/augur /bin/bash

Run augur

    docker run --volumes-from augur-data --volumes-from auspice-data trvrb/augur

Run shell

    docker run -t -i --volumes-from augur-data --volumes-from auspice-data trvrb/augur /bin/bash

## Auspice

From within the nextflu/auspice/ directory, run a shell from within the container

    docker run -t -i --volumes-from auspice-data -p 4000:4000 trvrb/auspice /bin/bash

Run auspice

    docker run --volumes-from auspice-data -p 4000:4000 trvrb/auspice
