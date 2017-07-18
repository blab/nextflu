#!/bin/bash
for LINEAGE in H3N2 H1N1pdm Vic Yam
do
  for RESOLUTION in 2y 3y 6y 12y
  do
    echo http://nextflu.org/data/${LINEAGE}_${RESOLUTION}_tree.json
    curl http://nextflu.org/data/${LINEAGE}_${RESOLUTION}_tree.json --compressed -o data/${LINEAGE}_${RESOLUTION}_tree.json
    curl http://nextflu.org/data/${LINEAGE}_${RESOLUTION}_sequences.json --compressed -o data/${LINEAGE}_${RESOLUTION}_sequences.json
    curl http://nextflu.org/data/${LINEAGE}_${RESOLUTION}_frequencies.json --compressed -o data/${LINEAGE}_${RESOLUTION}_frequencies.json
    curl http://nextflu.org/data/${LINEAGE}_${RESOLUTION}_HI.json --compressed -o data/${LINEAGE}_${RESOLUTION}_HI.json
    curl http://nextflu.org/data/${LINEAGE}_${RESOLUTION}_meta.json --compressed -o data/${LINEAGE}_${RESOLUTION}_meta.json
  done
done
