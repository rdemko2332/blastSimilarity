#!/usr/bin/env bash

cp $params.databaseFasta newdb.fasta
makeblastdb \
  -in newdb.fasta \
  -dbtype $params.databaseType
