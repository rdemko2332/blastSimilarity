#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process createDatabase {
  output:
    path 'newdb.fasta.*'

  script:
    template 'createDatabase.bash'
}


process blastSimilarity {
  input:
    val fastaName
    path 'subset.fa'
    path blastDatabase

  output:
    path 'blastSimilarity.out'
    path 'blastSimilarity.log'
    path '*.gz*' optional true

  script:
    template 'blastSimilarity.bash'
}


workflow nonConfiguredDatabase {
  take:
    seqs
  main:
    database = createDatabase()
    results = blastSimilarity("newdb.fasta", seqs, database)
    results[0] | collectFile(storeDir: params.outputDir, name: params.dataFile)
    results[1] | collectFile(storeDir: params.outputDir, name: params.logFile)
    results[2] | collectFile(storeDir: params.outputDir)
}


workflow preConfiguredDatabase {
  take:
    seqs
  main:
    database = file(params.databaseDir + "/*")
    results = blastSimilarity(params.databaseBaseName, seqs, database)
    results[0] | collectFile(storeDir: params.outputDir, name: params.dataFile)
    results[1] | collectFile(storeDir: params.outputDir, name: params.logFile)
    results[2] | collectFile(storeDir: params.outputDir)
}