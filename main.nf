nextflow.enable.dsl=2

process createDatabase {
  output:
    path 'newdb.fasta.*'

  """
  cp $params.databaseFasta newdb.fasta
  makeblastdb \
      -in newdb.fasta \
      -dbtype $params.databaseType	
  """
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

  """
  blastSimilarity.pl \
    --pValCutoff  $params.pValCutoff \
    --lengthCutoff $params.lengthCutoff \
    --percentCutoff $params.percentCutoff \
    --blastProgram $params.blastProgram \
    --database $fastaName \
    --seqFile subset.fa  \
    --blastParams $params.blastParamsFile \
    --doNotParse $params.doNotParse \
    --printSimSeqsFile $params.printSimSeqsFile \
    --saveAllBlastFiles $params.saveAllBlastFiles \
    --saveGoodBlastFiles $params.saveGoodBlastFiles \
    --remMaskedRes $params.adjustMatchLength
  """
}

workflow {
  seqs = Channel.fromPath(params.seqFile).splitFasta( by:params.fastaSubsetSize, file:true  )
  if (!params.preConfiguredDatabase) {
    database = createDatabase()
    results = blastSimilarity("newdb.fasta", seqs, database) 
  }
  if (params.preConfiguredDatabase) {
    database = file(params.databaseDir + "/*")
    results = blastSimilarity(params.databaseBaseName, seqs, database)
  }
  results[0] | collectFile(storeDir: params.outputDir, name: params.dataFile)
  results[1] | collectFile(storeDir: params.outputDir, name: params.logFile)
  results[2] | collectFile(storeDir: params.outputDir)
}

