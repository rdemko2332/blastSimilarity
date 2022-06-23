nextflow.enable.dsl=2

process createDatabase {
  output:
  path 'newdb.fasta.*'

  """
  createDatabase.pl --dbFile $params.databaseFasta --blastProgram $params.blastProgram 
  """
}

process blastSimilarity {
  publishDir params.outputDir, mode: 'copy', saveAs: {filename -> filename.endsWith(".out") ? params.dataFile : filename.endsWith(".log") ? params.logFile : filename }

  input:
  val fastaName
  path 'subset.fa'
  // Retaining names from the input into db_vch
  path blastDatabase

  output:
  path 'blastSimilarity.out'
  path 'blastSimilarity.log'
  path '*.gz*' optional true

  """
  blastSimilarity --pValCutoff  $params.pValCutoff --lengthCutoff $params.lengthCutoff --percentCutoff  $params.percentCutoff --blastProgram  $params.blastProgram --database $fastaName --seqFile subset.fa  --blastParams $params.blastParams --doNotParse $params.doNotParse --printSimSeqsFile $params.printSimSeqsFile --saveAllBlastFiles $params.saveAllBlastFiles --saveGoodBlastFiles $params.saveGoodBlastFiles --remMaskedRes $params.adjustMatchLength
  """
}

workflow {
  if (params.preConfiguredDatabase == "false") {
    database = createDatabase()
    seqs = Channel.fromPath(params.seqFile).splitFasta( by:1, file:true  )
    blastSimilarity("newdb.fasta", seqs, database)
  }
  if (params.preConfiguredDatabase == "true") {
    database = file(params.databaseDir + "/*")
    seqs = Channel.fromPath(params.seqFile).splitFasta( by:1, file:true  )
    blastSimilarity(params.databaseBaseName, seqs, database)
  }
}

