nextflow.enable.dsl=1
seq_qch = Channel.fromPath(params.seqFile).splitFasta( by:1, file:true  )
db_vch = Channel.value()
process createDatabase {
    output:
    path 'newdb.fasta.*' into db_vch
    """
    createDatabase.pl --dbFile $params.databaseFasta --blastProgram $params.blastProgram 
    """
}

process blastSimilarity {
   input:
   path 'subset.fa' from seq_qch
   // Retaining names from the input into db_vch
   path blastDatabase from db_vch
   output:
   path 'blastSimilarity.out' into output_qch
   path 'blastSimilarity.log' into log_qch
   path '*.gz*' optional true into zip_qch
   """
   blastSimilarity --pValCutoff  $params.pValCutoff --lengthCutoff $params.lengthCutoff --percentCutoff  $params.percentCutoff --blastProgram  $params.blastProgram --database newdb.fasta --seqFile subset.fa  --blastParams $params.blastParams --doNotParse $params.doNotParse --printSimSeqsFile $params.printSimSeqsFile --saveAllBlastFiles $params.saveAllBlastFiles --saveGoodBlastFiles $params.saveGoodBlastFiles --doNotExitOnBlastFailure $params.doNotExitOnBlastFailure --remMaskedRes $params.adjustMatchLength
   """
}

results = output_qch.collectFile(storeDir: params.outputDir, name: params.dataFile)
logResults = log_qch.collectFile(storeDir: params.outputDir, name: params.logFile)
zipResults = zip_qch.collectFile(storeDir: params.outputDir)

