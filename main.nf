nextflow.enable.dsl=1
seq_qch = Channel.fromPath(params.seqFile).splitFasta( by:1, file:true  )
db_vch = Channel.value()
process createDatabase {
    output:
    path 'newdb.fasta.*' into db_vch
    """
    createDatabase.pl --dbFile $params.databaseFasta --databaseType $params.databaseType --blastVendor $params.blastVendor
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
   path 'blastAnal.log' optional true into out_qch
   """
   blastSimilarity --regex $params.regex --pValCutoff  $params.pValCutoff --lengthCutoff $params.lengthCutoff --percentCutoff  $params.percentCutoff --blastBinDir /usr/bin/ncbi-blast-2.13.0+/bin --blastProgram  $params.blastProgram --database newdb.fasta --seqFile subset.fa  --blastParams $params.blastParams --blastVendor $params.blastVendor --doNotParse $params.doNotParse --printSimSeqsFile $params.printSimSeqsFile --saveAllBlastFiles $params.saveAllBlastFiles --saveGoodBlastFiles $params.saveGoodBlastFiles --doNotExitOnBlastFailure $params.doNotExitOnBlastFailure --databaseType $params.databaseType
   """
}

results = output_qch.collectFile(storeDir: params.outputDir, name: params.dataFile)
logResults = log_qch.collectFile(storeDir: params.outputDir, name: params.logFile)
zipResults = zip_qch.collectFile(storeDir: params.outputDir)
outResults = out_qch.collectFile(storeDir: params.outputDir)