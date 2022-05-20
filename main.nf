nextflow.enable.dsl=1
seq_qch = Channel.fromPath(params.seqFile).splitFasta( by:1, file:true  )

process blastSimilarity {
   //publishDir path: params.outputDir, pattern: "*.gz*"
   input:
   path 'subset.fa' from seq_qch
   output:
   path 'blastSimilarity.out' into output_qch
   path 'blastSimilarity.log' into log_qch
   path '*.gz*' into test_qch
   """
   cat $params.databaseFasta > newdb.fasta
   makeblastdb -in newdb.fasta -dbtype prot
   blastSimilarity --regex $params.regex --pValCutoff  $params.pValCutoff --lengthCutoff $params.lengthCutoff --percentCutoff  $params.percentCutoff --blastBinDir /usr/bin/ncbi-blast-2.13.0+/bin --blastProgram  $params.blastProgram --database newdb.fasta --seqFile subset.fa  --blastParams $params.blastParams --blastVendor $params.blastVendor --doNotParse $params.doNotParse --printSimSeqsFile $params.printSimSeqsFile --saveAllBlastFiles $params.saveAllBlastFiles --saveGoodBlastFiles $params.saveGoodBlastFiles --doNotExitOnBlastFailure $params.doNotExitOnBlastFailure --blastFileDir $params.blastFileDir --fileExtension $params.fileExtension
   """
}

results = output_qch.collectFile(storeDir: params.outputDir, name: params.dataFile)
logResults = log_qch.collectFile(storeDir: params.outputDir, name: params.logFile)
testResults = test_qch.collectFile(storeDir: params.outputDir)