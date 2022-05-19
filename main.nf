nextflow.enable.dsl=1
seq_qch = Channel.fromPath(params.seqFile).splitFasta( by:1, file:true  )

process blastSimilarity {
   input:
   path 'subset.fa' from seq_qch
   output:
   path 'blastSimilarity.out' into output_qch
   path 'blastSimilarity.log' into log_qch
   """
   cat $params.databaseFasta > newdb.fasta
   makeblastdb -in newdb.fasta -dbtype prot
   blastSimilarity --regex $params.regex --pValCutoff  $params.pValCutoff --lengthCutoff $params.lengthCutoff --percentCutoff  $params.percentCutoff --blastBinDir /usr/bin/ncbi-blast-2.13.0+/bin --blastProgram  $params.blastProgram --database newdb.fasta --seqFile subset.fa  --blastParams $params.blastParams --blastVendor $params.blastVendor
   """
}

results = output_qch.collectFile(storeDir: params.outputDir, name: "blastSimilarity.out")
results = log_qch.collectFile(storeDir: params.outputDir, name: "blastSimilarity.log")


