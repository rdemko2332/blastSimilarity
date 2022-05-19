nextflow.enable.dsl=1
seq_qch = Channel.fromPath(params.inputFilePath).splitFasta( by:1, file:true  )
db_vch = Channel.fromPath(params.databaseFasta)
process createDatabase {
    input:
    path 'newdb.fasta' from db_vch
    
    """
    cd /work
    makeblastdb -in newdb.fasta -dbtype prot
    """
}

process blastSimilarity {
   input:
   path 'subset.fa' from seq_qch
   output:
   path 'blastSimilarity.out' into output_qch
   """
   blastSimilarity --regex '(S+)' --pValCutoff 1e-5 --lengthCutoff 10 --percentCutoff 20 --blastBinDir /usr/bin/ncbi-blast-2.13.0+/bin --blastProgram blastp --database /work/newdb.fasta --seqFile subset.fa  --blastParams $params.blastParams --blastVendor ncbi
   """
}

results = output_qch.collectFile(storeDir: params.outputDir, name: "blastSimilarity.out")



