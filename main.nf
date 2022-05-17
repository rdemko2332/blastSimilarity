nextflow.enable.dsl=1
seq_qch = Channel.fromPath(params.inputFilePath).splitFasta( by:1, file:true  )

process blast {
    
    input:
    file 'subset.fa' from seq_qch
    output:
    file 'blast.stderr' into blast_qch
        
    """
    rpsblastn -d $params.database -i subset.fa $params.blastParams 2> blast.stderr 
    """
}

results = blast_qch.collectFile(storeDir: params.outputDir, name: params.outputFileName)
errors = error_qch.collectFile(storeDir: params.outputDir)


