# Nextflow Conversion of BlastSimilarityTask

**Explanation of nextflow.config file parameters:**
- **blastProgram** (*String*) Name of NCBI blast tool you want to run.
- **inputFile** (*String*) The location of your inputFile.
- **preConfiguredDatabase** (*Boolean*) If you have databasefiles generated from NCBI's makeblastdb, there is no need to generate these files. If this is set to true, you will need to supply databaseDir and databaseBaseName.
- **databaseDir** (*String*) The path to the directory containing the database files. There can be other files in this directory, but any file beginning with the databaseBaseName will be brought into the process.
- **databaseFasta** (*String*) The location of the fasta file that you would like to use to create your database. Needed if preConfiguredDatabase is false.
- **dataFile** (*String*) How you would like the main output file to be named.
- **logFile** (*String*) How you would like the log file to be named.
- **outputDir** (*String*) Directory where you would like the output fiels to be stored.
- **saveAllBlastFiles** (*Boolean*) If true, the blast output for each time blast is ran. If you have 9 sequences in your input file, and you have fastaSubsetSize as 1, you will recieve 9 zipped files. If fastaSubsetSize is equal to three, you will recieve 3 zipped files. Zipped file names will be the sequence identifier for the first sequence in the group being run that is put into the file (also will be the last in the zip file.
- **saveGoodBlastFiles (*Boolean*) Similar to saveAllBlastFiles, expect only files that contain a hit will be saved. saveGood and saveAll should not both be true.
- **doNotParse** (*Boolean*) This tool operates in two steps, running blast and grepping through the output to collect and return values. If doNotParse is true, only the blast output is generated and returned. If false, then the output will continue on to the processing step.
- **printSimSeqsFile** (*Boolean*) Changes the output format of dataFile. Returns sequence accession from seqFile, the taxon it matched with from the database, the p-value, the exponent for the p-Value, and some stats per identity and per match. 
- **blastParamsFile** (*String*) The file location of the file containing additional blast paramenters. These can just be written out in the file as if you were using them on the command line.
- **fastaSubsetSize** (*Int*) Number of sequences per split of seqFile passed to blastSimilarity process.
 
### Get Started
  * Install Nextflow
    
    `curl https://get.nextflow.io | bash`
  
  * Run the script
    
    `nextflow run VEuPathDB/blastSimilarity -with-trace -c  <config_file> -r main`
