#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------

if(!params.fastaSubsetSize) {
  throw new Exception("Missing params.fastaSubsetSize")
}

if(params.seqFile) {
  seqs = Channel.fromPath( params.seqFile )
           .splitFasta( by:params.fastaSubsetSize, file:true  )
}
else {
  throw new Exception("Missing params.seqFile")
}

if (params.preConfiguredDatabase) {
  if (!params.databaseDir) {
    throw new Exception("Missing params.databaseDir")
  }
  else if(!params.databaseBaseName) {
    throw new Exception("Missing params.databaseBaseName")
  }
}


//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { nonConfiguredDatabase } from './modules/blastSimilarity.nf'
include { preConfiguredDatabase } from './modules/blastSimilarity.nf'


//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  
  if (!params.preConfiguredDatabase) {
    nonConfiguredDatabase(seqs)
  }
   
  else if (params.preConfiguredDatabase) {
    proConfiguredDatabase(seqs)
  }

}

