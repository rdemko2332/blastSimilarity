#!/usr/bin/env bash

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
