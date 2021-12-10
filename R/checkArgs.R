checkArgBool = function(arg0, arg0name){
  if(!arg0 %in% c(TRUE, FALSE))
    stop(arg0name, " should be TRUE or FALSE.")
}	

checkArgNumeric = function(arg0, arg0name, minVal = NULL, maxVal = NULL){
 if(!is.numeric(arg0)){
     stop(arg0name, " should be a numeric value.")
  }	 

  if(!is.null(minVal)){
    if(arg0 < minVal){
      stop(arg0name, " should be a numeric value greater or equal to ", minVal,"\n")
    }	    
  }

  if(!is.null(maxVal)){
    if(arg0 > maxVal){
      stop(arg0name, " should be a numeric value less or equal to ", maxVal,"\n")
    }
  }
}



checkArgsListNumeric = function(start,
end,
maxMissingRate,
minMAC,
minMAF,
maxMAFforGroupTest,
minInfo,
SPAcutoff,
numLinesOutput,
MACCutoff_to_CollapseUltraRare,
DosageCutoff_for_UltraRarePresence,
weightMAFcutoff,
dosageZerodCutoff)
		     
{
	checkArgNumeric(start, deparse(substitute(start)), 1, 250000000)
	checkArgNumeric(end, deparse(substitute(end)), 1, 250000000)
	checkArgNumeric(maxMissingRate, deparse(substitute(maxMissingRate)), 0, 0.5)
	checkArgNumeric(minMAC, deparse(substitute(minMAC)), minVal = 20)
	checkArgNumeric(minMAF, deparse(substitute(minMAF)), 0, 0.5)
	checkArgNumeric(maxMAFforGroupTest, deparse(substitute(maxMAFforGroupTest)), 0, 0.5)
	checkArgNumeric(minInfo, deparse(substitute(start)), 0, 1)
	checkArgNumeric(SPAcutoff, deparse(substitute(SPAcutoff)), 0.5, 4)
	checkArgNumeric(numLinesOutput, deparse(substitute(numLinesOutput)), 1, 100000)
	checkArgNumeric(MACCutoff_to_CollapseUltraRare, deparse(substitute(MACCutoff_to_CollapseUltraRare)), 1, 100)
	checkArgNumeric(DosageCutoff_for_UltraRarePresence, deparse(substitute(DosageCutoff_for_UltraRarePresence)), 0, 2)
	checkArgNumeric(weightMAFcutoff, deparse(substitute(weightMAFcutoff)), 0, 0.5)
	
	checkArgNumeric(dosageZerodCutoff, deparse(substitute(dosageZerodCutoff)), 0, 0.5)
	cat("Any dosages <= ", dosageZerodCutoff, " for genetic variants with MAC <= 10 are set to be 0 in group tests\n")
}


checkArgsListBool = function( IsDropMissingDosages,
 IsSparse,
 IsOutputAFinCaseCtrl,
 IsOutputHetHomCountsinCaseCtrl,
 IsOutputNinCaseCtrl,
 IsOutputlogPforSingle,
 LOCO,
 IsSingleVarinGroupTest,
 IsOutputPvalueNAinGroupTestforBinary,
 IsAccountforCasecontrolImbalanceinGroupTest,
 IsOutputBETASEinBurdenTest,
 IsOutputMAFinCaseCtrlinGroupTest,
 is_rewrite_XnonPAR_forMales){

checkArgBool(IsDropMissingDosages, deparse(substitute(IsDropMissingDosages)))
checkArgBool(IsSparse, deparse(substitute(IsSparse)))
checkArgBool(IsOutputAFinCaseCtrl, deparse(substitute(IsOutputAFinCaseCtrl)))
checkArgBool(IsOutputHetHomCountsinCaseCtrl, deparse(substitute(IsOutputHetHomCountsinCaseCtrl)))
checkArgBool(IsOutputNinCaseCtrl, deparse(substitute(IsOutputNinCaseCtrl)))
checkArgBool(IsOutputlogPforSingle, deparse(substitute(IsOutputlogPforSingle)))
if (IsOutputlogPforSingle) {
  at("IsOutputlogPforSingle = TRUE. NOTE: log(Pvalue) will be output ONLY for single-variant assoc tests\n")
}

checkArgBool(LOCO, deparse(substitute(LOCO)))
checkArgBool(IsSingleVarinGroupTest, deparse(substitute(IsSingleVarinGroupTest)))
checkArgBool(IsOutputPvalueNAinGroupTestforBinary, deparse(substitute(IsOutputPvalueNAinGroupTestforBinary)))
checkArgBool(IsAccountforCasecontrolImbalanceinGroupTest, deparse(substitute(IsAccountforCasecontrolImbalanceinGroupTest)))
checkArgBool(IsOutputBETASEinBurdenTest, deparse(substitute(IsOutputBETASEinBurdenTest)))
checkArgBool(IsOutputMAFinCaseCtrlinGroupTest, deparse(substitute(IsOutputMAFinCaseCtrlinGroupTest)))
checkArgBool(is_rewrite_XnonPAR_forMales, deparse(substitute(is_rewrite_XnonPAR_forMales)))

}	
