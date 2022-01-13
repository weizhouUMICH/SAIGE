void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint32_t> & t_genoIndex,
                           std::string t_traitType,
                           unsigned int t_n)           // sample size
{
  unsigned int q = t_genoIndex.size();
  arma::mat P1Mat(q, t_n, fill:none);
  arma::mat P2Mat(t_n, q, fill:none);
  std::vector<double> TstatVec(q, arma::datum::nan);

  double Beta, seBeta, pval, pval_noSPA, Tstat, varT;
  bool isSPAConverge, is_gtilde;
  arma::vec P1Vec(t_n), P2Vec(t_n);

    for(unsigned int i = 0; i < q; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;

    uint32_t gIndex = t_genoIndex.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);

    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    //std::cout << "i " << i << std::endl;
    //std::cout << "info " << info << std::endl;

    flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, g_impute_method);


   Unified_getMarkerPval(
                    GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexForNonZero, Beta, seBeta, pval, pval_noSPA, Tstat, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec);
      TstatVec.at(i) = Tstat;

      P1Mat.row(i) = gtildeVec.t();
      P2Mat.col(i) = P2Vec;
      //std::cout << "here5" << std::endl;
   ptr_gSAIGEobj->assignConditionFactors(P1Mat,
		   			P2Mat,
					TstatVec)
}
      
