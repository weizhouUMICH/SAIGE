Rcpp::List assign_conditionMarkers_factors(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint32_t> & t_genoIndex,
                           std::string t_traitType,
                           unsigned int t_n)           // sample size
{
  unsigned int q = t_genoIndex.size();
  arma::mat P1Mat(q, t_n, fill:none);
  arma::mat P2Mat(t_n, q, fill:none);
  std::vector<std::string> markerVec(q);
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MACVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MAFVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> N_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> N_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> N_Vec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q, arma::datum::nan);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q, arma::datum::nan);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q, arma::datum::nan);
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);
  std::vector<double> pvalVec(q, arma::datum::nan);
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  std::vector<double> pvalNAVec(q, arma::datum::nan);

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

    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
     markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    altCountsVec.at(i) = altCounts;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    imputationInfoVec.at(i) = imputeInfo;

   Unified_getMarkerPval(
                    GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexForNonZero, Beta, seBeta, pval, pval_noSPA, Tstat, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec);
      BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
      seBetaVec.at(i) = seBeta;
      pvalVec.at(i) = pval;
      pvalNAVec.at(i) = pval_noSPA;
      TstatVec.at(i) = Tstat * (1 - 2*flip);
      varTVec.at(i) = varT;
      int n = GVec.size();

      P1Mat.row(i) = gtildeVec.t();
      P2Mat.col(i) = P2Vec;
      //std::cout << "here5" << std::endl;

      arma::vec dosage_case, dosage_ctrl;
      if(t_traitType == "binary"){
                        dosage_case = GVec.elem(g_case_indices);
                        dosage_ctrl = GVec.elem(g_ctrl_indices);
      }

      if(t_traitType == "binary"){
      AF_case = arma::mean(dosage_case) /2;
      AF_ctrl = arma::mean(dosage_ctrl) /2;
      if(flip){
         AF_case = 1-AF_case;
         AF_ctrl = 1-AF_ctrl;
      }
      AF_caseVec.at(i) = AF_case;
      AF_ctrlVec.at(i) = AF_ctrl;
      N_caseVec.at(i) = dosage_case.n_elem;
      N_ctrlVec.at(i) = dosage_ctrl.n_elem;
     }else if(t_traitType == "quantitative"){
      N_Vec.at(i) = n;
    }

   Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("P1Mat") = P1Mat,
		   			  Rcpp::Named("P2Mat") = P2Mat,
                                          Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("MAFVec") = MAFVec,
                                          Rcpp::Named("altCountsVec") = altCountsVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                         Rcpp::Named("imputationInfoVec") =imputationInfoVec,
                                          Rcpp::Named("pvalVec") = pvalVec,
                                          Rcpp::Named("BetaVec") = BetaVec,
                                          Rcpp::Named("seBetaVec") = seBetaVec,
                                          Rcpp::Named("TstatVec") = TstatVec,
                                          Rcpp::Named("varTVec") = varTVec);

  if(t_traitType == "binary"){
    OutList.push_back(pvalNAVec, "pvalNAVec");
    OutList.push_back(AF_caseVec, "AF_caseVec");
    OutList.push_back(AF_ctrlVec, "AF_ctrlVec");
    OutList.push_back(N_caseVec, "N_caseVec");
    OutList.push_back(N_ctrlVec, "N_ctrlVec");
  }else if(t_traitType == "quantitative"){
    OutList.push_back(N_Vec, "N_Vec");
  }
  return OutList;
}
      
