#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



#include "SAIGE_test.hpp"
#include "SPA.hpp"

#include "UTIL.hpp"
#include "getMem.hpp"
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

namespace SAIGE {

SAIGEClass::SAIGEClass(
	arma::mat & t_XVX,
	arma::mat  t_XXVX_inv,
	arma::mat & t_XV,
	arma::mat & t_XVX_inv_XV,
	arma::mat & t_X,
	arma::vec &  t_S_a,
	arma::vec & t_res,
	arma::vec & t_mu2,
	arma::vec & t_mu,
	arma::vec & t_varRatio,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
	double t_SPA_Cutoff,
	arma::vec & t_tauvec,
	std::string t_traitType,
	arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	arma::umat & t_locationMat,
	arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
        std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
        double t_pCutoffforFirth,
        arma::vec & t_offset){

    m_XVX = t_XVX;
    m_XV = t_XV;
    m_XXVX_inv = t_XXVX_inv;
    m_XVX_inv_XV = t_XVX_inv_XV;
    m_X = t_X;
    m_S_a = t_S_a;
    m_res = t_res;
    m_mu2 = t_mu2;
    m_mu = t_mu;
    m_varRatio = t_varRatio;
    m_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
    m_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
    m_tauvec = t_tauvec;
    m_traitType = t_traitType;
    m_y = t_y;

    m_case_indices = arma::find(m_y == 1);
    m_ctrl_indices = arma::find(m_y == 0);

    m_n = t_y.size();
    m_p = t_XV.n_rows;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_impute_method =  t_impute_method;
    m_isCondition = t_isCondition;
    m_condition_genoIndex = t_condition_genoIndex;
    if(m_isCondition){	
	        m_numMarker_cond = t_condition_genoIndex.size();      
    }else{
		m_numMarker_cond = 0;
    }	    
    //m_p = t_X.nrow();

    if(m_traitType == "binary"){
      //if(m_isOutputAFinCaseCtrl){
        m_case_indices = arma::find(m_y == 1);
        m_ctrl_indices = arma::find(m_y == 0);
	m_n_case = m_case_indices.n_elem;
	m_n_ctrl = m_ctrl_indices.n_elem;
	m_is_Firth_beta = t_is_Firth_beta;
	m_pCutoffforFirth = t_pCutoffforFirth;
	m_offset = t_offset;
      //}
    }
    m_dimNum = t_dimNum;
    m_flagSparseGRM = t_flagSparseGRM;
    if(m_dimNum != 0){
	m_locationMat = t_locationMat;
    	m_valueVec = t_valueVec;
    }

}    


// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void SAIGEClass::set_seed(unsigned int seed){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

void SAIGEClass::scoreTest(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
		     bool m_flagSparseGRM,
		     arma::vec & t_P2Vec,
		     double& t_gy, 
		     bool t_is_region){
    arma::vec Sm, var2m;
    double S, var2;
    getadjG(t_GVec, t_gtilde);
    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y);
    }


    S = dot(t_gtilde, m_res);
    S = S/m_tauvec[0];


    if(!m_flagSparseGRM){
      t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0];
    }else{
      arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat();
      t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
    }	      
    var2m = dot(t_P2Vec , t_gtilde);
    var2 = var2m(0,0);
    double var1 = var2 * m_varRatioVal;
    double stat = S*S/var1;
    double t_pval;


    if (var1 < std::pow(std::numeric_limits<double>::min(), 2)){
        t_pval = 1;
    }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }

    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;

    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(stat);
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    arma::vec Z = A1.t() * g1;
    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    if(m_traitType == "binary"){
      mu21  = m_mu2.elem(t_indexForNonZero);
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }
    var1 = var2 * m_varRatioVal;
    S1 = dot(res1, g1_tilde);
    arma::mat res1X1_temp = (res1.t()) * X1;
    arma::vec res1X1 = res1X1_temp.t();
    S_a2 = m_S_a - res1X1;
    S2 = - arma::dot(S_a2,  Z);
    S = S1 + S2;
    S = S/m_tauvec[0];


    double stat = S*S/var1;
    double t_pval;
    if (var1 < std::pow(std::numeric_limits<double>::min(), 2)){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }
    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(stat);
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;

}


void SAIGEClass::getadjG(arma::vec & t_GVec, arma::vec & g){
/*	  double mem1, mem2;
  process_mem_usage(mem1, mem2);
   std::cout << "VM scoreTest a 2 in getadjG: " << mem1/10000 << "; RSS scoreTest a 2 in getadjGFast: " << mem2/10000 << std::endl;
*/
   g = m_XV * t_GVec;
 // process_mem_usage(mem1, mem2);
 //  std::cout << "VM scoreTest a 2 1 in getadjG: " << mem1/10000 << "; RSS scoreTest a 2 in getadjGFast: " << mem2/10000 << std::endl;
    g = t_GVec - m_XXVX_inv * g;
 // process_mem_usage(mem1, mem2);
 //  std::cout << "VM scoreTest a 2 2 in getadjG: " << mem1/10000 << "; RSS scoreTest a 2 in getadjGFast: " << mem2/10000 << std::endl;
}


void SAIGEClass::getadjGFast(arma::vec & t_GVec, arma::vec & g)
{

  // To increase computational efficiency when lots of GVec elements are 0

 arma::vec m_XVG(m_p, arma::fill::zeros);
  for(int i = 0; i < m_n; i++){
    if(t_GVec(i) != 0){
      m_XVG += m_XV.col(i) * t_GVec(i);
    }
  }

   //g*=m_XXVX_inv
  //g = t_GVec - m_XXVX_inv * m_XVG;
  //arma::vec g1(m_n, arma::fill::zeros);
   for(int i = 0; i < m_n; i++){
  	g(i) = -dot(m_XXVX_inv.row(i), m_XVG);
  }

  g += t_GVec ;
  
}


void SAIGEClass::get_mu(arma::vec & t_mu){
    t_mu = m_mu;
}

void SAIGEClass::getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices){
     t_case_indices = m_case_indices;
     t_ctrl_indices = m_ctrl_indices;
  }


void SAIGEClass::setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    m_locationMat = locationMatinR;
    m_valueVec = valueVecinR;
    m_dimNum = r;
}



arma::sp_mat SAIGEClass::gen_sp_SigmaMat() {
    arma::sp_mat resultMat(m_locationMat, m_valueVec, m_dimNum, m_dimNum);
    return resultMat;
}

// revised
// need to add sparse Sigma version 
// This function only uses variance ratio and does not use sparse GRM
void SAIGEClass::getMarkerPval(arma::vec & t_GVec,
			       arma::uvec & iIndex,
			       arma::uvec & iIndexComVec,
                               double& t_Beta, 
                               double& t_seBeta, 
                               double& t_pval,
			       double& t_pval_noSPA, 
                               double t_altFreq,
                               double& t_Tstat,
			       double& t_gy,
			       double& t_var1, 
			       bool & t_isSPAConverge, 
			       arma::vec & t_gtilde, 
			       bool & is_gtilde,
			       bool  is_region,
                               arma::vec & t_P2Vec,
			       bool t_isCondition,
			       double& t_Beta_c,
                           	double& t_seBeta_c,
                           	double& t_pval_c,
                          	double& t_pval_noSPA_c,
                           	double& t_Tstat_c,
                           	double& t_varT_c,
			   	arma::rowvec & t_G1tilde_P_G2tilde)
{
  //arma::vec adjGVec = getadjGFast(t_GVec);
  std::string t_pval_str;
  double t_var2, t_SPApval;
  //iIndex = arma::find(t_GVec != 0);
  //arma::vec t_gtilde;
  bool isScoreFast = true;


  if((t_altFreq > 0.05 && t_altFreq < 0.95) || m_flagSparseGRM || is_region){
    isScoreFast = false;
  }  
 //arma::vec timeoutput3 = getTime();
  if(!isScoreFast){
  	is_gtilde = true;
  	scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, m_flagSparseGRM, t_P2Vec, t_gy, is_region);
  }else{
  	is_gtilde = false;
	//arma::uvec iIndexVec = arma::find(t_GVec > 0);
        scoreTestFast(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2);
  }

  double StdStat = std::abs(t_Tstat) / sqrt(t_var1);
  t_isSPAConverge = false;

  double pval_noadj;
  try {
        pval_noadj = std::stod(t_pval_str);
  } catch (const std::invalid_argument&) {
        pval_noadj = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
  } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
        pval_noadj = 0;
  }


 //arma::vec timeoutput3_a = getTime();
  double q, qinv, m1, NAmu, NAsigma, tol1, p_iIndexComVecSize;

  //arma::uvec iIndexComVec = arma::find(t_GVec == 0);
  //arma::uvec iIndexVec = arma::find(t_GVec != 0);
 
  
  unsigned int iIndexComVecSize = iIndexComVec.n_elem;
  unsigned int iIndexSize = iIndex.n_elem;
 
  arma::vec gNB(iIndexSize, arma::fill::none);
  arma::vec gNA(iIndexComVecSize, arma::fill::none);
  arma::vec muNB(iIndexSize, arma::fill::none);
  arma::vec muNA(iIndexComVecSize, arma::fill::none);
/*
    std::cout << "iIndexComVecSize " << iIndexComVecSize << std::endl;
    std::cout << "iIndexSize " << iIndexSize << std::endl;
    std::cout << "gNA.n_elem 1 " << gNA.n_elem << std::endl;
        std::cout << "gNB.n_elem 1 " << gNB.n_elem << std::endl;
        std::cout << "muNA.n_elem 1 " << muNA.n_elem << std::endl;
        std::cout << "muNB.n_elem 1 " << muNB.n_elem << std::endl;
*/


  double gmuNB;
  if(StdStat > m_SPA_Cutoff && m_traitType != "quantitative"){

       if(!is_gtilde){
          t_gtilde.resize(m_n);
          getadjGFast(t_GVec, t_gtilde);
	  is_gtilde = true;
       }
	//int t_gtilden = t_gtilde.n_elem;
        p_iIndexComVecSize = double(iIndexComVecSize)/m_n;
   	m1 = dot(m_mu, t_gtilde);

	if(p_iIndexComVecSize >= 0.5){
		unsigned int j1 = 0;
		unsigned int j2 = 0;
/*
		for(unsigned int j = 0; j < m_n ; j++){	
			//std::cout << "j " << j << std::endl;
			if(t_GVec(j) != 0){
			//std::cout << "j1 " << j1 << std::endl;
				gNB(j1) = t_gtilde(j);
				muNB(j1) = m_mu(j);
				j1 = j1 + 1;	
			}else{
			//std::cout << "j2 " << j2 << std::endl;
				gNA(j2) = t_gtilde(j);
				muNA(j2) = m_mu(j);
				j2 = j2 + 1;
          //process_mem_usage(mem1, mem2);
//   std::cout << "VM 4 a 1.3c: " << mem1/1000 << "; RSS 4 a 1.3: " << mem2/1000 << std::endl;
			}	
		}
		*/
//	std::cout << "gNB.n_elem " <<  gNB.n_elem << std::endl;	
//	std::cout << "gNA.n_elem " <<  gNA.n_elem << std::endl;	
	gNB = t_gtilde(iIndex);
	gNA = t_gtilde(iIndexComVec);
   	muNB = m_mu(iIndex);
   	muNA = m_mu(iIndexComVec);

	/*
	    std::cout << "gNA.n_elem 2 " << gNA.n_elem << std::endl;
        std::cout << "gNB.n_elem 2 " << gNB.n_elem << std::endl;
        std::cout << "muNA.n_elem 2 " << muNA.n_elem << std::endl;
        std::cout << "muNB.n_elem 2 " << muNB.n_elem << std::endl;
*/


  	gmuNB = dot(gNB,muNB);	 
   	NAmu= m1-gmuNB;
	}
	/*else{
		gNA.clear();
		gNB.clear();
		muNA.clear();
		gNB.clear();

	}*/

   	if(m_traitType == "binary"){
                q = t_Tstat/sqrt(t_var1/t_var2) + m1;

                if((q-m1) > 0){
                        qinv = -1 * std::abs(q-m1) + m1;
                }else if ((q-m1) == 0){
                        qinv =  m1;
                }else{
                        qinv = std::abs(q-m1) + m1;
                }
		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % (1-muNB) % arma::pow(gNB,2));
		}
        }else if(m_traitType == "survival"){
                q = t_Tstat/sqrt(t_var1/t_var2);
                qinv = -q;
		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % arma::pow(gNB,2));
		}
        }
    	bool logp=false;
	double tol0 = std::numeric_limits<double>::epsilon();
	tol1 = std::pow(tol0, 0.25);
	if(p_iIndexComVecSize >= 0.5){
		//std::cout << "SPA_fast" << std::endl;
        	SPA_fast(m_mu, t_gtilde, q, qinv, pval_noadj, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);
	}else{
		//std::cout << "SPA" << std::endl;
		SPA(m_mu, t_gtilde, q, qinv, pval_noadj, tol1, logp, m_traitType, t_SPApval, t_isSPAConverge);	
	}


    	boost::math::normal ns;
	t_pval = t_SPApval;
    	double t_qval;
        try {
           t_qval = boost::math::quantile(ns, t_pval/2);
           t_qval = fabs(t_qval);
           t_seBeta = fabs(t_Beta)/t_qval;
        }catch (const std::overflow_error&) {
          t_qval = std::numeric_limits<double>::infinity();
          t_seBeta = 0;
        } 
  }
   t_pval_noSPA = pval_noadj; 
   if(m_traitType!="quantitative"){
        if(t_isSPAConverge){
                t_pval = t_SPApval;
        }else{
                t_pval = pval_noadj;
        }
   }else{
        t_pval = t_pval_noSPA;
   }

   if(m_traitType!="quantitative" & m_is_Firth_beta & t_pval <= m_pCutoffforFirth){
	if(!is_gtilde){
                getadjG(t_GVec, t_gtilde);
                is_gtilde = true;
        }
	arma::mat x(t_GVec.n_elem, 2, arma::fill::zeros);	
	x.col(1) = t_gtilde;
	arma::vec init(2, arma::fill::zeros);
	fast_logistf_fit_simple(x, m_y, m_offset, true, init, 50, 15, 15, 1e-5, 1e-5, 1e-5, t_Beta ,t_seBeta);	
   }
   
 //arma::vec timeoutput4 = getTime();
 //printTime(timeoutput3, timeoutput3_a, "Test Marker  ScoreTest");
//printTime(timeoutput3, timeoutput4, "Test Marker 3 to 4");
//printTime(timeoutput3_a, timeoutput4, "Test Marker SPA");

   //condition
   if(t_isCondition){
	if(!is_gtilde){
        	getadjG(t_GVec, t_gtilde);
        	is_gtilde = true;
        }
        t_G1tilde_P_G2tilde = sqrt(m_varRatioVal) * t_gtilde.t() * m_P2Mat_cond;
        arma::vec t_Tstat_ctemp =  t_G1tilde_P_G2tilde * m_VarInvMat_cond * m_Tstat_cond;
	arma::mat tempgP2 = t_gtilde.t() * m_P2Mat_cond;

    	t_Tstat_c = t_Tstat - t_Tstat_ctemp(0);
    	arma::vec t_varT_ctemp = t_G1tilde_P_G2tilde * m_VarInvMat_cond * (t_G1tilde_P_G2tilde.t());
    	t_varT_c = t_var1 - t_varT_ctemp(0);

    double S_c = t_Tstat_c;

    double stat_c = S_c*S_c/t_varT_c;
     if (t_varT_c < std::pow(std::numeric_limits<double>::min(), 2)){
        t_pval_noSPA_c = 1;
	stat_c = 0;
     }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval_noSPA_c = boost::math::cdf(complement(chisq_dist, stat_c));
     }

    char pValueBuf_c[100];
    if (t_pval_noSPA_c != 0)
        sprintf(pValueBuf_c, "%.6E", t_pval_noSPA_c);
    else {
        double log10p_c = log10(2.0) - M_LOG10E*stat_c/2 - 0.5*log10(stat_c*2*M_PI);
        int exponent_c = floor(log10p_c);
        double fraction_c = pow(10.0, log10p_c - exponent_c);
        if (fraction_c >= 9.95) {
          fraction_c = 1;
           exponent_c++;
         }
        sprintf(pValueBuf_c, "%.1fE%d", fraction_c, exponent_c);
    }
    std::string buffAsStdStr_c = pValueBuf_c;
    std::string& t_pval_noSPA_str_c = buffAsStdStr_c;

    t_Beta_c = S_c/t_varT_c;
    t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
    t_Tstat_c = S_c;
/*
*/

  double pval_noSPA_c;  
  try {
        pval_noSPA_c = std::stod(t_pval_noSPA_str_c);
  } catch (const std::invalid_argument&) {
        pval_noSPA_c = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
  } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
        pval_noSPA_c = 0;
  }
  t_pval_noSPA_c = pval_noSPA_c; 

  //std::cout << "stat_c " << stat_c << std::endl;
  //std::cout << "t_varT_c " << t_varT_c << std::endl;

    if(m_traitType != "quantitative" && stat_c > std::pow(m_SPA_Cutoff,2)){
	bool t_isSPAConverge_c;
	double q_c, qinv_c, pval_noadj_c, SPApval_c;    
	if(m_traitType == "binary"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2) + m1;

                if((q_c-m1) > 0){
                        qinv_c = -1 * std::abs(q_c-m1) + m1;
                }else if ((q_c-m1) == 0){
                        qinv_c =  m1;
                }else{
                        qinv_c = std::abs(q_c-m1) + m1;
                }
        }else if(m_traitType == "survival"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2);
                qinv = -q_c;
        }

        bool logp=false;

        if(p_iIndexComVecSize >= 0.5){
                SPA_fast(m_mu, t_gtilde, q_c, qinv_c, pval_noadj_c, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, SPApval_c, t_isSPAConverge_c);
        }else{
                SPA(m_mu, t_gtilde, q_c, qinv_c, pval_noadj_c, tol1, logp, m_traitType, SPApval_c, t_isSPAConverge_c);
        }

        boost::math::normal ns;
        t_pval_c = SPApval_c;
        double t_qval_c;
        try {
           t_qval_c = boost::math::quantile(ns, t_pval_c/2);
           t_qval_c = fabs(t_qval_c);
           t_seBeta_c = fabs(t_Beta_c)/t_qval_c;
        }catch (const std::overflow_error&) {
          t_qval_c = std::numeric_limits<double>::infinity();
          t_seBeta_c = 0;
        }
    }else{
    	t_pval_c = t_pval_noSPA_c;	    
    }	    
   }

    gNA.clear();
    gNB.clear();
    muNA.clear();
    gNB.clear();
}


void SAIGEClass::assignVarianceRatio(double MAC){
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){    	    
		m_varRatioVal = m_varRatio(i);
	}	
    }    
}

void SAIGEClass::assignSingleVarianceRatio(){ 
	m_varRatioVal = m_varRatio(0);
}


void SAIGEClass::assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      arma::vec & t_p_cond
      ){
	m_P2Mat_cond = t_P2Mat_cond;
	m_VarInvMat_cond = t_VarInvMat_cond;
	m_VarMat_cond = t_VarMat_cond;
	m_Tstat_cond = t_Tstat_cond;
	m_MAF_cond = t_MAF_cond;
	m_qsum_cond = t_qsum_cond;
	m_gsum_cond = t_gsum_cond;
	m_G2_Weight_cond = t_G2_Weight_cond;
	m_p_cond = t_p_cond;
}

void SAIGEClass::assignConditionFactors_scalefactor(
	arma::vec & t_scalefactor_G2_cond	
		){
	m_scalefactor_G2_cond = t_scalefactor_G2_cond;
	arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
	arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t(); 
	arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * m_VarMat_cond * scalefactor_G2_cond_Mat;
	arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
       
	m_VarInvMat_cond_scaled_weighted = VarMat_cond_scaled_weighted.i();
	//m_VarInvMat_cond_region_binary = (1/scalefactor_G2_cond_Mat) * m_VarInvMat_cond	* (1/scalefactor_G2_cond_Mat);
	
}

void SAIGEClass::extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv){
	t_XV = m_XV;
	t_XXVX_inv = m_XXVX_inv;	
}



void SAIGEClass::fast_logistf_fit_simple(arma::mat & x,
                arma::vec & y,
                arma::vec & offset,
                bool firth,
        arma::vec init,
        int maxit,
        int maxstep,
        int maxhs,
        double lconv,
        double gconv,
        double xconv,
        double & beta_G,
        double & sebeta_G){
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  arma::vec pi_0 = -x * beta - offset;
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  while(iter <= maxit){
        beta_old = beta;
        arma::vec wpi = pi % (1 - pi);
        arma::vec W2 = arma::sqrt(wpi);
        //arma::vec wpi_sqrt = arma::sqrt(wpi);
        //arma::vec W2 = weight % wpi_sqrt;
        arma::mat XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                XW2.col(j) = x.col(j) % W2;
        }

        arma::mat Q;
        arma::mat R;
        arma::qr_econ(Q, R, XW2);
        arma::vec h = Q % Q * oneVec;
        arma::vec U_star(2, arma::fill::zeros);
        arma::vec ypih;
        if(firth){
                ypih = (y - pi) + (h % (0.5 - pi));
        }else{
                ypih = (y - pi);
        }
        //ypih.print();
        arma::vec xcol(n, arma::fill::zeros);
        U_star = x.t() * ypih;

        arma::mat XX_XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                xcol = x.col(j);
                XX_XW2.col(j) = xcol % W2;
        }
        arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
        bool isinv = arma::inv_sympd (XX_covs, XX_Fisher);
        if(!isinv){
                break;
        }
        //}
        arma::vec delta = XX_covs * U_star;
        delta.replace(arma::datum::nan, 0);

        double mx = arma::max(arma::abs(delta))/maxstep;
        if(mx > 1){
                delta = delta/mx;
        }
        evals = evals + 1;
        iter = iter + 1;
        beta = beta + delta;
        pi_0 = -x * beta - offset;
        pi_0 = arma::exp(pi_0) + 1;
        pi = 1/pi_0;
        if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
                break;
        }
  }
        arma::mat var;
        if(XX_covs.has_nan()){
                var = XX_covs;
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = sqrt(XX_covs(1,1));
        }
        //return beta;
}









}
