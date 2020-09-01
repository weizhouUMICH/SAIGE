#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "../thirdParty/bgen/genfile/include/genfile/bgen/bgen.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/View.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/IndexQuery.hpp"
#include <sstream>
#include <time.h>
#include <Rcpp.h>
#include <stdint.h>


// #define DEBUG 1

namespace {
        template< typename Integer >
        std::string atoi( Integer const value ) {
                std::ostringstream stream ;
                stream << value ;
                return stream.str() ;
        }

}


//A global variable for the dosage file to test

std::auto_ptr< std::istream > gm_stream;

// 0 for additive, 1 for recessive, 2 for dominant
int dosage_type;

uint32_t  gm_offset ;
genfile::bgen::Context gm_context ;
bool gm_have_sample_ids ;
//std::vector< std::string > gm_sample_ids ;

// Added by SLEE for parsing 09/07/2017
Rcpp::IntegerVector gm_sample_idx;
Rcpp::IntegerVector cc_idx;

int gmtest_samplesize;

//genoToTest_bgenDosage
genfile::bgen::View::UniquePtr genoToTest_bgenDosage;
bool isQuery;
bool isReadVariantBgen = true;
//bool isOutputHetHomCountsinCaseCtrl = false;
bool isDropMissingDosages_bgen = false;

//double bgenMinMAF = 0;
//double bgenMinINFO = 0;
double markerInfo;
int numSamples_bgen;
// bool isDropMissingDosages_bgen = false;


// // [[Rcpp::export]]
//void setIsOutputHetHomCountsinCaseCtrl( bool isoutputhethom ) {
//  isOutputHetHomCountsinCaseCtrl = isoutputhethom;
//  std::cout << "IsOutputHetHomCountsinCaseCtrl = true, so the Heterozygous and homozygous counts in cases and controls will be output. " << std::endl; 
//}

// [[Rcpp::export]]
void setIsDropMissingDosages_bgen (bool isdropmissingdosages){
  isDropMissingDosages_bgen = isdropmissingdosages;
}

// // [[Rcpp::export]]
//void setIsDropMissingDosages_bgen(bool isDropMissing){
//  isDropMissingDosages_bgen = isDropMissing;
//}


// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct ProbSetter {
        typedef std::vector< std::vector< double > > Data ;
        ProbSetter( Data* result ):
                m_result( result ),
                m_sample_i(0)
        {}

        // Called once allowing us to set storage.
        void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
                m_result->clear() ;
                m_result->resize( number_of_samples ) ;
        }

        // If present with this signature, called once after initialise()
        // to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
        // This enables us to set up storage for the data ahead of time.
        void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
                for( std::size_t i = 0; i < m_result->size(); ++i ) {
                        m_result->at( i ).reserve( max_entries ) ;
                }
        }

        // Called once per sample to determine whether we want data for this sample
        bool set_sample( std::size_t i ) {
                m_sample_i = i ;
                // Yes, here we want info for all samples.
                return true ;
        }

        // Called once per sample to set the number of probabilities that are present.
        void set_number_of_entries(
                std::size_t ploidy,
                std::size_t number_of_entries,
                genfile::OrderType order_type,
                genfile::ValueType value_type
        ) {
                assert( value_type == genfile::eProbability ) ;
                m_result->at( m_sample_i ).resize( number_of_entries ) ;
                m_entry_i = 0 ;
        }

        // Called once for each genotype (or haplotype) probability per sample.
        void set_value( uint32_t, double value ) {
                m_result->at( m_sample_i ).at( m_entry_i++ ) = value ;
        }

        // Ditto, but called if data is missing for this sample.
        void set_value( uint32_t, genfile::MissingValue value ) {
                // Here we encode missing probabilities with -1
                m_result->at( m_sample_i ).at( m_entry_i++ ) = -1 ;
        }

        // If present with this signature, called once after all data has been set.
        void finalise() {
                // nothing to do in this implementation.
        }

private:
        Data* m_result ;
        std::size_t m_sample_i ;
        std::size_t m_entry_i ;
} ;


// Read genotype probability data for the SNP just read using read_variant()
// After calling this method it should be safe to call read_variant() to fetch
// the next variant from the file.
void read_probs( std::vector< std::vector< double > >* probs ) {
        ProbSetter setter( probs ) ;
	genoToTest_bgenDosage->read_genotype_data_block(setter);
}




// [[Rcpp::export]]
int setgenoTest_bgenDosage(std::string & filename,
	std::string & index_filename,
	Rcpp::DataFrame & ranges_to_include,
	Rcpp::DataFrame & ranges_to_exclude,
	std::vector< std::string > const& ids_to_include,
	std::vector< std::string > const& ids_to_exclude,
	std::string & analysis_type
){
  
  if (analysis_type == "recessive") {
    dosage_type = 1;
    std::cout << "running recessive analysis" << std::endl;
  } else if (analysis_type == "dominant") {
    dosage_type = 2;
    std::cout << "running dominant analysis" << std::endl;
  } else {
    dosage_type = 0;
    std::cout << "running additive analysis" << std::endl;
  }

   //bgenMinMAF = bgenMinMaf;
   //bgenMinINFO = bgenMinInfo;
   int numMarkers;
   {
     using namespace genfile::bgen ;
     using namespace Rcpp;
     isQuery = false; 
     if(index_filename == ""){
       isQuery = false;
       std::cout << "no index file for bgen is provided" << std::endl; 
     }else if(ranges_to_include.nrow() == 0 &&
	  ranges_to_exclude.nrow() == 0 &&
	  ids_to_include.size() == 0 &&
	  ids_to_exclude.size() == 0){
       std::cout << "no query list is provided" << std::endl;
       isQuery = false;
     }else{
       isQuery = true;
     }


     if(isQuery){
	std::cout << "TEST " << ids_to_include.size() << std::endl;
        genoToTest_bgenDosage = View::create( filename ) ;
	std::cout << "TEST 1 OK" << std::endl;
	IndexQuery::UniquePtr query = IndexQuery::create( index_filename ) ;
	std::cout << "TEST2 " << ids_to_include.size() << std::endl;

	std::cout << "ranges_to_include.nrow() " << ranges_to_include.nrow() << std::endl;
	std::cout << "ranges_to_exclude.nrow() " << ranges_to_exclude.nrow() << std::endl;
	std::cout << "ids_to_include.size() " << ids_to_include.size() << std::endl;
	std::cout << "ids_to_exclude.size() " << ids_to_exclude.size() << std::endl;
	//check the query list
	if (ranges_to_include.nrow() > 0){
		StringVector const& chromosome = ranges_to_include["chromosome"] ;
                IntegerVector const& start = ranges_to_include["start"] ;
                IntegerVector const& end = ranges_to_include["end"] ;
		for( int i = 0; i < ranges_to_include.nrows(); ++i ) {
                        if( end[i] < start[i] ) {
                                throw std::invalid_argument( "Range (" + chromosome[i] + ":" + atoi( start[i] ) + "-" + atoi( end[i] ) + ") is malformed." ) ;
                         }
                        query->include_range( IndexQuery::GenomicRange( std::string( chromosome[i] ), start[i], end[i] )) ;
                }

	}

	if (ranges_to_exclude.nrow() > 0){

                StringVector const& chromosome_exclude = ranges_to_exclude["chromosome"] ;
                IntegerVector const& start_exclude = ranges_to_exclude["start"] ;
                IntegerVector const& end_exclude = ranges_to_exclude["end"] ;
                for( int i = 0; i < ranges_to_exclude.nrows(); ++i ) {
                        if( end_exclude[i] < start_exclude[i] ) {
                                throw std::invalid_argument( "Range (" + chromosome_exclude[i] + ":" + atoi( start_exclude[i] ) + "-" + atoi( end_exclude[i] ) + ") is malformed." ) ;
                         }
                        query->exclude_range( IndexQuery::GenomicRange( std::string( chromosome_exclude[i] ), start_exclude[i], end_exclude[i] )) ;
                }

        }

	if (ids_to_include.size() != 0){
		query->include_rsids(ids_to_include);

	}

	if (ids_to_exclude.size() != 0){
                query->exclude_rsids(ids_to_exclude);

        }

        query->initialise() ;
	if(query->number_of_variants() > 0){
          genoToTest_bgenDosage->set_query( query ) ;
	  numMarkers = genoToTest_bgenDosage->number_of_variants() ;
	  numSamples_bgen = genoToTest_bgenDosage->number_of_samples(); 
	  std::cout << numMarkers << " markers will be analyzed " << std::endl;
          return numMarkers ;
	}else{
          std::cout << "No queried variant is found in the bgen file! All variants bgen file will be analyzed" << std::endl;
          isQuery = false;
        }
      }

      if(!isQuery){  

          gm_stream.reset(
            new std::ifstream( filename.c_str(), std::ifstream::binary )
          ) ;

          if( !*gm_stream ) {
            throw std::invalid_argument( filename ) ;
          }

          //printf("1\n");fflush(NULL);
          gm_stream->seekg( 0, std::ios::beg ) ;
          genfile::bgen::read_offset( *gm_stream, &gm_offset ) ;
          //printf("2\n");fflush(NULL);   
          genfile::bgen::read_header_block( *gm_stream, &gm_context ) ;
	  
          uint Nbgen = gm_context.number_of_samples;
          numSamples_bgen = int(Nbgen);
          std::cout << numSamples_bgen << " samples are found in the bgen file" << std::endl;
	
          // Jump to the first variant data block.
          gm_stream->seekg( gm_offset + 4 ) ;
          //printf("4\n");fflush(NULL);
          uint Mbgen = gm_context.number_of_variants;
          numMarkers = int(Mbgen);
          //std::cout << "All " << numMarkers << " markers will be analyzed " << std::endl;
          std::cout << numMarkers << " markers are found in the bgen file " << std::endl;
          return numMarkers ;
		//numMarkers = genoToTest_bgenDosage->number_of_variants() ;
		//std::cout << "All " << numMarkers << " markers in the bgen file will be analyzed " << std::endl;
	}
   }

   //return numMarkers ;
 }
  

// [[Rcpp::export]]
Rcpp::List getDosage_inner_bgen_withquery(){

  using namespace genfile::bgen ;
  using namespace Rcpp ;

	//std::size_t max_entries_per_sample = 3;

	//std::size_t const number_of_variants = genoToTest_bgenDosage->number_of_variants() ;
  std::size_t const number_of_samples = genoToTest_bgenDosage->number_of_samples() ;

//        std::cout << "number_of_samples: " << number_of_samples << std::endl;

//        StringVector sampleNames( number_of_samples ) ;

//	view->get_sample_ids( set_sample_names( &sampleNames ) ) ;

  std::string SNPID, rsid, chromosome ;
  genfile::bgen::uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< double > dosages;
  genoToTest_bgenDosage->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;
  double dosage;

	//size_t variant_i = 0;
  int k;
  read_probs( &probs ) ;

/* with missing dosages
  for(std::size_t i = 0; i < probs.size(); ++i ) {
    dosage = 0;
    k = -1;
    for(std::size_t j = 0; j < probs[i].size(); ++j ) {
      k = k + 1;	
//                	std::cout << ( j > 0 ? "," : "" ) ;
      if( probs[i][j] == -1 ) {
        dosage = -1;	
      }else {
        dosage = dosage + probs[i][j]*k;
      }
    }
        dosages.push_back(dosage);
  }
*/
//assume no missing dosages, biallelic
  int N = probs.size();
  double sum_eij = 0, sum_fij_minus_eij2 = 0, fij, p10, p11;// for INFO
  for(std::size_t i = 0; i < N; ++i ) {  
    p10 = probs[i][1];
    p11 = probs[i][2];
    dosage = p10+p11*2;
    dosages.push_back(dosage);
/*
    if(i < 3){
        std::cout << "p11: " << p11 << std::endl;
        std::cout << "p10: " << p10 << std::endl;
        std::cout << "dosage: " << dosage << std::endl;
    }
*/

    sum_eij = sum_eij + dosage;
    fij = 4*p11 + p10;
    sum_fij_minus_eij2 = sum_fij_minus_eij2 + fij - dosage*dosage;
  }
  double thetaHat = sum_eij / (2*N);
  double info = thetaHat==0 || thetaHat==1 ? 1 :
  1 - sum_fij_minus_eij2 / (2*N*thetaHat*(1-thetaHat));
  markerInfo = info;

//  std::cout << "info: " << info << std::endl;
  DataFrame variants = DataFrame::create(
                Named("chromosome") = chromosome,
                Named("position") = position,
                Named("rsid") = rsid,
                Named("SNPID") = SNPID,
        //        Named("number_of_alleles") = number_of_allele,
                Named("allele0") = alleles[0],
                Named("allele1") = alleles[1],
		_["stringsAsFactors"] = false
  );


  List result ;
  result[ "variants" ] = variants ;
  result[ "dosages" ] = dosages ;

  dosages.clear();

  return( result ) ;
}

/*
// [[Rcpp::export]]
Rcpp::List getDosage_bgen_withquery()
{
        try {
		//return(getDosage_inner_bgen_withquery());
		return(getDosage_inner_bgen_withquery_new());
        }
        catch( std::exception const& e ) {
                forward_exception_to_r( e ) ;
        }
        catch( ... ) {
                ::Rf_error("A C++ exception occurred (unknown reason)") ;
        }
        return Rcpp::List() ;
}
*/

/*
        Nbgen: number of samples
*/


/**************************************
        This function is revised based on the Parse function in BOLT-LMM v2.3 source code
*************************************/

double  Parse(unsigned char * buf, size_t bufLen,  std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<int> & indexforMissing, double & homN_cases, double & hetN_cases, double & homN_ctrls, double & hetN_ctrls, double & ref_homN_cases, double & ref_homN_ctrls){

    size_t destLen = bufLen;

    unsigned char * bufAt = buf;
    uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;

    if (N != Nbgen) {
      std::cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << std::endl;
      exit(1);
    }
    uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
    if (K != 2U) {
      std::cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << std::endl;
      exit(1);
    }
    uint Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
      std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << std::endl;
      exit(1);
    }
    uint Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
      std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << std::endl;
      exit(1);
    }

    //deal with missing dosages
    std::vector <bool> missingIdxVec;
    missingIdxVec.clear();
    missingIdxVec.reserve(N);
    missingIdxVec.resize(N); 
    int missingSamplesize = 0;

    for (uint i = 0; i < N; i++) {
      uint ploidyMiss = *bufAt; bufAt++;
      bool const missing = (ploidyMiss & 0x80) ;
      missingIdxVec[i] = missing;
	if(missing){
          missingSamplesize = missingSamplesize + 1;
        }	
    }


    //for (uint i = 0; i < N; i++) {
    //  uint ploidyMiss = *bufAt; bufAt++;
    //  if (ploidyMiss != 2U) {
    //    std::cerr << "ERROR: " << snpName << " has ploidy/missingness byte = " << ploidyMiss
    //         << " (not 2)" << std::endl;
    //    exit(1);
    //  }
    // }
    uint Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
      std::cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << std::endl;
      exit(1);
    }
    uint B = *bufAt; bufAt++;
    if (B != 8U) {
      std::cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << std::endl;
      exit(1);
    }

        // Parse 
    double lut[256];
    for (int i = 0; i <= 255; i++)
      lut[i] = i/255.0;

    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0, sum_fij_minus_eij2_sub = 0; // for INFO
    double p11,p10,p00,dosage,eij,fij, eijsub, fijsub;
    dosages.clear();
    dosages.reserve(gmtest_samplesize);
    dosages.resize(gmtest_samplesize);
    std::size_t missing_cnt = 0;

    homN_cases = 0;
    hetN_cases = 0;
    ref_homN_cases = 0;
    homN_ctrls = 0;
    hetN_ctrls = 0;
    ref_homN_ctrls = 0;

   for (uint i = 0; i < N; i++) {
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      p00 = 1 - p11 - p10;
      if (dosage_type == 0) {
	dosage = p10+p00*2;
      } else if (dosage_type == 1) {
        dosage = p00;
      } else if (dosage_type == 2 ) {
        dosage = p10+p00;
      }      

      if(!missingIdxVec[i]){
        eij = p10+p00*2;
        fij = 4*p00 + p10;
        sum_eij += eij;
        sum_fij_minus_eij2 += fij - eij*eij;
        if(gm_sample_idx[i] >= 0){
          dosages[gm_sample_idx[i]] = dosage;
	  if(cc_idx[gm_sample_idx[i]] == 0) {
	    hetN_ctrls = hetN_ctrls + p10;
	    homN_ctrls = homN_ctrls + p00;
	    ref_homN_ctrls = ref_homN_ctrls + p11;
	  }
	  if (cc_idx[gm_sample_idx[i]] == 1) {
	    hetN_cases = hetN_cases + p10;
	    homN_cases = homN_cases + p00;
	    ref_homN_cases = ref_homN_cases + p11;
	  }
	  sum_eij_sub += dosage;
        }
     }else{
        if(gm_sample_idx[i] >= 0){        
          indexforMissing.push_back(gm_sample_idx[i]);
          ++missing_cnt;
          dosages[gm_sample_idx[i]] = -1; 
        }
     }
    //std::cout << "i: " <<  i << std::endl;
    }

   // counts can be negative zero so take abs
    homN_cases = std::abs(homN_cases); 
    hetN_cases = std::abs(hetN_cases);
    ref_homN_cases = std::abs(ref_homN_cases);
    homN_ctrls = std::abs(homN_ctrls);
    hetN_ctrls = std::abs(hetN_ctrls);    
    ref_homN_ctrls = std::abs(ref_homN_ctrls);

 
     AC = sum_eij_sub;


     if(gmtest_samplesize == missing_cnt){
       AF = 0;
     }else{
       AF = AC/ 2/ ((double) (gmtest_samplesize - missing_cnt)) ;
     }

     double thetaHat = sum_eij / (2* (N - missingSamplesize));
     double info = thetaHat==0 || thetaHat==1 ? 1 :
     1 - sum_fij_minus_eij2 / (2*N*thetaHat*(1-thetaHat));

     if(missing_cnt > 0){
       double imputeDosage = 2*AF;
       for (unsigned int i = 0; i < indexforMissing.size(); i++)
       {
          dosages[indexforMissing[i]] = imputeDosage;
       }
       if(!isDropMissingDosages_bgen){
	 std::cout << "AC is " << AC << std::endl;
         AC = AC + missing_cnt * imputeDosage;
	 std::cout << "AC_new is " << AC << std::endl;
       }	
     }



/*
    for (uint i = 0; i < gm_sample_idx.size(); i++) {
      int index = gm_sample_idx[i] *2;
      unsigned char * bufAt1 = &(bufAt[index]) ;	
      p11 = lut[*bufAt1]; bufAt1++;
      p10 = lut[*bufAt1];

      dosage = 2*p11 + p10;
      dosages.push_back(2 - dosage);
      eij = dosage;
      fij = 4*p11 + p10;
      sum_eij += eij;
      sum_fij_minus_eij2 += fij - eij*eij;
    }
*/


    return(info);
}


// [[Rcpp::export]]
Rcpp::List getDosage_inner_bgen_withquery_new(){

  using namespace genfile::bgen ;
  using namespace Rcpp ;

        //std::size_t max_entries_per_sample = 3;

        //std::size_t const number_of_variants = genoToTest_bgenDosage->number_of_variants() ;
  //std::size_t const number_of_samples = genoToTest_bgenDosage->number_of_samples() ;

//        std::cout << "number_of_samples: " << number_of_samples << std::endl;

//        StringVector sampleNames( number_of_samples ) ;

//      view->get_sample_ids( set_sample_names( &sampleNames ) ) ;

  std::string SNPID, rsid, chromosome ;
  genfile::bgen::uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< double > dosages;
  double AC, AF, homN_cases, hetN_cases, homN_ctrls, hetN_ctrls, ref_homN_cases, ref_homN_ctrls;
   
  //clock_t t1,t2;
  //t1=clock();
  isReadVariantBgen = genoToTest_bgenDosage->read_variant(&SNPID, &rsid, &chromosome, &position, &alleles ) ;
  //t2=clock();
  //float diff = ((float)t2-(float)t1);
  //float seconds = diff / CLOCKS_PER_SEC;
  //std::cout << seconds << std::endl;
  

  std::vector< genfile::byte_t > buffer2;


  //t1=clock();

  buffer2 = genoToTest_bgenDosage->read_and_uncompress_genotype_data_block();
  //t2=clock();
  //diff = ((float)t2-(float)t1);
  //seconds = diff / CLOCKS_PER_SEC;
  //std::cout << seconds << std::endl;


  unsigned char * buf  = (unsigned char *) buffer2.data();
  uint Nbgen = genoToTest_bgenDosage->number_of_samples();
  std::vector< int > indexforMissing;

  AC=0; AF=0; 
  homN_cases=0; hetN_cases=0; homN_ctrls=0; hetN_ctrls=0; ref_homN_cases=0; ref_homN_ctrls=0;
  markerInfo = Parse(buf, buffer2.size(), SNPID, Nbgen, dosages, AC, AF, indexforMissing, homN_cases, hetN_cases, homN_ctrls, hetN_ctrls, ref_homN_cases, ref_homN_ctrls);

  //t1=clock();
  //t2=clock();
  //diff =  ((float)t2-(float)t1);
  //seconds = diff / CLOCKS_PER_SEC;
  //std::cout << seconds << std::endl;

  //std::cout << "new way for probs" << std::endl;
	  DataFrame variants = DataFrame::create(
                Named("chromosome") = chromosome,
                Named("position") = position,
                Named("rsid") = rsid,
                Named("SNPID") = SNPID,
        //        Named("number_of_alleles") = number_of_allele,
                Named("allele0") = alleles[0],
                Named("allele1") = alleles[1],
                _["stringsAsFactors"] = false,
		Named("AC") = AC,
                Named("AF") = AF,
	        Named("ref_homN_cases") = ref_homN_cases,
                Named("hetN_cases") = hetN_cases,
	        Named("homN_cases") = homN_cases,
	        Named("ref_homN_ctrls") = ref_homN_ctrls,
                Named("hetN_ctrls") = hetN_ctrls,
		Named("homN_ctrls") = homN_ctrls
        ) ;


  List result ;
  result[ "variants" ] = variants ;
  result[ "dosages" ] = dosages ;
  result["indexforMissing"] = indexforMissing;

  indexforMissing.clear();
  dosages.clear();
  
  return(result);

}


// [[Rcpp::export]]
Rcpp::List getDosage_bgen_withquery()
{
        try {
//                return(getDosage_inner_bgen_withquery());
                return(getDosage_inner_bgen_withquery_new());
        }
        catch( std::exception const& e ) {
                forward_exception_to_r( e ) ;
        }
        catch( ... ) {
                ::Rf_error("A C++ exception occurred (unknown reason)") ;
        }
        return Rcpp::List() ;
}





/**************************************
        By SLEE 09/06/17
*************************************/

Rcpp::List getDosage_inner_bgen_noquery(){
  using namespace genfile;
  using namespace Rcpp ;

  bool temp;
  std::string SNPID, RSID, chromosome, first_allele,second_allele ;
  genfile::bgen::uint32_t position;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< byte_t > buffer1;
  std::vector< byte_t > buffer2;
  std::vector< double > dosages;
  double AC, AF, homN_cases, hetN_cases, homN_ctrls, hetN_ctrls, ref_homN_cases, ref_homN_ctrls;

  isReadVariantBgen = genfile::bgen::read_snp_identifying_data(
                        *gm_stream,
                        gm_context,
                        &SNPID,
                        &RSID,
                        &chromosome,
                        &position,
                        &first_allele,
                        &second_allele
                ) ;



  genfile::bgen::read_genotype_data_block(
                        *gm_stream,
                        gm_context,
                        &buffer1
                ) ;



  genfile::bgen::uncompress_probability_data(
                        gm_context,
                        buffer1,
                        &buffer2
                ) ;



  unsigned char * buf  = (unsigned char *) buffer2.data();
  uint Nbgen = gm_context.number_of_samples;
  AC=0; AF=0; 
  homN_cases=0; hetN_cases=0; homN_ctrls=0; hetN_ctrls=0; ref_homN_cases=0; ref_homN_ctrls=0;
  std::vector< int > indexforMissing;
  markerInfo = Parse(buf, buffer2.size(), SNPID, Nbgen, dosages, AC, AF, indexforMissing, homN_cases, hetN_cases, homN_ctrls, hetN_ctrls, ref_homN_cases, ref_homN_ctrls);

  DataFrame variants = DataFrame::create(
                Named("chromosome") = chromosome,
                Named("position") = position,
                Named("rsid") = RSID,
		Named("SNPID") = SNPID,
        //        Named("number_of_alleles") = number_of_allele,
                Named("allele0") = first_allele,
                Named("allele1") = second_allele,
                _["stringsAsFactors"] = false,
		Named("AC") = AC,
                Named("AF") = AF,
	        Named("ref_homN_cases") = ref_homN_cases,
                Named("hetN_cases") = hetN_cases,
	        Named("homN_cases") = homN_cases,
	        Named("ref_homN_ctrls") = ref_homN_ctrls,
                Named("hetN_ctrls") = hetN_ctrls,
		Named("homN_ctrls") = homN_ctrls
        ) ;

  List result ;
  result[ "variants" ] = variants ;
  result[ "dosages" ] = dosages ;
  result["indexforMissing"] = indexforMissing;

  dosages.clear();
  indexforMissing.clear();

  return(result);
}


// [[Rcpp::export]]
Rcpp::List getDosage_bgen_noquery()
{
  try {
    return(getDosage_inner_bgen_noquery());
  }
  catch( std::exception const& e ) {
    forward_exception_to_r( e ) ;
  }
  catch( ... ) {
    ::Rf_error("A C++ exception occurred (unknown reason)") ;
  }
  return Rcpp::List() ;
}

// [[Rcpp::export]]
bool getQueryStatus()
{
  return(isQuery);
}

// [[Rcpp::export]]
bool getisReadVariantBgen()
{
  return(isReadVariantBgen);
}



// [[Rcpp::export]]
double getMarkerInfo()
{
  return(markerInfo);
}

/**************************************
        By SLEE 09/06/17
*************************************/

// [[Rcpp::export]]
void SetSampleIdx(Rcpp::IntegerVector sample_idx, Rcpp::IntegerVector cc_index, int Ntest){
	gmtest_samplesize = Ntest;
	gm_sample_idx = sample_idx;
	cc_idx = cc_index;
}


// [[Rcpp::export]]
void closetestGenoFile_bgenDosage() //needs further check
{
  
//  gm_sample_idx.erase();

  printf("closed the genofile!\n");

}


//for gene-based test

// [[Rcpp::export]]
int setgenoTest_bgenDosage_v2(std::string & filename,
        std::string & index_filename,
        Rcpp::DataFrame & ranges_to_include,
        Rcpp::DataFrame & ranges_to_exclude,
        std::vector< std::string > const& ids_to_include,
        std::vector< std::string > const& ids_to_exclude
){

   //bgenMinMAF = bgenMinMaf;
   //bgenMinINFO = bgenMinInfo;
   int numMarkers;
   {
     using namespace genfile::bgen ;
     using namespace Rcpp;
     isQuery = false;
     if(index_filename == ""){
       isQuery = false;
       std::cout << "no index file for bgen is provided" << std::endl;
     }else if(ranges_to_include.nrow() == 0 &&
          ranges_to_exclude.nrow() == 0 &&
          ids_to_include.size() == 0 &&
          ids_to_exclude.size() == 0){
       std::cout << "no query list is provided" << std::endl;
       isQuery = false;
     }else{
       isQuery = true;
     }


     if(isQuery){
        std::cout << "TEST " << ids_to_include.size() << std::endl;
        genoToTest_bgenDosage = View::create( filename ) ;
        std::cout << "TEST 1 OK" << std::endl;
        IndexQuery::UniquePtr query = IndexQuery::create( index_filename ) ;
        std::cout << "TEST2 " << ids_to_include.size() << std::endl;

        std::cout << "ranges_to_include.nrow() " << ranges_to_include.nrow() << std::endl;
        std::cout << "ranges_to_exclude.nrow() " << ranges_to_exclude.nrow() << std::endl;
        std::cout << "ids_to_include.size() " << ids_to_include.size() << std::endl;
        std::cout << "ids_to_exclude.size() " << ids_to_exclude.size() << std::endl;
        //check the query list
        if (ranges_to_include.nrow() > 0){
                StringVector const& chromosome = ranges_to_include["chromosome"] ;
                IntegerVector const& start = ranges_to_include["start"] ;
                IntegerVector const& end = ranges_to_include["end"] ;
                for( int i = 0; i < ranges_to_include.nrows(); ++i ) {
                        if( end[i] < start[i] ) {
                                throw std::invalid_argument( "Range (" + chromosome[i] + ":" + atoi( start[i] ) + "-" + atoi( end[i] ) + ") is malformed." ) ;
                         }
                        query->include_range( IndexQuery::GenomicRange( std::string( chromosome[i] ), start[i], end[i] )) ;
                }

        }

        if (ranges_to_exclude.nrow() > 0){

                StringVector const& chromosome_exclude = ranges_to_exclude["chromosome"] ;
                IntegerVector const& start_exclude = ranges_to_exclude["start"] ;
                IntegerVector const& end_exclude = ranges_to_exclude["end"] ;
                for( int i = 0; i < ranges_to_exclude.nrows(); ++i ) {
                        if( end_exclude[i] < start_exclude[i] ) {
                                throw std::invalid_argument( "Range (" + chromosome_exclude[i] + ":" + atoi( start_exclude[i] ) + "-" + atoi( end_exclude[i] ) + ") is malformed." ) ;
                         }
                        query->exclude_range( IndexQuery::GenomicRange( std::string( chromosome_exclude[i] ), start_exclude[i], end_exclude[i] )) ;
                }

        }

        if (ids_to_include.size() != 0){
                query->include_rsids(ids_to_include);

        }

        if (ids_to_exclude.size() != 0){
                query->exclude_rsids(ids_to_exclude);

        }

        query->initialise() ;
        if(query->number_of_variants() > 0){
          genoToTest_bgenDosage->set_query( query ) ;
          numMarkers = genoToTest_bgenDosage->number_of_variants() ;
          std::cout << numMarkers << " markers will be analyzed " << std::endl;
          return numMarkers ;
        }else{
          std::cout << "No queried variant is found in the bgen file! All variants bgen file will be analyzed" << std::endl;
//          isQuery = false;
	 genoToTest_bgenDosage->set_query( query ) ;
	 numMarkers = genoToTest_bgenDosage->number_of_variants() ;
	 return numMarkers ;
        }
      }

      if(!isQuery){

          gm_stream.reset(
            new std::ifstream( filename.c_str(), std::ifstream::binary )
          ) ;

          if( !*gm_stream ) {
            throw std::invalid_argument( filename ) ;
          }

          //printf("1\n");fflush(NULL);
          gm_stream->seekg( 0, std::ios::beg ) ;
          genfile::bgen::read_offset( *gm_stream, &gm_offset ) ;
          //printf("2\n");fflush(NULL);
          genfile::bgen::read_header_block( *gm_stream, &gm_context ) ;

          uint Nbgen = gm_context.number_of_samples;
          int numSamples = int(Nbgen);
          std::cout << numSamples << " samples are found in the bgen file" << std::endl;

          // Jump to the first variant data block.
          gm_stream->seekg( gm_offset + 4 ) ;
          //printf("4\n");fflush(NULL);
          uint Mbgen = gm_context.number_of_variants;
          numMarkers = int(Mbgen);
          //std::cout << "All " << numMarkers << " markers will be analyzed " << std::endl;
          std::cout << numMarkers << " markers are found in the bgen file " << std::endl;
          return numMarkers ;
                //numMarkers = genoToTest_bgenDosage->number_of_variants() ;
                //std::cout << "All " << numMarkers << " markers in the bgen file will be analyzed " << std::endl;
        }
   }

   //return numMarkers ;
 }


// [[Rcpp::export]]
int getSampleSizeinBgen(){
        return(numSamples_bgen);
}
