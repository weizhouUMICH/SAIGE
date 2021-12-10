
// All C++ codes related to PLINK file manipulation

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <boost/algorithm/string.hpp>
#include "PLINK.hpp"
#include "UTIL.hpp"

#include <RcppArmadilloExtensions/sample.h> // sample
#include <string>

// make a global variable for future usage
// static PLINK::PlinkClass* ptr_gPLINKobj = NULL;

// ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
//                                       t_famFile,
//                                       t_bedFile,
//                                       t_SampleInModel);

namespace PLINK {

PlinkClass::PlinkClass(std::string t_bimFile,
                       std::string t_famFile,
                       std::string t_bedFile,
                       std::vector<std::string> t_SampleInModel,
                       std::string t_AlleleOrder)
{
  setPlinkobj(t_bimFile, t_famFile, t_bedFile);
  setPosSampleInPlink(t_SampleInModel);
  m_AlleleOrder = t_AlleleOrder;
}

// set PlinkClass by reading plink files
void PlinkClass::setPlinkobj(std::string t_bimFile,
                             std::string t_famFile,
                             std::string t_bedFile)
{
  m_bimFile = t_bimFile;
  m_famFile = t_famFile;
  m_bedFile = t_bedFile;
  
  // setChrMaps();
  readBimFile();
  readFamFile();
  m_ibedFile.open(m_bedFile.c_str(), std::ios::binary);
  
  m_ibedFile.seekg(2);
  char magicNumber3;
  m_ibedFile.read(&magicNumber3, 1);
  
  if(magicNumber3 != 1)
    Rcpp::stop("The third magic number of the plink bed file is not 00000001. Please use SNP-major plink (plink version >= 1.9) files.");
}

// not used after 2021-02-28
// void PlinkClass::setChrMaps()
// {
//   uint8_t last_chr_autosome = 22;
//   // Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
//   for(uint8_t index = 0; index <= 26; index++){
//     m_chrMaps[std::to_string(index)] = index;
//   }
//   m_chrMaps["X"] = last_chr_autosome + 1;
//   m_chrMaps["x"] = last_chr_autosome + 1;
//   m_chrMaps["Y"] = last_chr_autosome + 2;
//   m_chrMaps["y"] = last_chr_autosome + 2;
//   m_chrMaps["XY"] = last_chr_autosome + 3;
//   m_chrMaps["xy"] = last_chr_autosome + 3;
//   m_chrMaps["MT"] = last_chr_autosome + 4;
//   m_chrMaps["mt"] = last_chr_autosome + 4;
// }    

void PlinkClass::readBimFile()
{
  std::cout << "Reading bim file...." << std::endl;
  std::ifstream bim(m_bimFile);
  m_M0 = 0;
  std::string line;
  // uint8_t chr_item;
  while(getline(bim, line)){
    m_M0++;
    std::vector<std::string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
    // chr_item = m_chrMaps.at(line_elements[0]);
    m_chr.push_back(line_elements[0]);
    m_MarkerInPlink.push_back(line_elements[1]);
    m_gd.push_back(std::stof(line_elements[2]));
    m_pd.push_back(std::stoi(line_elements[3]));
    std::transform(line_elements[4].begin(), line_elements[4].end(), line_elements[4].begin(), toupper);
    std::transform(line_elements[5].begin(), line_elements[5].end(), line_elements[5].begin(), toupper);
    m_alt.push_back(line_elements[4]);  // allele 1, usually minor allele
    m_ref.push_back(line_elements[5]);  // allele 2, usually major allele
  }
  m_M = m_M0;
}

void PlinkClass::readFamFile()
{
  std::cout << "Reading fam file...." << std::endl;
  std::ifstream fam(m_famFile);
  m_N0 = 0;
  std::string line;
  while(getline(fam, line)){
    m_N0 ++;
    std::vector<std::string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    // boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
    m_SampleInPlink.push_back(line_elements[1]);
  }
  m_N = m_N0;
  m_numBytesofEachMarker0 = (m_N0 + 3) / 4;
  m_OneMarkerG4.reserve(m_numBytesofEachMarker0);  
  m_OneMarkerG4.resize(m_numBytesofEachMarker0);
} 


void PlinkClass::setPosSampleInPlink(std::vector<std::string> t_SampleInModel)
{
  std::cout << "Setting position of samples in plink files...." << std::endl;
  m_N = t_SampleInModel.size();
  m_numBytesofEachMarker = (m_N + 3) / 4;
  
  Rcpp::CharacterVector SampleInModel(m_N);
  for(uint32_t i = 0; i < m_N; i++)
    SampleInModel(i) = t_SampleInModel.at(i);
  
  uint32_t N_plink = m_SampleInPlink.size();
  Rcpp::CharacterVector SampleInPlink(N_plink);
  for(uint32_t i = 0; i < N_plink; i++)
    SampleInPlink(i) = m_SampleInPlink.at(i);
    
  Rcpp::IntegerVector posSampleInPlink = Rcpp::match(SampleInModel, SampleInPlink);
  m_posSampleInPlink.resize(m_N);
  for(uint32_t i = 0; i < m_N; i++){
    if(Rcpp::IntegerVector::is_na(posSampleInPlink.at(i)))
      Rcpp::stop("At least one subject requested is not in Plink file.");
    m_posSampleInPlink.at(i) = posSampleInPlink.at(i) - 1;   // convert "starting from 1" to "starting from 0"
  }
  
  // m_posSampleInPlink.clear();
  // for(uint32_t i = 0; i < m_N; i++){
  //   std::string sample = t_SampleInModel.at(i);
  //   auto pos = std::find(m_SampleInPlink.begin(), m_SampleInPlink.end(), sample);
  //   if(pos != m_SampleInPlink.end()){
  //     m_posSampleInPlink.push_back(pos - m_SampleInPlink.begin());
  //   }else{
  //     Rcpp::stop("At least one subject requested is not in Plink file.");
  //   }
  // }
}

std::vector<uint32_t> PlinkClass::getPosMarkerInPlink(std::vector<std::string> t_MarkerReqstd)
{
  int M = t_MarkerReqstd.size();
  std::vector<uint32_t> posMarkerInPlink;
  for(int i = 0; i < M; i++){
    std::string marker = t_MarkerReqstd.at(i);
    auto pos = std::find(m_MarkerInPlink.begin(), m_MarkerInPlink.end(), marker);
    if(pos != m_MarkerInPlink.end()){
      posMarkerInPlink.push_back(pos - m_MarkerInPlink.begin());
    }else{
      Rcpp::warning("Marker %s is not found in plink file.", marker);
    }
  }
  return posMarkerInPlink;
}

bool PlinkClass::getOneMarker(uint64_t t_gIndex,        // different meanings for different genoType
                                   std::string& t_ref,       // REF allele
                                   std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                   std::string& t_marker,    // marker ID extracted from genotype file
                                   uint32_t& t_pd,           // base position
                                   std::string& t_chr,       // chromosome
                                   double& t_altFreq,        // frequency of ALT allele
                                   double& t_altCounts,      // counts of ALT allele
                                   double& t_mac, 
                                   double& t_missingRate,    // missing rate
                                   double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                                   bool t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                                   std::vector<uint32_t>& t_indexForMissing,     // index of missing genotype data
                                   bool t_isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                   std::vector<uint32_t>& t_indexForNonZero,
                                   bool t_isTrueGenotype,    // only used in PLINK. check m_genoMaps for details about the genotype mapping in PLINK.
				   std::vector<double> & t_dosage)	
{
  int sum = 0;
  int numMissing = 0;
  
  if(!t_isOnlyOutputNonZero){
    t_dosage.resize(m_N);
  }
  
  uint64_t posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex;
  m_ibedFile.seekg(posSeek);
  m_ibedFile.read((char*)(&m_OneMarkerG4[0]), m_numBytesofEachMarker0);
  
  t_indexForMissing.clear();
  t_ref = m_ref[t_gIndex];
  t_alt = m_alt[t_gIndex];
  t_marker = m_MarkerInPlink[t_gIndex];
  t_pd = m_pd[t_gIndex];
  t_chr = m_chr[t_gIndex];
  
  for(uint32_t i = 0; i < m_N; i++){
    uint32_t ind = m_posSampleInPlink[i];             // C++ start from 0
    unsigned char bufferG4 = m_OneMarkerG4[ind/4];    // unsigned char: 1 byte for 4 genotypes (4 samples)
    int bufferG1;                                     // int: 1 genotype (1 sample)
    getGenotype(&bufferG4, ind%4, bufferG1);          // bufferG4 -> bufferG1
    switch(bufferG1){
    case HOM_REF: break;
    case HET: sum+=1; break;
    case HOM_ALT: sum+=2; break;
    case MISSING: numMissing++; if(t_isOutputIndexForMissing){t_indexForMissing.push_back(i);} break;  
    }
    
    // OneMarkerG1[i] = bufferG1;
    if(t_isTrueGenotype){
      bufferG1 = m_genoMaps[bufferG1];
    }
    
    if(t_isOnlyOutputNonZero){
      t_dosage.push_back(bufferG1);
    }else{
      t_dosage.at(i) = bufferG1;
    }

  }
  
  int count = m_N - numMissing;
  t_missingRate = (double)numMissing / (double)m_N;
  t_imputeInfo = 1;
  t_altCounts = (double)sum;
  t_altFreq = t_altCounts/ (double)count / 2;
  
  // updated on 03/14/2021
  if(m_AlleleOrder == "ref-first"){
    t_ref = m_alt[t_gIndex];
    t_alt = m_ref[t_gIndex];
    t_altFreq = 1 - t_altFreq;
    t_altCounts = 2 * (double)count * t_altFreq;
    for(unsigned int i = 0; i < t_dosage.size(); i++)
      t_dosage.at(i) = 2 - t_dosage.at(i);
  }

  t_mac = t_altCounts;
  if(t_altFreq > 0.5){
    t_mac = (double)count * 2 - t_altCounts;
  }	  
  
  return true; //need a bool variabke to show whether the function is successfuly run
}

std::vector<unsigned int> whichCPP(std::vector<std::string> strVec, 
                                   std::string strValue)
{
  std::vector<unsigned int> indexVec;
  for(unsigned int i = 0; i < strVec.size(); i++){
    if(strVec.at(i)==strValue)
      indexVec.push_back(i);
  }
  return(indexVec);
}

