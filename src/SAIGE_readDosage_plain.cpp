// Important: this definition ensures Armadillo enables SuperLU
#define ARMA_USE_SUPERLU 1

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <bitset>

#include "../thirdParty/bgen/3rd_party/boost_1_55_0/boost/iostreams/filtering_stream.hpp"
#include "../thirdParty/bgen/3rd_party/boost_1_55_0/boost/iostreams/filter/gzip.hpp"
#include "../thirdParty/bgen/3rd_party/boost_1_55_0/boost/iostreams/filter/bzip2.hpp"


using namespace Rcpp;
using namespace std;
using namespace RcppParallel;



//This is a class with attritbutes about the genotype to test
class genoTestClass_plainDosage{
public:

  boost::iostreams::filtering_istream test_genoGZfile;
  std::ifstream test_genofile;
  bool isGZfile;
  int Mtest;
  int NSampleTest;
  arma::fvec  oneSNPGenoTestVec;
  std::vector<std::string> rowHeaderVec;


  int setGenoTestObjGZ(std::string testGenoFile, int testGenofileNrowSkip, int testGenofileNcolSkip){
    //printf("setGenoTestObjGZ!\n");
    printf("set GenoTestObj from dosage gz file!\n");	
    std::ifstream test_genoGZfile0;
    test_genoGZfile0.open(testGenoFile.c_str(), std::ios_base::in | std::ios_base::binary); //add

    if (!test_genoGZfile0.is_open()){
      printf("Error! genofile to test cannot be openned!");
    }

    boost::iostreams::filtering_istream readfile;
    readfile.push(boost::iostreams::gzip_decompressor());
    readfile.push(test_genoGZfile0);



    std::string line;
    std::getline(readfile,line);
    //std::getline(test_genofile,line);
    std::istringstream iss(line);

    int colNum = 0;
     do {
        std::string sub;
        iss >> sub;
        //cout << "sub is " << sub << endl;
        if (sub.length())
          ++colNum;
        } while(iss);


     readfile.clear(); //add
     test_genoGZfile0.close();  //add
     //test_genofile.clear();  //add
     //number of samples N = colNum - genofileNcolSkip;
     NSampleTest = colNum - testGenofileNcolSkip; 
     cout << "NSampleTest " << NSampleTest << endl;

     //printf("setGenoTestObjGZ2!\n");
     //test_genofile.open(testGenoFile.c_str());
     

     std::ifstream test_genoGZfile2a;	
     test_genoGZfile2a.open(testGenoFile.c_str(), ios_base::in | ios_base::binary); //add
     if (!test_genoGZfile2a.is_open()){
       printf("Error! genofile to test cannot be openned!");
     }

     boost::iostreams::filtering_istream test_genoGZfile2b;
     test_genoGZfile2b.push(boost::iostreams::gzip_decompressor()); //add
     test_genoGZfile2b.push(test_genoGZfile2a); //add     
     //printf("setGenoTestObjGZ3!\n");
     std::string junk;
     int indexRow = 0;
     int buffer;
     int TotalRead=0;

     /////////////////////////////
      // Added for reserve
      //while (std::getline(test_genofile,junk)){
     while (std::getline(test_genoGZfile2b,junk)){  //add
     //while (std::getline(test_genofile,junk)){
        indexRow = indexRow + 1;
        junk.clear();
     }
     TotalRead= indexRow;
     Mtest = indexRow - testGenofileNrowSkip;
     cout << "Mtest " << Mtest << endl;
      //oneSNPGenoTestVec.set_size(Mtest);
     oneSNPGenoTestVec.set_size(NSampleTest);
     test_genoGZfile2b.clear(); //add
     test_genoGZfile2a.close();  //add
      //test_genofile.clear();  //add

      //test_genofile.close();
     // test_genofile.open(testGenoFile.c_str());

     test_genofile.open(testGenoFile.c_str(), ios_base::in | ios_base::binary); //add
     //printf("setGenoTestObjGZ4!\n");
     if (!test_genofile.is_open()){
       printf("Error! genofile to test cannot be openned!");
     }

     test_genoGZfile.push(boost::iostreams::gzip_decompressor()); //add
     test_genoGZfile.push(test_genofile); 

     cout << "\nTesting markers: " << Mtest << " , samples: " << NSampleTest<< endl;

     std::string sbuffer;
   
     if(testGenofileNrowSkip > 0){
       for (int i = 0; i < testGenofileNrowSkip; i++){
          sbuffer.clear();
          std::getline(test_genoGZfile, sbuffer);
       }
     }
     sbuffer.clear();
     //std::getline(test_genoGZfile, sbuffer);
     //std::string str2 = sbuffer.substr (10,500);
     //cout << "sbuffer: " << sbuffer << endl;	

     //test_genoGZfile.empty(); //add
//     test_genofile.close();  //add
     return(Mtest);
  }




  int setGenoTestObj(std::string testGenoFile, int testGenofileNrowSkip, int testGenofileNcolSkip)  {

    //printf("setGenoTestObj!\n");
    printf("set GenoTestObj from dosage file!\n");


    //test_genofile.open(testGenoFile.c_str(), std::ios_base::in | std::ios_base::binary); //add
    test_genofile.open(testGenoFile.c_str());
    if (!test_genofile.is_open()){
      printf("Error! genofile to test cannot be openned!");
    }
    //boost::iostreams::filtering_istream readfile;  //add
    //readfile.push(boost::iostreams::gzip_decompressor()); //add
    //readfile.push(test_genofile); //add
    std::string line;
    //std::getline(readfile,line); //add
    std::getline(test_genofile,line);
    std::istringstream iss(line);
    int colNum = 0;
     do {
        std::string sub;
        iss >> sub;
        //cout << "sub is " << sub << endl;
        if (sub.length())
          ++colNum;
        } while(iss);

     //readfile.clear(); //add
     test_genofile.close();  //add
     //test_genofile.clear();  //add
     //number of samples N = colNum - genofileNcolSkip;
     NSampleTest = colNum - testGenofileNcolSkip;
     cout << "NSampleTest " << NSampleTest << endl;

     //printf("setGenoTestObj2!\n");    
     test_genofile.open(testGenoFile.c_str());
     //test_genofile.open(testGenoFile.c_str(), ios_base::in | ios_base::binary); //add
     if (!test_genofile.is_open()){
       printf("Error! genofile to test cannot be openned!");
     }
     //boost::iostreams::filtering_istream readfile2;
     //readfile2.push(boost::iostreams::gzip_decompressor()); //add
     //readfile2.push(test_genofile); //add     
     //printf("setGenoTestObj3!\n");     
     std::string junk;
     int indexRow = 0;
     int buffer;
     int TotalRead=0;

      /////////////////////////////
      // Added for reserve
      //while (std::getline(test_genofile,junk)){
      //while (std::getline(readfile2,junk)){  //add
      while (std::getline(test_genofile,junk)){
	indexRow = indexRow + 1;
        junk.clear();
      }
      TotalRead= indexRow;
      Mtest = indexRow - testGenofileNrowSkip;
      cout << "Mtest " << Mtest << endl;
      //oneSNPGenoTestVec.set_size(Mtest);
      oneSNPGenoTestVec.set_size(NSampleTest);
      //readfile2.clear(); //add
      test_genofile.close();  //add
      //test_genofile.clear();  //add
      
      //test_genofile.close();
      test_genofile.open(testGenoFile.c_str());
      //test_genofile.open(testGenoFile.c_str(), ios_base::in | ios_base::binary); //add
      //printf("setGenoTestObj4!\n");
      if (!test_genofile.is_open()){
       printf("Error! genofile to test cannot be openned!\n");
     }
     //test_genoGZfile.push(boost::iostreams::gzip_decompressor()); //add
     //test_genoGZfile.push(test_genofile); 
     cout << "\nTesting markers: " << Mtest << " , samples: " << NSampleTest<< endl;


     std::string sbuffer;

     for (int i = 0; i < testGenofileNrowSkip; i++){
	sbuffer.clear();
        std::getline(test_genofile, sbuffer);
     }
     
     return(Mtest);
  }

};


//A global variable for the geno file to test
genoTestClass_plainDosage genoToTest_plainDosage;
// Added by SLEE for parsing 09/07/2017
Rcpp::IntegerVector gm_sample_idx_plainDosage;
int gmtest_samplesize_plainDosage;

// [[Rcpp::export]]
int setgenoTest_plainDosage(std::string testGenoFile, int testGenofileNrowSkip, int testGenofileNcolSkip)
{
  int Mmarkers;
  cout << "setgenoTest here\n" << endl;

  //check the file extension. Is it a gz file?
    std::string fileExt;
    if (testGenoFile.find_last_of(".") != std::string::npos){
        fileExt = testGenoFile.substr(testGenoFile.find_last_of(".")+1);
        cout << "file extenstion is " << fileExt << endl;
    }

    if (fileExt != "gz"){
        genoToTest_plainDosage.isGZfile = false;
  	Mmarkers = genoToTest_plainDosage.setGenoTestObj(testGenoFile, testGenofileNrowSkip, testGenofileNcolSkip);
    }else{
        genoToTest_plainDosage.isGZfile = true;
        Mmarkers = genoToTest_plainDosage.setGenoTestObjGZ(testGenoFile, testGenofileNrowSkip, testGenofileNcolSkip);
    }	
  return(Mmarkers);
}



// [[Rcpp::export]]
void closetestGenoFile_plainDosage()
{
  //genoToTest_plainDosage.test_genoGZfile.close();

  if(!genoToTest_plainDosage.isGZfile){
    genoToTest_plainDosage.test_genofile.close();//add

  }else{
    genoToTest_plainDosage.test_genoGZfile.clear();
    genoToTest_plainDosage.test_genofile.close();

  }

  genoToTest_plainDosage.oneSNPGenoTestVec.clear();
  genoToTest_plainDosage.rowHeaderVec.clear();
  printf("closed the genofile!\n");

}

//arma::fvec getGenoOfnthVar_plainDosage(int mth, int testGenofileNrowSkip, int testGenofileNcolSkip, std::string dosageFileMissing) {



// [[Rcpp::export]]
arma::fvec getGenoOfnthVar_plainDosage(int mth, int testGenofileNrowSkip, int testGenofileNcolSkip) {
  	genoToTest_plainDosage.oneSNPGenoTestVec.clear();
  	genoToTest_plainDosage.oneSNPGenoTestVec.set_size(gmtest_samplesize_plainDosage);


  	float buffer;
  	std::string junk;
  	std::string line;
  	std::string subbuffer;
  	int colCount = 0;
  	genoToTest_plainDosage.rowHeaderVec.clear();

  	if(!genoToTest_plainDosage.isGZfile){
   		for (int col=0; col < testGenofileNcolSkip; col++){
      			junk.clear();
      			genoToTest_plainDosage.test_genofile >> junk;
      			genoToTest_plainDosage.rowHeaderVec.push_back(junk);
    		}
   	 	for (int col=0; col< genoToTest_plainDosage.NSampleTest; col++){
			junk.clear();
			genoToTest_plainDosage.test_genofile >> junk;
			std:istringstream iss(junk);
        		iss >> buffer;
			if(gm_sample_idx_plainDosage[col] >= 0){
				genoToTest_plainDosage.oneSNPGenoTestVec(gm_sample_idx_plainDosage[col]) = buffer;
			}
    		}

  	}else{
    		std::getline(genoToTest_plainDosage.test_genoGZfile,line);

    		stringstream ss(line);

    		for (int col=0; col < testGenofileNcolSkip; col++){
      			ss >> junk;
      			if( ss.fail() ) break;
      			genoToTest_plainDosage.rowHeaderVec.push_back(junk);
    		} 

    		for (int col=0; col< genoToTest_plainDosage.NSampleTest; col++){

      			ss >> buffer;
      			if( ss.fail() ) break;
			if(gm_sample_idx_plainDosage[col] >= 0){
                		genoToTest_plainDosage.oneSNPGenoTestVec(gm_sample_idx_plainDosage[col]) = buffer;
        		}
    		}

    	line.clear();
    	junk.clear();
    	buffer = 0;

   	}

  	return(genoToTest_plainDosage.oneSNPGenoTestVec);
}


// [[Rcpp::export]]
std::vector<std::string> getrowHeaderVec_plainDosage(){
	return(genoToTest_plainDosage.rowHeaderVec);
}


// [[Rcpp::export]]
int SetSampleIdx_plainDosage(Rcpp::IntegerVector sample_idx, int Ntest){
        gmtest_samplesize_plainDosage = Ntest;
        gm_sample_idx_plainDosage = sample_idx;

}
