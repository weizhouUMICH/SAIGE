//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h> 
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>// include this header for calculating execution time 
#include <cassert>
using namespace Rcpp;
using namespace std;
using namespace RcppParallel;


//This is a class with attritbutes about the genotype informaiton 
class genoClass{
private:
        //COPY from RVTEST:
        // we reverse the two bits as defined in PLINK format
        const static unsigned char HOM_REF = 0x0;  // 0b00 ;
        const static unsigned char HET = 0x2;      // 0b10 ;
        const static unsigned char HOM_ALT = 0x3;  // 0b11 ;
        const static unsigned char MISSING = 0x1;  // 0b01 ;


public:
        //to chunk the geno vector to avoid large continuous memory usage 
	int numMarkersofEachArray;
        int numofGenoArray;
        int numMarkersofLastArray;
        std::vector< std::vector<unsigned char>* > genoVecofPointers;
        ///////////

	//vector<unsigned char> genoVec; 	 
  	size_t M;
  	size_t N;
	size_t Nnomissing;
  	vector<float>	invstdvVec;
	vector<int>	ptrsubSampleInGeno;
	
  	arma::fvec 	alleleFreqVec;
  	arma::fvec	m_OneSNP_Geno;
  	arma::fvec	m_OneSNP_StdGeno;
  	arma::fvec	m_DiagStd;

	unsigned char m_genotype_buffer[4];
	int geno_idx;
	int m_size_of_esi;
	unsigned char m_bits_val[8];


        void setBit(unsigned char & ch, int ii, int aVal, int bVal){

                if (bVal == 1 && aVal == 1){
			ch ^= char(1 << ((ii*2) + 1)); //set a to be 1

                }else if(bVal == 0){
			ch ^= char(1 << (ii*2)); //change b to 0

                        if(aVal == 1){
				ch ^= char(1 << ((ii*2) + 1)); //change a to 1
                        }
                }
        }



	//COPY from RVTEST:
	void setGenotype(unsigned char* c, const int pos, const int geno) {
    		(*c) |= (geno << (pos << 1));
  	}





	void Init_OneSNP_Geno(){
		m_size_of_esi = (Nnomissing+3)/4;
		int k = 8;
		while (k > 0){
			-- k;
			m_bits_val[k] = 1 << k;
		}
	}
	

        arma::fvec * Get_OneSNP_Geno(size_t SNPIdx){
                m_OneSNP_Geno.zeros(Nnomissing);

		//avoid large continuous memory usage
		int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
		////////////////

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind= 0;
                unsigned char geno1;
                int bufferGeno;
                for(size_t i=Start_idx; i< Start_idx+m_size_of_esi; i++){
                        //geno1 = genoVec[i];
			geno1 = genoVecofPointers[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
                        for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
				bufferGeno = 2-(a+b);
				m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                if(ind >= Nnomissing){

                                //printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                        return & m_OneSNP_Geno;
                                }
                        }
                }

                return & m_OneSNP_Geno;
       }
   
	arma::fvec * Get_OneSNP_Geno_atBeginning(size_t SNPIdx, vector<int> & indexNA, vector<unsigned char> & genoVecOneMarkerOld){

		arma::fvec m_OneSNP_GenoTemp;
		m_OneSNP_GenoTemp.zeros(N);
		m_OneSNP_Geno.zeros(Nnomissing);
		int m_size_of_esi_temp = (N+3)/4;
		size_t ind= 0;
		unsigned char geno1;
		int bufferGeno;
		for(size_t i=0; i< m_size_of_esi_temp; i++){
			geno1 = genoVecOneMarkerOld[i];
			for(int j=0; j<4; j++){
				int b = geno1 & 1 ;
				geno1 = geno1 >> 1;
				int a = geno1 & 1 ;
				if (b == 1 && a == 0){
                                        bufferGeno = 3;
                                }else if(b == 0 && a == 0){
                                        bufferGeno = 2;
                                }else if(b == 0 && a == 1){
                                        bufferGeno = 1;
                                }else if(b == 1 && a == 1){
                                        bufferGeno = 0;
                                }else{
                                        cout << "Error GENO!!\n";
                                        break;
                                }

				m_OneSNP_GenoTemp[ind] = bufferGeno;
				ind++;
                                geno1 = geno1 >> 1;
                                if(ind >= N){

                                        int indxInOut = 0;
                                //cout << "HERE4\n";
                                        for(int indx=0; indx < Nnomissing; indx++){
                                                //cout << "HERE5\n";
                                                m_OneSNP_Geno[indxInOut] = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];

                                                if(m_OneSNP_Geno[indxInOut] == 3){
                                                        indexNA.push_back(indxInOut);
                                                }
                                                indxInOut = indxInOut + 1;
                                	}

                                        return & m_OneSNP_Geno;
                                }
		

			}		
		}

		 return & m_OneSNP_Geno;

	}



 	int Get_OneSNP_StdGeno(size_t SNPIdx, arma::fvec * out ){
		//avoid large continuous memory usage
                int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
                ////////////////

 		out->zeros(Nnomissing);

 		size_t Start_idx = m_size_of_esi * SNPIdxinVec;
		size_t ind= 0;
		unsigned char geno1;
		
		float freq = alleleFreqVec[SNPIdx];
 		float invStd = invstdvVec[SNPIdx];
		for(size_t i=Start_idx; i< Start_idx+m_size_of_esi; i++){
//			geno1 = genoVec[i];
			geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

			for(int j=0; j<4; j++){
    			int b = geno1 & 1 ;
    			geno1 = geno1 >> 1;
    			int a = geno1 & 1 ;
    			(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
			ind++;
    			geno1 = geno1 >> 1;
    			
    			if(ind >= Nnomissing){
    				return 1;
    			}
    		}
		}
		
		return 1;

		
 	}


	arma::fvec * Get_Diagof_StdGeno(){
	
		arma::fvec * temp = &m_OneSNP_StdGeno;
		// Not yet calculated
		if(size(m_DiagStd)[0] != Nnomissing){
			m_DiagStd.zeros(Nnomissing);
			for(size_t i=0; i< M; i++){
				Get_OneSNP_StdGeno(i, temp);
				m_DiagStd = m_DiagStd + (*temp) % (*temp);
			}
		
		}
	
		return & m_DiagStd;
	}
 	
 
  	//Function to assign values to all attributes
  	//This function is used instead of using a constructor because using constructor can not take
  	//genofile as an argument from runModel.R 
        //genofile is the predix for plink bim, bed, fam, files   
  	void setGenoObj(std::string genofile, std::vector<int> subSampleInGeno, float memoryChunk){
		//cout << "OK1\n";   
		ptrsubSampleInGeno = subSampleInGeno;
		Nnomissing = subSampleInGeno.size(); 
    		// reset
    		//genoVec.clear();
    		//alleleFreqVec.clear();
  		invstdvVec.clear();
   		M=0;
  		N=0;
   	
		std::string bedfile = genofile+".bed";
		std::string bimfile = genofile+".bim"; 
		std::string famfile = genofile+".fam"; 
		std::string junk;
		//cout << "OK2\n";
		//count the number of individuals
		ifstream test_famfile;
		test_famfile.open(famfile.c_str());
        	if (!test_famfile.is_open()){
                	printf("Error! fam file not open!");
                	return ;
        	}
		int indexRow = 0;
		while (std::getline(test_famfile,junk)){
                	indexRow = indexRow + 1;
                	junk.clear();
        	}
		N = indexRow;
		test_famfile.clear();	
		//cout << "OK3\n";
		//count the number of markers
		ifstream test_bimfile;
        	test_bimfile.open(bimfile.c_str());
        	if (!test_bimfile.is_open()){
                	printf("Error! bim file not open!");
                	return ;
        	}
        	indexRow = 0;
        	while (std::getline(test_bimfile,junk)){
                	indexRow = indexRow + 1;
                	junk.clear();
        	}
        	M = indexRow;
        	test_bimfile.clear(); 

    		junk.clear();
		//cout << "OK3b\n";
    		// Init OneSNP Geno
    		Init_OneSNP_Geno();
		//cout << "OK3c\n";
    
    		//std::string junk;
    		indexRow = 0;
    		int buffer;
    		int TotalRead=0;

		std::vector<unsigned char> genoVecOneMarkerOld;
		std::vector<unsigned char> genoVecOneMarkerNew;
		/////////////////////////////
		// Added for reserve for genoVec
		size_t nbyteOld = ceil(float(N)/4);
		size_t nbyteNew = ceil(float(Nnomissing)/4);
		size_t reserve = ceil(float(Nnomissing)/4) * M + M*2;
		cout << "nbyte: " << nbyteOld << endl;
		cout << "nbyte: " << nbyteNew << endl;		
		cout << "reserve: " << reserve << endl;		

    		genoVecOneMarkerOld.reserve(nbyteOld);
    		genoVecOneMarkerOld.resize(nbyteOld);
 		//genoVec.reserve(reserve);
		
		//cout << "OK4\n";

		ifstream test_bedfile;
        	test_bedfile.open(bedfile.c_str(), ios::binary);
        	if (!test_bedfile.is_open()){
                	printf("Error! file open!");
                	return;
        	}		
		//printf("\nM: %zu, N: %zu, Reserve: %zu\n", M, N, reserveTemp);
		printf("\nM: %zu, N: %zu\n", M, N);

        	//test_bedfile.seekg(3);
		//set up the array of vectors for genotype
		numMarkersofEachArray = floor((memoryChunk*pow (10.0, 9.0))/(ceil(float(N)/4)));
		//cout << "numMarkersofEachArray: " << numMarkersofEachArray << endl;
		if(M % numMarkersofEachArray == 0){
                        numofGenoArray = M / numMarkersofEachArray;
			genoVecofPointers.resize(numofGenoArray);
                        //cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
                        for (int i = 0; i < numofGenoArray ; i++){
                                genoVecofPointers[i] = new vector<unsigned char>;
                                genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(N)/4));
                        }

                }else{
                        numofGenoArray = M/numMarkersofEachArray + 1;
                        genoVecofPointers.resize(numofGenoArray);
			numMarkersofLastArray = M - (numofGenoArray-1)*numMarkersofEachArray;
                        cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
			try{	
                        for (int i = 0; i < numofGenoArray-1; i++){
				//cout << "i = " << i << endl;
                                genoVecofPointers[i] = new vector<unsigned char>;
                                genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(N)/4));
				//cout <<((*genoVecofPointers[i]).capacity()==numMarkersofEachArray*ceil(float(N)/4))<< endl;
                        }
			cout << "here\n";
			genoVecofPointers[numofGenoArray-1] = new vector<unsigned char>;
			genoVecofPointers[numofGenoArray-1]->reserve(numMarkersofLastArray*ceil(float(N)/4));
			}
			catch(std::bad_alloc& ba)
                        {
                                std::cerr << "bad_alloc caught1: " << ba.what() << '\n';
                                exit(EXIT_FAILURE);
                        }
			/*
			numMarkersofLastArray = M - (numofGenoArray-1)*numMarkersofEachArray;
			cout << "numMarkersofLastArray " << numMarkersofLastArray << endl;
    			cout << "numMarkersofLastArray*ceil(float(N)/4) " << numMarkersofLastArray*ceil(float(N)/4) << endl;
			genoVecofPointers[numofGenoArray - 1] = new vector<unsigned char>;
			cout << "setgeno mark0" << endl;
			try{
			genoVecofPointers[numofGenoArray - 1]->reserve(numMarkersofLastArray*ceil(float(N)/4));
			}
			 catch(std::bad_alloc& ba)
  			{
    				std::cerr << "bad_alloc caught1: " << ba.what() << '\n';
    				exit(EXIT_FAILURE);
  			}
			*/
			//cout << "setgeno mark0b" << endl;
		}


		alleleFreqVec.zeros(M);
        	float freq, Std, invStd;
        	std::vector<int> indexNA;
        	int lengthIndexNA;
        	int indexGeno;
        	int indexBit;
        	int fillinMissingGeno;
        	int b2;
        	int a2;

		size_t ind= 0;
                unsigned char geno1 = 0;
                int bufferGeno;
                int u;
		//std::vector<int> genoVec4Markers(4);
		//test_bedfile.read((char*)(&genoVecTemp[0]),nbyteTemp*M);


		for(int i = 0; i < M; i++){
			genoVecOneMarkerOld.clear();
			genoVecOneMarkerOld.reserve(nbyteOld);
                        genoVecOneMarkerOld.resize(nbyteOld);

			test_bedfile.seekg(3+nbyteOld*i);
			test_bedfile.read((char*)(&genoVecOneMarkerOld[0]),nbyteOld);
 			//printf("\nFile read is done: M: %zu, N: %zu, TotalByte %zu\n", M, N, genoVecTemp.size());
			//cout << "Imputing missing genotypes and extracting the subset of samples with nonmissing genotypes and phenotypes\n";  
	//		cout << "i is " << i << endl;  

      			indexNA.clear();
		//}	
        		Get_OneSNP_Geno_atBeginning(i, indexNA, genoVecOneMarkerOld);

			ind = 0;
			geno1 = 0;

                	for(unsigned int j=0; j< Nnomissing; j++){
				u = j & (4 - 1);
				bufferGeno = m_OneSNP_Geno[j];
				if(bufferGeno == 0){
					setGenotype(&geno1, u, HOM_ALT);
				}else if(bufferGeno == 1){
					setGenotype(&geno1, u, HET);
				}else if(bufferGeno == 2){	
					setGenotype(&geno1, u, HOM_REF);	
				}else{
					setGenotype(&geno1, u, MISSING);
					m_OneSNP_Geno[j] = 0;  //12-18-2017 	
				}

				if(u == 3){
					//genoVec.push_back(geno1);
					genoVecofPointers[i/numMarkersofEachArray]->push_back(geno1); //avoid large continuous memory usage
					geno1 = 0;
				}
			}
		//cout << "size of genoVec: " << genoVec.size() << endl;
				
			if(Nnomissing%4 != 0){
				//genoVec.push_back(geno1);
				genoVecofPointers[i/numMarkersofEachArray]->push_back(geno1); //avoid large continuous memory usage
				geno1 = 0;
			}

		
			lengthIndexNA = indexNA.size();
      			freq = sum(m_OneSNP_Geno)/(2*(Nnomissing-lengthIndexNA));


			//cout << "setgeno mark3" << endl;
		
			if (lengthIndexNA > 0){

				fillinMissingGeno = int(round(2*freq));
				if(fillinMissingGeno == 0){
                                        b2 = 1;
                                        a2 = 1;
                        	}else if(fillinMissingGeno == 1){
                                        b2 = 0;
                                        a2 = 1;
                        	}else{
                                        b2 = 0;
                                        a2 = 0;
                        	}

		

				for (int k=0; k<lengthIndexNA; k++){
					indexGeno = indexNA[k];
					m_OneSNP_Geno[indexGeno] = fillinMissingGeno;
					//genoVecofPointers[i/numMarkersofEachArray]
					setBit(genoVecofPointers[i/numMarkersofEachArray]->at((i%numMarkersofEachArray)*nbyteNew+(indexGeno/4)),indexGeno%4, a2, b2);
					//setBit(genoVec[i*nbyteNew+(indexGeno/4)], indexGeno%4, a2, b2);
				}


			}


			freq = sum(m_OneSNP_Geno)/(2*Nnomissing);
      			Std = std::sqrt(2*freq*(1-freq));
      			if(Std == 0){
      				invStd= 0;
      			} else {
      				invStd= 1/Std;
      			}
			alleleFreqVec[i] = freq;
      			invstdvVec.push_back(invStd);
			m_OneSNP_Geno.clear();

    		}//end for(int i = 0; i < M; i++){

        	test_bedfile.close();
		printAlleleFreqVec();
		printGenoVec();
   		Get_Diagof_StdGeno();
  	}//End Function
 

  	void printFromgenoVec(unsigned char genoBinary0){
		unsigned char genoBinary = genoBinary0;
  		for(int j=0; j<4; j++){
        		int b = genoBinary & 1 ;
                	genoBinary = genoBinary >> 1;
                	int a = genoBinary & 1 ;
			genoBinary = genoBinary >> 1;
			cout << 2-(a+b) << " " << endl;
		}
		cout << endl;
  	}
 
  
  	int getM() const{
    		return(M);
  	}
  
  	int getN() const{
    		return(N);
  	}
 
  	int getNnomissing() const{
    		return(Nnomissing);
  	}


  	float getAC(int m){
    		return(alleleFreqVec[m]*2*Nnomissing);
  	}

  	float getMAC(int m){
    		if(alleleFreqVec[m] > 0.5){
      			return((1-alleleFreqVec[m])*2*Nnomissing);
    		}else{
      			return(alleleFreqVec[m]*2*Nnomissing);
    		}
  	}


 
  	//print out the vector of genotypes
  	void printGenoVec(){
    		//for(unsigned int i=0; i<M; ++i)
    		for(unsigned int i=0; i<2; ++i)
    		{
	
    			Get_OneSNP_Geno(i);
    			//for(unsigned int j=0; j< Nnomissing; j++){
    			for(unsigned int j=0; j< 100; j++){
      				cout << m_OneSNP_Geno[j] << ' ';
      			}
      			cout << endl;
    		}
    		//cout << "genoVec.size()" << genoVec.size() << endl;
    		cout << "M = " << M << endl;
    		cout << "N = " << N << endl;
  	}
  
  	//print out the vector of allele frequency
  	void printAlleleFreqVec(){
    		//for(int i=0; i<alleleFreqVec.size(); ++i)
    		for(int i=(M-100); i<M; ++i)
    		{
      			cout << alleleFreqVec[i] << ' ';
    		}
    		cout << endl;
  	}

  
};



// //create a geno object as a global variable
genoClass geno;


// [[Rcpp::export]]
void closeGenoFile_plink()
{
  //genoToTest_plainDosage.test_genoGZfile.close();
	for (int i = 0; i < geno.numofGenoArray; i++){
		(*geno.genoVecofPointers[i]).clear();	
    		delete geno.genoVecofPointers[i];
  	}

  	geno.genoVecofPointers.clear();

  	//geno.genoVec.clear();
  	geno.invstdvVec.clear();
  	geno.ptrsubSampleInGeno.clear();
  	geno.alleleFreqVec.clear();
  	geno.m_OneSNP_Geno.clear();
  	geno.m_OneSNP_StdGeno.clear();
  	geno.m_DiagStd.clear();
  	printf("closed the plinkFile!\n");
}


// [[Rcpp::export]]
int gettotalMarker (){
  	int numMarker = geno.getM();
  	return(numMarker); 
}

// [[Rcpp::export]]
arma::fvec getAlleleFreqVec(){
  	return(geno.alleleFreqVec);
}




//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd : public Worker
{   
  	// source vectors
	arma::fcolvec & m_bVec;
	unsigned int m_N;
	unsigned int m_M;
  
  	// product that I have accumulated
  	arma::fvec m_bout;
  	
  
  	// constructors
  	CorssProd(arma::fcolvec & y)
  		: m_bVec(y) {
  		
  		m_M = geno.getM();
  		m_N = geno.getNnomissing();
  		m_bout.zeros(m_N);
  	} 
  	CorssProd(const CorssProd& CorssProd, Split)
  		: m_bVec(CorssProd.m_bVec)
  	{

  		m_N = CorssProd.m_N;
  		m_M = CorssProd.m_M;
  		m_bout.zeros(m_N);
  	
  	}  
  	// process just the elements of the range I've been asked to
  	void operator()(std::size_t begin, std::size_t end) {
  	  	arma::fvec vec;
  	  	for(unsigned int i = begin; i < end; i++){
			
			geno.Get_OneSNP_StdGeno(i, &vec);
			float val1 = dot(vec,  m_bVec);
			m_bout += val1 * (vec);
  		}
  	}
  
  	// join my value with that of another InnerProduct
  	void join(const CorssProd & rhs) { 
    	m_bout += rhs.m_bout; 
  	}
};

// [[Rcpp::export]]
arma::fvec parallelCrossProd(arma::fcolvec & bVec) {
  
  // declare the InnerProduct instance that takes a pointer to the vector data
  	int M = geno.getM();
  
  	CorssProd CorssProd(bVec);
  
  // call paralleReduce to start the work
  	parallelReduce(0, M, CorssProd);
  
  // return the computed product
  	return CorssProd.m_bout/M;
}


// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin(arma::fcolvec& bVec){
  
  	arma::fvec crossProdVec = parallelCrossProd(bVec) ;
  
  	return(crossProdVec);
}


// [[Rcpp::export]]
void setgeno(std::string genofile, std::vector<int> & subSampleInGeno, float memoryChunk)
{
	
	int start_s=clock();
        geno.setGenoObj(genofile, subSampleInGeno, memoryChunk);
	//geno.printAlleleFreqVec();
	//geno.printGenoVec();
	int stop_s=clock();
	cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
}




// [[Rcpp::export]]
arma::fvec Get_OneSNP_Geno(int SNPIdx)
{

	arma::fvec temp = * geno.Get_OneSNP_Geno(SNPIdx);
	return(temp);

}


// [[Rcpp::export]]
arma::fvec Get_OneSNP_Geno_forVarianceRatio(int SNPIdx)
{
       
        arma::fvec temp = * geno.Get_OneSNP_Geno(SNPIdx);
        return(temp);

}



// [[Rcpp::export]]
arma::fvec Get_OneSNP_StdGeno(int SNPIdx)
{

	arma::fvec temp; 
	geno.Get_OneSNP_StdGeno(SNPIdx, & temp);
	return(temp);

}
  
    
  

//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
// [[Rcpp::export]]
arma::fvec getDiagOfSigma(arma::fvec& wVec, arma::fvec& tauVec){
  
	int Nnomissing = geno.getNnomissing();
	int M = geno.getM();
	
	//cout << "N=" << N << endl;
	arma::fvec diagVec(Nnomissing);
	float diagElement;
	float floatBuffer;
  	//float minvElement;
  
   
	diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /M + tauVec(0)/wVec;
	for(unsigned int i=0; i< Nnomissing; i++){
  		if(diagVec(i) < 1e-4){
  			diagVec(i) = 1e-4 ;
  		}
  	}
  
    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
  	return(diagVec);

}

// [[Rcpp::export]]
arma::fcolvec getCrossprod(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;

        return(crossProdVec);
}




//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
//This function needs the function getDiagOfSigma and function getCrossprod
// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG){
  	arma::fvec rVec = bVec;
  	arma::fvec r1Vec;
  	int Nnomissing = geno.getNnomissing();

  	arma::fvec crossProdVec(Nnomissing);
  	arma::fvec minvVec = 1/getDiagOfSigma(wVec, tauVec);
  	float sumr2 = sum(rVec % rVec);

  	arma::fvec zVec = minvVec % rVec;
  	arma::fvec z1Vec;
 	arma::fvec pVec = zVec;
  
  	arma::fvec xVec(Nnomissing);
  	xVec.zeros();
  
  	int iter = 0;
  	while (sumr2 > tolPCG && iter < maxiterPCG) {
    		iter = iter + 1;
    		arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
    		arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

    		float a = preA(0);
    
    		xVec = xVec + a * pVec;
    
    		r1Vec = rVec - a * ApVec;

    		z1Vec = minvVec % r1Vec;
    		arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
    		float bet = Prebet(0);
    		pVec = z1Vec+ bet*pVec;
    		zVec = z1Vec;
    		rVec = r1Vec;
    
    		sumr2 = sum(rVec % rVec);
  	}
  
  	if (iter >= maxiterPCG){
    		cout << "pcg did not converge. You may increase maxiter number." << endl;
     
  	}
  	cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
  	return(xVec);
}


//http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
	Rcpp::Environment base_env("package:base");
  	Rcpp::Function set_seed_r = base_env["set.seed"];
  	set_seed_r(seed);  
}

// [[Rcpp::export]]
Rcpp::NumericVector nb(int n) {
  	return(rbinom(n,1,0.5));
}


//This function needs the function getPCG1ofSigmaAndVector and function getCrossprodMatAndKin
// [[Rcpp::export]]
float GetTrace(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& cov1, int nrun, int maxiterPCG, float tolPCG){
  
  	set_seed(200);
  	int Nnomissing = geno.getNnomissing();
  	arma::fmat Sigma_iXt = Sigma_iX.t();  

  	arma::fvec tempVec(nrun);
  	tempVec.zeros();
  	for(int i = 0; i < nrun; i++){
    		Rcpp::NumericVector uVec0;
    		uVec0 = nb(Nnomissing);
    		arma::fvec uVec = as<arma::fvec>(uVec0);
    		uVec = uVec*2 - 1;
    		arma::fvec Sigma_iu;
    		Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG);
    
    
    		arma::fcolvec Pu;
    		Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
    		arma::fvec Au;

    		Au = getCrossprodMatAndKin(uVec);
    		tempVec(i) = dot(Au, Pu);
  	}
  	float tra = sum(tempVec)/nrun;
  	return(tra);
}	


// Added by SLEE, 04/16/2017
//      This function calculate fixed and random effect coefficients
// [[Rcpp::export]]
Rcpp::List getCoefficients(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG){

  	int Nnomissing = geno.getNnomissing();
  	arma::fvec Sigma_iY;

  	Sigma_iY = getPCG1ofSigmaAndVector(wVec, tauVec, Yvec, maxiterPCG, tolPCG);
  	int colNumX = Xmat.n_cols;
  	arma::fmat Sigma_iX(Nnomissing,colNumX);
  	arma::fvec XmatVecTemp;
  	for(int i = 0; i < colNumX; i++){
    		XmatVecTemp = Xmat.col(i);

    		Sigma_iX.col(i) = getPCG1ofSigmaAndVector(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

  	}

  	arma::fmat Xmatt = Xmat.t();
  	arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
 	arma::fmat Sigma_iXt = Sigma_iX.t();
  	arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
  	arma::fvec alpha = cov * SigmaiXtY;

  	arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
  	return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
}


// Modified by SLEE, 04/16/2017
// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//      This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace
// [[Rcpp::export]]
Rcpp::List getAIScore(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG){

	arma::fmat Sigma_iXt = Sigma_iX.t();

  	arma::fvec PY1 = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));
  	arma::fvec APY = getCrossprodMatAndKin(PY1);
  	float YPAPY = dot(PY1, APY);

  	float Trace = GetTrace(Sigma_iX, Xmat, wVec, tauVec, cov, nrun, maxiterPCG, tolPCG);
  	arma::fvec PAPY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG, tolPCG);
  	arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
  	float AI = dot(APY, PAPY);

  	return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY1, Named("AI") = AI);
}

// Modified by SLEE, 04/16/2017
// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore
// [[Rcpp::export]]
Rcpp::List fitglmmaiRPCG(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float tol){

  	Rcpp::List re = getAIScore(Yvec, Xmat,wVec,  tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG);
  	float YPAPY = re["YPAPY"];
  	float Trace = re["Trace"];
  	float score1 = YPAPY - Trace;
  	float AI1 = re["AI"];
  	float Dtau = score1/AI1;
  	arma::fvec tau0 = tauVec;
  	tauVec(1) = tau0(1) + Dtau;

  	for(int i=0; i<tauVec.n_elem; ++i) {
    		if (tauVec(i) < tol){
      			tauVec(i) = 0;
    		}
  	}

  	float step = 1.0;
  	while (tauVec(1) < 0.0){

    		step = step*0.5;
    		tauVec(1) = tau0(1) + step * Dtau;

  	}

  //cout << "DEBUG4" << endl;
  	for(int i=0; i<tauVec.n_elem; ++i) {
    		if (tauVec(i) < tol){
      			tauVec(i) = 0;
    		}
  	}
  	return List::create(Named("tau") = tauVec);
}



/*add for SPA by Wei 04222017*/
// [[Rcpp::export]]
arma::fmat getSigma_X(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG){


  	int Nnomissing = Xmat.n_rows;
  	int colNumX = Xmat.n_cols;

  	cout << colNumX << endl;
  	cout << size(wVec) << endl;
  	cout << size(tauVec) << endl;


  	arma::fmat Sigma_iX1(Nnomissing,colNumX);
  	arma::fvec XmatVecTemp;

  	for(int i = 0; i < colNumX; i++){
    		XmatVecTemp = Xmat.col(i);
    		Sigma_iX1.col(i) = getPCG1ofSigmaAndVector(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);
  	}
  	return(Sigma_iX1);
}


// [[Rcpp::export]]
arma::fvec  getSigma_G(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG){
  	arma::fvec Sigma_iG;
  	Sigma_iG = getPCG1ofSigmaAndVector(wVec, tauVec, Gvec, maxiterPCG, tolPCG);
  	return(Sigma_iG);
}



//This function needs the function getPCG1ofSigmaAndVector and function getCrossprodMatAndKin
// [[Rcpp::export]]
arma::fvec GetTrace_q(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& cov1,  int nrun, int maxiterPCG, float tolPCG){

  	set_seed(200);

  	arma::fmat Sigma_iXt = Sigma_iX.t();
  	int Nnomissing = geno.getNnomissing();
  	arma::fvec tempVec(nrun);
  	tempVec.zeros();
  	arma::fvec tempVec0(nrun);
  	tempVec0.zeros();
  	for(int i = 0; i < nrun; i++){

    		Rcpp::NumericVector uVec0;
    		uVec0 = nb(Nnomissing);
    		arma::fvec uVec = as<arma::fvec>(uVec0);
    		uVec = uVec*2 - 1;
    		arma::fvec Sigma_iu;
    		Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG);


    		arma::fcolvec Pu;

    		Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
    		arma::fvec Au;
    		Au = getCrossprodMatAndKin(uVec);
    		tempVec(i) = dot(Au, Pu);
    		tempVec0(i) = dot(uVec, Pu);

  	}
  	arma::fvec traVec(2);
  	traVec(1) = sum(tempVec)/nrun;
  	traVec(0) = sum(tempVec0)/nrun;
  	return(traVec);
}

//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace
// [[Rcpp::export]]
Rcpp::List getAIScore_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG){

  	int Nnomissing = geno.getNnomissing();
  	arma::fvec Sigma_iY1;
  	Sigma_iY1 = getPCG1ofSigmaAndVector(wVec, tauVec, Yvec, maxiterPCG, tolPCG);

  	int colNumX = Xmat.n_cols;
  	arma::fmat Sigma_iX1(Nnomissing,colNumX);
  	arma::fvec XmatVecTemp;

 	for(int i = 0; i < colNumX; i++){
    		XmatVecTemp = Xmat.col(i);

    		Sigma_iX1.col(i) = getPCG1ofSigmaAndVector(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

  	}


  	arma::fmat Sigma_iX1t = Sigma_iX1.t();
  	arma::fmat Xmatt = Xmat.t();

  	arma::fmat cov1 = inv_sympd(Xmatt * Sigma_iX1);
  	arma::fvec PY1 = Sigma_iY1 - Sigma_iX1 * (cov1 * (Sigma_iX1t * Yvec));
  	arma::fvec APY = getCrossprodMatAndKin(PY1);

  	float YPAPY = dot(PY1, APY);

  	arma::fvec A0PY = PY1; ////Quantitative


  	float YPA0PY = dot(PY1, A0PY); ////Quantitative

  	arma::fvec Trace = GetTrace_q(Sigma_iX1, Xmat, wVec, tauVec, cov1, nrun, maxiterPCG, tolPCG);

  	arma::fmat AI(2,2);
  	arma::fvec PA0PY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, A0PY, maxiterPCG, tolPCG);
  	arma::fvec PA0PY = PA0PY_1 - Sigma_iX1 * (cov1 * (Sigma_iX1t * PA0PY_1));

  	AI(0,0) =  dot(A0PY, PA0PY);

  	cout << "A1(0,0) " << AI(0,0)  << endl;
  	arma::fvec PAPY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG, tolPCG);
  	arma::fvec PAPY = PAPY_1 - Sigma_iX1 * (cov1 * (Sigma_iX1t * PAPY_1));
  	AI(1,1) = dot(APY, PAPY);

  	AI(0,1) = dot(A0PY, PAPY);

  	AI(1,0) = AI(0,1);

  	cout << "AI " << AI << endl;
  	cout << "Trace " << Trace << endl;
  	cout << "YPAPY " << YPAPY << endl;
  	cout << "cov " << cov1 << endl;

  	return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace,Named("Sigma_iY") = Sigma_iY1, Named("Sigma_iX") = Sigma_iX1, Named("PY") = PY1, Named("AI") = AI, Named("cov") = cov1,  Named("YPA0PY") = YPA0PY);
}




//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore_q
// [[Rcpp::export]]
Rcpp::List fitglmmaiRPCG_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float tol){

  	arma::uvec zeroVec = (tauVec < tol); //for Quantitative, GMMAT
  	Rcpp::List re = getAIScore_q(Yvec, Xmat, wVec, tauVec, nrun, maxiterPCG, tolPCG);

  	arma::fmat cov = re["cov"];
  	arma::fmat Sigma_iX = re["Sigma_iX"];
 	arma::fmat Sigma_iXt = Sigma_iX.t();

  	arma::fvec alpha1 = cov * (Sigma_iXt * Yvec);
  	arma::fvec Sigma_iY = re["Sigma_iY"];
  	arma::fvec eta1 = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha1) / wVec;
  	float YPAPY = re["YPAPY"];
  	float YPA0PY = re["YPA0PY"]; //for Quantitative
  	arma::fvec Trace = re["Trace"]; //for Quantitative

  	float score0 = YPA0PY - Trace(0); //for Quantitative
  	float score1 = YPAPY - Trace(1); //for Quantitative
  	arma::fvec scoreVec(2); //for Quantitative
  	scoreVec(0) = score0; //for Quantitative
  	scoreVec(1) = score1; //for Quantitative

    //for Quantitative
  	arma::fmat AI = re["AI"];
  	cout << "0,0" << AI(0,0) << endl;
  	cout << "0,1" << AI(0,1) << endl;
  	cout << "1,0" << AI(1,0) << endl;
  	cout << "1,1" << AI(1,1) << endl;

  	arma::fvec Dtau = solve(AI, scoreVec);


  	arma::fvec tau0 = tauVec;
  	tauVec = tau0 + Dtau;


  	tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros(); //for Quantitative Copied from GMMAT  

  	float step = 1.0;


  	cout << "tau2 " << tauVec(0) << " " << tauVec(1) << endl;
  	while (tauVec(0) < 0.0 || tauVec(1)  < 0.0){ //for Quantitative
     		cout << "tauVec Here: " << tauVec << endl;
    		step = step*0.5;
    		tauVec = tau0 + step * Dtau; //for Quantitative
    		cout << "tau_4: " << tauVec << endl;
    		tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros(); //for Quantitative Copied from GMMAT
    		cout << "tau_5: " << tauVec << endl;
 	}



  	tauVec.elem( find(tauVec < tol) ).zeros();
  	return List::create(Named("tau") = tauVec, Named("cov") = cov, Named("alpha") = alpha1, Named("eta") = eta1);
}

