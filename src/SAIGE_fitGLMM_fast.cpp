#define ARMA_USE_SUPERLU 1
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h> 
#include <omp.h>
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


float minMAFtoConstructGRM = 0;
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
  	arma::fvec	invstdvVec;
	vector<int>	ptrsubSampleInGeno;
	
  	arma::fvec 	alleleFreqVec;
  	arma::ivec	m_OneSNP_Geno;
  	arma::fvec	m_OneSNP_StdGeno;
  	arma::fvec	m_DiagStd;

	arma::ivec	MACVec; //for variance ratio based on different MAC categories
	arma::ivec	subMarkerIndex; //for sparse GRM
	arma::fmat      stdGenoMultiMarkersMat;	
	std::vector<float> stdGenoforSamples; //for sparse GRM
	std::vector<float>     kinValueVecFinal;
        float relatednessCutoff;

	tbb::concurrent_vector< std::pair<int, int> > indiceVec;
	arma::ivec xout;
        arma::ivec yout;
	//int Mmafge1perc;
	bool setKinDiagtoOne;
	int numberofMarkerswithMAFge_minMAFtoConstructGRM = 0;
//	arma::SpMat<float> sparseGRMinC(2,2);
	


        //std::vector<float> stdGenoVec;
	//for LOCO
	//bool LOCO = false;
	//vector<int> chromosomeStartIndex;
	//vector<int> chromosomeEndIndex;
	//vector<int> chromosomeVec;
        size_t Msub;
        int startIndex;
        int endIndex;
	int Msub_MAFge_minMAFtoConstructGRM;
	//end for LOCO

	unsigned char m_genotype_buffer[4];
	int geno_idx;
	int m_size_of_esi;
	unsigned char m_bits_val[8];

	
	//look-up table for std geno
	//float stdGenoLookUpArr[3] = {0};
	void setStdGenoLookUpArr(float mafVal, float invsdVal, arma::fvec & stdGenoLookUpArr){
	//	arma::fvec stdGenoLookUpArr(3);
		float mafVal2 = 2*mafVal;
		stdGenoLookUpArr(0) = (0-mafVal2)*invsdVal;
		stdGenoLookUpArr(1) = (1-mafVal2)*invsdVal;
		stdGenoLookUpArr(2) = (2-mafVal2)*invsdVal;
	//	return(stdGenoLookUpArr)
	}


        //look-up table in a 2D array for sparseKin 
        float sKinLookUpArr[3][3] = {{0}};
	//(g - 2*freq)* invStd;;
        void setSparseKinLookUpArr(float mafVal, float invsdVal){
		float mafVal2 = 2*mafVal;
		float a0 = (0-mafVal2)*invsdVal;
		float a1 = (1-mafVal2)*invsdVal;
		float a2 = (2-mafVal2)*invsdVal;
		
		sKinLookUpArr[0][0] = a0*a0;
		sKinLookUpArr[0][1] = a0*a1;
		sKinLookUpArr[0][2] = a0*a2;
		sKinLookUpArr[1][0] = sKinLookUpArr[0][1];
		sKinLookUpArr[1][1] = a1*a1;
		sKinLookUpArr[1][2] = a1*a2;
		sKinLookUpArr[2][0] = sKinLookUpArr[0][2];
		sKinLookUpArr[2][1] = sKinLookUpArr[1][2];
		sKinLookUpArr[2][2] = a2*a2;

	}




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
	

        arma::ivec * Get_OneSNP_Geno(size_t SNPIdx){
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
   
	arma::ivec * Get_OneSNP_Geno_atBeginning(size_t SNPIdx, vector<int> & indexNA, vector<unsigned char> & genoVecOneMarkerOld){

		arma::ivec m_OneSNP_GenoTemp;
		m_OneSNP_GenoTemp.zeros(N);
		m_OneSNP_Geno.zeros(Nnomissing);
		int m_size_of_esi_temp = (N+3)/4;
		size_t ind= 0;
		unsigned char geno1;
		int bufferGeno;
		for(int i=0; i< m_size_of_esi_temp; i++){
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
//		cout << "Get_OneSNP_StdGeno here" << endl; 
		float invStd = invstdvVec[SNPIdx];

		arma::fvec stdGenoLookUpArr(3);
		setStdGenoLookUpArr(freq, invStd, stdGenoLookUpArr);

		//setStdGenoLookUpArr(freq, invStd);
		//std::cout << "stdGenoLookUpArr[0]: " << stdGenoLookUpArr[0] << std::endl;
		//std::cout << "stdGenoLookUpArr[1]: " << stdGenoLookUpArr[1] << std::endl;
		//std::cout << "stdGenoLookUpArr[2]: " << stdGenoLookUpArr[2] << std::endl;
//		cout << "Get_OneSNP_StdGeno here2"  << endl;
		for(size_t i=Start_idx; i< Start_idx+m_size_of_esi; i++){
//			geno1 = genoVec[i];
			geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

			for(int j=0; j<4; j++){
    			int b = geno1 & 1 ;
    			geno1 = geno1 >> 1;
    			int a = geno1 & 1 ;
    			//(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
    			(*out)[ind] = stdGenoLookUpArr(2-(a+b));
			ind++;
    			geno1 = geno1 >> 1;
    			
    			if(ind >= Nnomissing){
//				cout << "Get_OneSNP_StdGeno " << SNPIdx << endl; 
//				cout << "Nnomissing " << Nnomissing << endl; 
				stdGenoLookUpArr.clear();
    				return 1;
    			}
    		}
		}
		stdGenoLookUpArr.clear();
		return 1;
		
		
 	}


	arma::fvec * Get_Diagof_StdGeno(){
	
		arma::fvec * temp = &m_OneSNP_StdGeno;
		// Not yet calculated
		//cout << "size(m_DiagStd)[0] " << size(m_DiagStd)[0] << endl;
		if(size(m_DiagStd)[0] != Nnomissing){
			m_DiagStd.zeros(Nnomissing);
			for(size_t i=0; i< M; i++){
				if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
				Get_OneSNP_StdGeno(i, temp);

				/*if(i == 0){
					cout << "setgeno mark7 " << i <<  endl;
					for(int j=0; j<10; ++j)
					{
                				cout << (*temp)[j] << ' ';
                			}
                			cout << endl;
				}
				*/
				m_DiagStd = m_DiagStd + (*temp) % (*temp);
				}
			}
		
		}

		//std::cout << "test\n";
		//for(int i=0; i<10; ++i)
        	//{
        	//  cout << m_DiagStd[i] << ' ';
        	//}
		//cout << endl;
	
		return & m_DiagStd;
	}

	
	

 	

	arma::fvec * Get_Diagof_StdGeno_LOCO(){

                arma::fvec * temp = &m_OneSNP_StdGeno;
		Msub_MAFge_minMAFtoConstructGRM = 0;	
                // Not yet calculated
                if(size(m_DiagStd)[0] != Nnomissing){
                        m_DiagStd.zeros(Nnomissing);
                        for(size_t i=0; i< M; i++){
				if(i < startIndex || i > endIndex){
					if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
                                  Get_OneSNP_StdGeno(i, temp);
                                  m_DiagStd = m_DiagStd + (*temp) % (*temp);
					Msub_MAFge_minMAFtoConstructGRM = Msub_MAFge_minMAFtoConstructGRM + 1;
					}
				}
                        }

                }

                return & m_DiagStd;
        }


 
  	//Function to assign values to all attributes
  	//This function is used instead of using a constructor because using constructor can not take
  	//genofile as an argument from runModel.R 
        //genofile is the predix for plink bim, bed, fam, files   
  	void setGenoObj(std::string genofile, std::vector<int> subSampleInGeno, float memoryChunk, bool  isDiagofKinSetAsOne){
		//cout << "OK1\n";
		setKinDiagtoOne = isDiagofKinSetAsOne;   
		ptrsubSampleInGeno = subSampleInGeno;
		Nnomissing = subSampleInGeno.size(); 
    		// reset
    		//genoVec.clear();
    		alleleFreqVec.clear();
		MACVec.clear();
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
			//	cout << "i = " << i << endl;
                                genoVecofPointers[i] = new vector<unsigned char>;
                                genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(N)/4));
				//cout <<((*genoVecofPointers[i]).capacity()==numMarkersofEachArray*ceil(float(N)/4))<< endl;
                        }
			//cout << "here\n";
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

		cout << "setgeno mark1" << endl;
		alleleFreqVec.zeros(M);
		invstdvVec.zeros(M);
		MACVec.zeros(M);
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

		cout << "setgeno mark2" << endl;
		//Mmafge1perc = 0;
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
				
			if(Nnomissing%4 != 0){
				//genoVec.push_back(geno1);
				genoVecofPointers[i/numMarkersofEachArray]->push_back(geno1); //avoid large continuous memory usage
				geno1 = 0;
			}

		
			lengthIndexNA = indexNA.size();
      			freq = float(sum(m_OneSNP_Geno))/(2*(Nnomissing-lengthIndexNA));


//			cout << "setgeno mark3" << endl;
		
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
			
			//cout << "setgeno mark4" << endl;

			freq = float(sum(m_OneSNP_Geno))/(2*Nnomissing);
      			Std = std::sqrt(2*freq*(1-freq));
      			if(Std == 0){
      				invStd= 0;
      			} else {
      				invStd= 1/Std;
      			}
			alleleFreqVec[i] = freq;
			//if(freq >= 0.01 && freq <= 0.99){
			//	Mmafge1perc = Mmafge1perc + 1;
			//}
			if(minMAFtoConstructGRM > 0){
				if(freq >= minMAFtoConstructGRM && freq <= (1-minMAFtoConstructGRM)){
					numberofMarkerswithMAFge_minMAFtoConstructGRM = numberofMarkerswithMAFge_minMAFtoConstructGRM + 1;
				}
			}else{
				numberofMarkerswithMAFge_minMAFtoConstructGRM = M;
			}
			
			if(freq > 0.5){
			  MACVec[i] = (2*Nnomissing) - sum(m_OneSNP_Geno);
			}else{
			  MACVec[i] = sum(m_OneSNP_Geno);
			}


      			invstdvVec[i] = invStd;
			m_OneSNP_Geno.clear();

    		}//end for(int i = 0; i < M; i++){

		if( minMAFtoConstructGRM > 0){
			cout << numberofMarkerswithMAFge_minMAFtoConstructGRM << " markers with MAF >= " << minMAFtoConstructGRM << " are used for GRM." << endl;
		}else{
			cout << M << " markers with MAF >= " << minMAFtoConstructGRM << " are used for GRM." << endl;

		}
        	test_bedfile.close();
		cout << "setgeno mark5" << endl;
//		printAlleleFreqVec();
		//printGenoVec();
   		//Get_Diagof_StdGeno();
		cout << "setgeno mark6" << endl;
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

        int getnumberofMarkerswithMAFge_minMAFtoConstructGRM() const{
                return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
        }

	int getMsub_MAFge_minMAFtoConstructGRM_in() const{
		return(Msub_MAFge_minMAFtoConstructGRM);
	}
 
	//int getMmafge1perc() const{
	//	return(Mmafge1perc);
 	//}

	int getMsub() const{
                return(Msub);
        }

	int getStartIndex() const{
		return(startIndex);
	}

	int getEndIndex() const{
                return(endIndex);
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

	//int getnumberofMarkerswithMAFge_minMAFtoConstructGRM(){
 	//	return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
	//}
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


	void Get_Samples_StdGeno(arma::ivec SampleIdsVec){
        	int indexOfVectorPointer;
        	int SNPIdxinVec;

        	int numSamples = SampleIdsVec.n_elem;
        	//stdGenoVec.zeros(Nnomissing*numSamples);
        	stdGenoforSamples.clear();
        	stdGenoforSamples.resize(M*numSamples);

        	arma::ivec sampleGenoIdxVec;
        	sampleGenoIdxVec.zeros(numSamples);
        	arma::ivec sampleGenoIdxSubVec;
        	sampleGenoIdxSubVec.zeros(numSamples);

        	for(int j=0; j < numSamples; j++){
                	sampleGenoIdxVec[j] = SampleIdsVec[j] / 4;
                	sampleGenoIdxSubVec[j] = SampleIdsVec[j] % 4;
        	}


        	int startidx;
        	unsigned char geno1;

        	for(int i=0; i < M; i++){
                	indexOfVectorPointer = i/numMarkersofEachArray;
                	SNPIdxinVec = i % numMarkersofEachArray;
                	startidx = m_size_of_esi * SNPIdxinVec;

                	float freq = alleleFreqVec[i];
                	float invStd = invstdvVec[i];

                	for(int j=0; j < numSamples; j++){
                        	int k = startidx + sampleGenoIdxVec[j];
                        	geno1 = genoVecofPointers[indexOfVectorPointer]->at(k);
                        	for(int q=0; q<4; q++){
                                	if(q == sampleGenoIdxSubVec[j]){
                                        	int b = geno1 & 1 ;
                                        	geno1 = geno1 >> 1;
                                        	int a = geno1 & 1 ;
                                        	stdGenoforSamples[i*(numSamples)+j] = ((2-(a+b)) - 2*freq)* invStd;
                                //(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
                                //ind++;
                                        	geno1 = geno1 >> 1;
                                	}else{
                                        	geno1 = geno1 >> 1;
                                        	geno1 = geno1 >> 1;
                                	}
                        	}
                	}
        	}

        //return(stdGenoVec);
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
int gettotalMarker(){
  	int numMarker = geno.getM();
  	return(numMarker); 
}

// [[Rcpp::export]]
arma::fvec getAlleleFreqVec(){
  	return(geno.alleleFreqVec);
}

// [[Rcpp::export]]
arma::ivec getMACVec(){
        return(geno.MACVec);
}

// [[Rcpp::export]]
arma::ivec getSubMarkerIndex(){
	return(geno.subMarkerIndex);
}

// [[Rcpp::export]]
int getSubMarkerNum(){
        return(geno.subMarkerIndex.n_elem);
}


void initKinValueVecFinal(int ni){
	geno.kinValueVecFinal.resize(ni);
        std::fill(geno.kinValueVecFinal.begin(), geno.kinValueVecFinal.end(), 0);
};

// [[Rcpp::export]]
int getNnomissingOut(){
	return(geno.getNnomissing());
}

//  // [[Rcpp::export]]
//int getMmafge1perc(){
//	return(geno.getMmafge1perc());
//}
//arma::fmat Get_MultiMarkersBySample_StdGeno_Mat(arma::fvec& markerIndexVec){
//arma::fmat Get_MultiMarkersBySample_StdGeno_Mat(){

// [[Rcpp::export]]
void Get_MultiMarkersBySample_StdGeno_Mat(){
	//geno.subMarkerIndex
	//int m_M_Submarker = markerIndexVec.n_elem;
	int m_M_Submarker = getSubMarkerNum();
        //arma::fvec stdGenoMultiMarkers;
        int Nnomissing = geno.getNnomissing();
	  //int nSubMarker = markerIndexVec.n_elem;
          //int Ntotal = geno.getNnomissing();
        //std::vector<float> stdGenoMultiMarkers;
        //stdGenoMultiMarkers.resize(Nnomissing*m_M_Submarker);

        int indexOfVectorPointer;
        int SNPIdxinVec;
        size_t Start_idx;
        size_t ind= 0;
        size_t indtotal = 0;
        unsigned char geno1;
        float freq;
        float invStd;
        int flag;
        int SNPIdx;

//      std::cout << "createSparseKin1d" << std::endl;
        for(size_t k=0; k< m_M_Submarker; k++){
                ind = 0;
                flag = 0;
                //SNPIdx = markerIndexVec[k];
		SNPIdx = (geno.subMarkerIndex)[k];
                indexOfVectorPointer = SNPIdx/(geno.numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (geno.numMarkersofEachArray);
                Start_idx = (geno.m_size_of_esi) * SNPIdxinVec;
                freq = (geno.alleleFreqVec)[SNPIdx];
                invStd = (geno.invstdvVec)[SNPIdx];
                if(k == 0){
                        std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
                }

                while(flag == 0){
//              std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(geno.m_size_of_esi); i++){
                        geno1 = (geno.genoVecofPointers)[indexOfVectorPointer]->at(i);
                        //std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
			(geno.stdGenoMultiMarkersMat)(k, ind) = ((2-(a+b)) - 2*freq)* invStd;
//			std::cout << "k,ind " << k << " " << ind << std::endl;
//			std::cout << "(geno.stdGenoMultiMarkersMat)(k, ind) " << (geno.stdGenoMultiMarkersMat)(k, ind) << std::endl;

//                        stdGenoMultiMarkers[ind*m_M_Submarker+k] = ((2-(a+b)) - 2*freq)* invStd;;
//                      if(k == 0){
    //                    std::cout << "ind*m_M_Submarker+k: " << ind*m_M_Submarker+k << " stdGenoMultiMarkers[ind*m_M_Submarker+k]: " << stdGenoMultiMarkers[ind*m_M_Submarker+k] <<  std::endl;
  //              }


                        indtotal++;
                        ind++;
                        geno1 = geno1 >> 1;

                                if(ind == Nnomissing){
                                        flag = 1;
                                        break;
                                }
                        }// end of for(int j=0; j<4; j++){
                    }// end of for(size_t i=Start_idx
                } //end of while(flag == 0){

        }

        std::cout << "stdGenoMultiMarkersMat.n_rows: " << geno.stdGenoMultiMarkersMat.n_rows << std::endl;
        std::cout << "stdGenoMultiMarkersMat.n_cols: " << geno.stdGenoMultiMarkersMat.n_cols << std::endl;
//	arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), m_M_Submarker, Nnomissing);

//	return(stdGenoMultiMarkersMat);
        //std::cout << "stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] " << stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] << std::endl;

}





// [[Rcpp::export]]
void Get_MultiMarkersBySample_StdGeno(arma::fvec& markerIndexVec, std::vector<float> &stdGenoMultiMarkers){

//	std::cout << "createSparseKin1c" << std::endl;
        int indexOfVectorPointer;
        int SNPIdxinVec;
        size_t Start_idx;
        size_t ind= 0;
        size_t indtotal = 0;
        unsigned char geno1;
        float freq;
        float invStd;
	int flag;
	int SNPIdx;

        int m_M_Submarker = markerIndexVec.n_elem;
        //arma::fvec stdGenoMultiMarkers;
	int Nnomissing = geno.getNnomissing();
	

//	std::cout << "createSparseKin1d" << std::endl;
        for(size_t k=0; k< m_M_Submarker; k++){
                ind = 0;
                flag = 0;
                SNPIdx = markerIndexVec[k];
                indexOfVectorPointer = SNPIdx/(geno.numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (geno.numMarkersofEachArray);
                Start_idx = (geno.m_size_of_esi) * SNPIdxinVec;
		freq = (geno.alleleFreqVec)[SNPIdx];
                invStd = (geno.invstdvVec)[SNPIdx];
		//if(k == 0){
		//	std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
		//}

                while(flag == 0){
//		std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(geno.m_size_of_esi); i++){
                        geno1 = (geno.genoVecofPointers)[indexOfVectorPointer]->at(i);
			//std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        stdGenoMultiMarkers[ind*m_M_Submarker+k] = ((2-(a+b)) - 2*freq)* invStd;;
//			stdGenoMultiMarkers[ind*m_M_Submarker+k] = 2-(a+b);
//			if(k == 0){
    //                    std::cout << "ind*m_M_Submarker+k: " << ind*m_M_Submarker+k << " stdGenoMultiMarkers[ind*m_M_Submarker+k]: " << stdGenoMultiMarkers[ind*m_M_Submarker+k] <<  std::endl;
  //              }


                        indtotal++;
                        ind++;
                        geno1 = geno1 >> 1;

                                if(ind == Nnomissing){
                                        flag = 1;
					break;	
                                }
                        }// end of for(int j=0; j<4; j++){
                    }// end of for(size_t i=Start_idx
                } //end of while(flag == 0){

        }

	//std::cout << "stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] " << stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] << std::endl;

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
        int Msub_mafge1perc;	
  
  	// constructors
  	CorssProd(arma::fcolvec & y)
  		: m_bVec(y) {
  		
  		m_M = geno.getM();
  		m_N = geno.getNnomissing();
  		m_bout.zeros(m_N);
		Msub_mafge1perc=0;
  	} 
  	CorssProd(const CorssProd& CorssProd, Split)
  		: m_bVec(CorssProd.m_bVec)
  	{

  		m_N = CorssProd.m_N;
  		m_M = CorssProd.m_M;
  		m_bout.zeros(m_N);
		Msub_mafge1perc=0;
  	
  	}  
  	// process just the elements of the range I've been asked to
  	void operator()(std::size_t begin, std::size_t end) {
  	  	arma::fvec vec;
  	  	for(unsigned int i = begin; i < end; i++){
			if(geno.alleleFreqVec[i] >= minMAFtoConstructGRM && geno.alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
				geno.Get_OneSNP_StdGeno(i, &vec);
				float val1 = dot(vec,  m_bVec);
				m_bout += val1 * (vec) ;
				Msub_mafge1perc += 1;
			}
		/*	std::cout << "i: " << i << std::endl;
			for(unsigned int j = 0; j < 10; j++){
				std::cout << "m_bVec[j] " << m_bVec[j] << std::endl;
				std::cout << "vec[j] " << vec[j] << std::endl;
			}
		*/
			//m_bout += val1 * (vec) / m_M;
  		}
  	}
  
  	// join my value with that of another InnerProduct
  	void join(const CorssProd & rhs) { 
    		m_bout += rhs.m_bout;
		Msub_mafge1perc += rhs.Msub_mafge1perc; 
  	}
};



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_LOCO : public Worker
{
        // source vectors
        arma::fcolvec & m_bVec;
        unsigned int m_N;
        unsigned int m_Msub;
        unsigned int m_M;
	int startIndex;
	int endIndex;
        // product that I have accumulated
        arma::fvec m_bout;
	unsigned int m_Msub_mafge1perc;

        // constructors
        CorssProd_LOCO(arma::fcolvec & y)
                : m_bVec(y) {

                m_Msub = geno.getMsub(); //LOCO
		startIndex = geno.getStartIndex();
		endIndex = geno.getEndIndex();
                m_M = geno.getM(); //LOCO
                m_N = geno.getNnomissing();
                m_bout.zeros(m_N);
		m_Msub_mafge1perc=0;
        }
        CorssProd_LOCO(const CorssProd_LOCO& CorssProd_LOCO, Split)
                : m_bVec(CorssProd_LOCO.m_bVec)
        {

                m_N = CorssProd_LOCO.m_N;
                m_M = CorssProd_LOCO.m_M;
		m_Msub = CorssProd_LOCO.m_Msub;
		startIndex = geno.getStartIndex();
                endIndex = geno.getEndIndex();
                m_bout.zeros(m_N);
		m_Msub_mafge1perc=0;
        }
	
	   // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
		float val1;
                for(unsigned int i = begin; i < end; i++){
                        geno.Get_OneSNP_StdGeno(i, &vec);
			if(geno.alleleFreqVec[i] >= minMAFtoConstructGRM && geno.alleleFreqVec[i] <= 1-minMAFtoConstructGRM){	
				if(i >= startIndex && i <= endIndex){
					val1 = 0;
					//if(endIndex == 4){
					//		cout << "i: " << i << endl;
					//}
				}else{
                        		val1 = dot(vec,  m_bVec);
				}
                        	m_bout += val1 * (vec);
				m_Msub_mafge1perc += 1;
			}

                }
        }

        // join my value with that of another InnerProduct
        void join(const CorssProd_LOCO & rhs) {
        m_bout += rhs.m_bout;
	m_Msub_mafge1perc += rhs.m_Msub_mafge1perc;	
        }
};






// [[Rcpp::export]]
arma::fvec parallelCrossProd(arma::fcolvec & bVec) {
  
  // declare the InnerProduct instance that takes a pointer to the vector data
  	int M = geno.getM();
 	//int Msub_mafge1perc = geno.getMmafge1perc();
  	CorssProd CorssProd(bVec);
  
  // call paralleReduce to start the work
  	parallelReduce(0, M, CorssProd);
 	
	//cout << "print test; M: " << M << endl;
        //for(int i=0; i<10; ++i)
        //{
        //        cout << (CorssProd.m_bout)[i] << ' ' << endl;
        //        cout << bVec[i] << ' ' << endl;
        //        cout << (CorssProd.m_bout/M)[i] << ' ' << endl;
        //}
        ////cout << endl; 
  // return the computed product
	//std::cout << "number of markers with maf ge " << minMAFtoConstructGRM << " is " << CorssProd.Msub_mafge1perc << std::endl;
  	return CorssProd.m_bout/(CorssProd.Msub_mafge1perc);
  	//return CorssProd.m_bout;
}

// [[Rcpp::export]]
float innerProductFun(std::vector<float> &x, std::vector<float> & y) {
   return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}



// [[Rcpp::export]]
arma::fvec parallelCrossProd_LOCO(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        //int Msub = geno.getMsub();
	int M = geno.getM();

        CorssProd_LOCO CorssProd_LOCO(bVec);

  // call paralleReduce to start the work
        parallelReduce(0, M, CorssProd_LOCO);

  // return the computed product
	//cout << "Msub: " << Msub << endl;
        //for(int i=0; i<100; ++i)
        //{
        //	cout << (CorssProd_LOCO.m_bout/Msub)[i] << ' ';
        //}
	//cout << endl;
        //return CorssProd_LOCO.m_bout/Msub;
        return CorssProd_LOCO.m_bout/(CorssProd_LOCO.m_Msub_mafge1perc);
}


arma::umat locationMat;
arma::vec valueVec;
int dimNum = 0;

// [[Rcpp::export]]
void setupSparseGRM(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    // sparse x sparse -> sparse
    //arma::sp_mat result(a);
    //int r = a.n_rows;
        locationMat.zeros(2,r);
        valueVec.zeros(r);

    locationMat = locationMatinR;
    valueVec = valueVecinR;
    dimNum = r;

    std::cout << locationMat.n_rows << " locationMat.n_rows " << std::endl;
    std::cout << locationMat.n_cols << " locationMat.n_cols " << std::endl;
    std::cout << valueVec.n_elem << " valueVec.n_elem " << std::endl;
    for(size_t i=0; i< 10; i++){
        std::cout << valueVec(i) << std::endl;
        std::cout << locationMat(0,i) << std::endl;
        std::cout << locationMat(1,i) << std::endl;
    }

    //arma::vec y = arma::linspace<arma::vec>(0, 5, r);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    //arma::vec x = arma::spsolve( result, y );

    //return x;
}

bool isUsePrecondM = false;
bool isUseSparseSigmaforInitTau = false;
bool isUseSparseSigmaforModelFitting = false;

// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin(arma::fcolvec& bVec){
       arma::fvec crossProdVec;

	//std::cout << "getCrossprodMatAndKin" << std::endl;

if(isUseSparseSigmaforInitTau | isUseSparseSigmaforModelFitting){
        //cout << "use sparse kinship to estimate initial tau and for getCrossprodMatAndKin" <<  endl;
	arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
	arma::vec x = result * arma::conv_to<arma::dcolvec>::from(bVec);


//double wall3in = get_wall_time();
// double cpu3in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;


    crossProdVec = arma::conv_to<arma::fvec>::from(x);

}else{ 
  	crossProdVec = parallelCrossProd(bVec) ;
}  
  	return(crossProdVec);
}






// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin_LOCO(arma::fcolvec& bVec){

        arma::fvec crossProdVec = parallelCrossProd_LOCO(bVec) ;

        return(crossProdVec);
}


// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
struct indicesRelatedSamples : public RcppParallel::Worker {

  int  Ntotal;
  tbb::concurrent_vector< std::pair<int, int> > &output;

  indicesRelatedSamples(int Ntotal, tbb::concurrent_vector< std::pair<int, int> > &output) : 
    Ntotal(Ntotal), output(output) {} 


  void operator()(std::size_t begin, size_t end) {
    int m_M_Submarker = getSubMarkerNum();
    for(std::size_t k=begin; k < end; k++) {
      int i = (int)(k / Ntotal);
      int j = (int)(k % Ntotal);
      if((j <= i)){
                        i = Ntotal - i - 2;
                        j = Ntotal - j - 1;
      }
      //std::cout << "i,j,k debug: " << i << " " << j << " " << k << std::endl;  
      float kinValueTemp = arma::dot((geno.stdGenoMultiMarkersMat).col(i), (geno.stdGenoMultiMarkersMat).col(j));
      kinValueTemp = kinValueTemp/m_M_Submarker;
      if(kinValueTemp >=  geno.relatednessCutoff) {
        output.push_back( std::pair<int, int>(i, j) );
      }
    }
  }

};


// [[Rcpp::export]]
void printComb(int N){
  int x = N*(N-1)/2 - 1;
  for(std::size_t k=0; k < x; k++) {
      int i = k / N;
      int j = k % N;
      if((j < i)){
                        i = N - i - 2;
                        j = N - j - 1;
      }
     std::cout << "i,j " << i << "," << j << std::endl;
  }

}


//arma::fmat findIndiceRelatedSample(){
//arma::fmat findIndiceRelatedSample(){

// [[Rcpp::export]]
void findIndiceRelatedSample(){

  int Ntotal = geno.getNnomissing(); 
//  tbb::concurrent_vector< std::pair<float, float> > output;

//  indicesRelatedSamples indicesRelatedSamples(Ntotal,output);
  indicesRelatedSamples indicesRelatedSamples(Ntotal,geno.indiceVec);

  long int Ntotal2 = (long int)Ntotal;

  long int totalCombination = Ntotal2*(Ntotal2-1)/2 - 1;
  std::cout << "Ntotal: " << Ntotal << std::endl;
  std::cout << std::numeric_limits<int>::max() << std::endl;
  std::cout << std::numeric_limits<long int>::max() << std::endl;
  std::cout << std::numeric_limits<long long int>::max() << std::endl;
  std::cout << "totalCombination: " << totalCombination << std::endl;
  long int x = 1000001;
  int b = (int)(x / Ntotal);
  int a = (int)(x % Ntotal);
  std::cout << "a " << a << std::endl;
  std::cout << "b " << b << std::endl;
  
  parallelFor(0, totalCombination, indicesRelatedSamples);

//  arma::fmat xout(output.size()+Ntotal,2);

//  for(int i=0; i<output.size(); i++) {
//    xout(i,0) = output[i].first;
//    xout(i,1) = output[i].second;
//  }
//  for(int i=output.size(); i < output.size()+Ntotal; i++) {
//    xout(i,0) = i - output.size();
//    xout(i,1) = xout(i,0);
//  }

/*
  for(int i=0; i < Ntotal; i++){
    (geno.indiceVec).push_back( std::pair<int, int>(i, i) );
  }
*/

//  return(xout);
}



struct sparseGRMUsingOneMarker : public Worker {
   // input matrix to read from
  // arma::imat & iMat;
   // output matrix to write to
   arma::fvec & GRMvec;

   int M = geno.getM();
   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
//   sparseGRMUsingOneMarker(arma::imat & iMat, arma::fvec &GRMvec)
//      : iMat(iMat), GRMvec(GRMvec) {}


  sparseGRMUsingOneMarker(arma::fvec &GRMvec)
      : GRMvec(GRMvec) {}


   // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
            // rows we will operate on
//            int iint = iMat(i,0);
//            int jint = iMat(i,1);
	   int iint = (geno.indiceVec)[i].first;	
	   int jint = (geno.indiceVec)[i].second;	
/*
            float ival = geno.m_OneSNP_StdGeno(iint);
            float jval = geno.m_OneSNP_StdGeno(jint);
            // write to output matrix
            //rmat(i,j) = sqrt(.5 * (d1 + d2));
            GRMvec(i) = ival*jval/M;
*/
	//use Look-Up table for calucate GRMvec(i)
	    int ival = geno.m_OneSNP_Geno(iint);	
	    int jval = geno.m_OneSNP_Geno(jint);
	    GRMvec(i) = geno.sKinLookUpArr[ival][jval]; 

      }
   }
};

//void parallelcalsparseGRM(arma::imat & iMat, arma::fvec &GRMvec) {

// [[Rcpp::export]]
void parallelcalsparseGRM(arma::fvec &GRMvec) {

//  int n1 = geno.indiceVec.size();
  // allocate the output matrix
  //GRMvec.set_size(n1);
//  std::cout << "OKKK3: "  << std::endl;
//  sparseGRMUsingOneMarker sparseGRMUsingOneMarker(iMat, GRMvec);
  sparseGRMUsingOneMarker sparseGRMUsingOneMarker(GRMvec);
//  std::cout << "OKKK4: "  << std::endl;

//  std::cout << "n1 " << n1 << std::endl;
//  std::cout << "iMat.n_cols " << iMat.n_cols << std::endl;
  // call parallelFor to do the work
//  parallelFor(0, iMat.n_rows, sparseGRMUsingOneMarker);
  parallelFor(0, (geno.indiceVec).size(), sparseGRMUsingOneMarker);

  // return the output matrix
  // return GRMvec;
}


struct sumTwoVec : public Worker
{   
   // source vectors
   arma::fvec &x;
   
   arma::fvec &sumVec;
  
   int M = geno.getM(); 
   // constructors
   sumTwoVec(arma::fvec &x,arma::fvec &sumVec) 
      : x(x), sumVec(sumVec) {}
   
     // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
            // rows we will operate on
            sumVec(i) = x(i)+(geno.kinValueVecFinal)[i];
	    (geno.kinValueVecFinal)[i] = sumVec(i);	
      }
   }
   
};

// [[Rcpp::export]]
void  parallelsumTwoVec(arma::fvec &x) {

  int n1 = x.n_elem;
  // allocate the output matrix
  arma::fvec sumVec;
  sumVec.set_size(n1);

  sumTwoVec sumTwoVec(x, sumVec);

  // call parallelFor to do the work
  parallelFor(0, x.n_elem, sumTwoVec);

}



// [[Rcpp::export]]
void setgeno(std::string genofile, std::vector<int> & subSampleInGeno, float memoryChunk, bool isDiagofKinSetAsOne)
{
	int start_s=clock();
        geno.setGenoObj(genofile, subSampleInGeno, memoryChunk, isDiagofKinSetAsOne);
	//geno.printAlleleFreqVec();
	//geno.printGenoVec();
	int stop_s=clock();
	cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
}




// [[Rcpp::export]]
arma::ivec Get_OneSNP_Geno(int SNPIdx)
{

	arma::ivec temp = * geno.Get_OneSNP_Geno(SNPIdx);
	return(temp);

}


// [[Rcpp::export]]
arma::ivec Get_OneSNP_Geno_forVarianceRatio(int SNPIdx)
{
       
        arma::ivec temp = * geno.Get_OneSNP_Geno(SNPIdx);
        return(temp);

}



// [[Rcpp::export]]
arma::fvec Get_OneSNP_StdGeno(int SNPIdx)
{

	arma::fvec temp; 
	geno.Get_OneSNP_StdGeno(SNPIdx, & temp);
//	for(int j = 0; j < 100; j++){
//                std::cout << "temp(j): " << j << " " << temp(j) << std::endl;

 //       }


	return(temp);

}
  
    
  

//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
// [[Rcpp::export]]
arma::fvec getDiagOfSigma(arma::fvec& wVec, arma::fvec& tauVec){
  
	int Nnomissing = geno.getNnomissing();
	int M = geno.getM();
	int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
		
	//cout << "N=" << N << endl;
	arma::fvec diagVec(Nnomissing);
	float diagElement;
	float floatBuffer;
  	//float minvElement;
  
  	if(!(geno.setKinDiagtoOne)){ 
	  diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /MminMAF + tauVec(0)/wVec;

	}else{
	  diagVec = tauVec(1) + tauVec(0)/wVec;
	}
	//cout << "wVec: " << endl;
	//wVec.print();	
	//std::cout << "M " << M << std::endl;
	//std::cout << "tauVec(0) " << tauVec(0) << std::endl;
	//std::cout << "tauVec(1) " << tauVec(1) << std::endl;
       // for(unsigned int i=0; i< Nnomissing; i++){
	//	std::cout << "(*geno.Get_Diagof_StdGeno()) /M: " << (*geno.Get_Diagof_StdGeno()) /M << std::endl;
	//}

	//make diag of kin to be 1 to compare results of emmax and gmmat
	//diagVec = tauVec(1) + tauVec(0)/wVec;
	for(unsigned int i=0; i< Nnomissing; i++){
//	if(i < 100){
//		std::cout << i << "th element of diag of sigma and wVec " << diagVec(i) << " " << wVec(i) << std::endl;
//	}
  		if(diagVec(i) < 1e-4){
  			diagVec(i) = 1e-4 ;
  		}
  	}
  


    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
  	return(diagVec);
}

// [[Rcpp::export]]
arma::fvec getDiagOfSigma_LOCO(arma::fvec& wVec, arma::fvec& tauVec){

        int Nnomissing = geno.getNnomissing();
        int Msub = geno.getMsub();

        //cout << "N=" << N << endl;
        arma::fvec diagVec(Nnomissing);
        float diagElement;
        float floatBuffer;
        //float minvElement;
	int Msub_MAFge_minMAFtoConstructGRM = geno.getMsub_MAFge_minMAFtoConstructGRM_in();
        diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno_LOCO()) /(Msub_MAFge_minMAFtoConstructGRM) + tauVec(0)/wVec;
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
arma::fvec getDiagOfSigma_surv(arma::fvec& diagofWminusUinv, arma::fvec& tauVec){
  
	int Nnomissing = geno.getNnomissing();
	int M = geno.getM();
	int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
		
	//cout << "N=" << N << endl;
	arma::fvec diagVec(Nnomissing);
	float diagElement;
	float floatBuffer;
  	//float minvElement;
  
  	if(!(geno.setKinDiagtoOne)){ 
	  diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /MminMAF + tauVec(0) * diagofWminusUinv;

	}else{
	  diagVec = tauVec(1) + tauVec(0) * diagofWminusUinv;
	}
	//cout << "wVec: " << endl;
	//wVec.print();	
	//std::cout << "M " << M << std::endl;
	//std::cout << "tauVec(0) " << tauVec(0) << std::endl;
	//std::cout << "tauVec(1) " << tauVec(1) << std::endl;
       // for(unsigned int i=0; i< Nnomissing; i++){
	//	std::cout << "(*geno.Get_Diagof_StdGeno()) /M: " << (*geno.Get_Diagof_StdGeno()) /M << std::endl;
	//}

	//make diag of kin to be 1 to compare results of emmax and gmmat
	//diagVec = tauVec(1) + tauVec(0)/wVec;
	for(unsigned int i=0; i< Nnomissing; i++){
//	if(i < 100){
//		std::cout << i << "th element of diag of sigma and wVec " << diagVec(i) << " " << wVec(i) << std::endl;
//	}
  		if(diagVec(i) < 1e-4){
  			diagVec(i) = 1e-4 ;
  		}
  	}
  


    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
  	return(diagVec);
}

// [[Rcpp::export]]
arma::fvec getDiagOfSigma_surv_LOCO(arma::fvec& diagofWminusUinv, arma::fvec& tauVec){

        int Nnomissing = geno.getNnomissing();
        int Msub = geno.getMsub();

        //cout << "N=" << N << endl;
        arma::fvec diagVec(Nnomissing);
        float diagElement;
        float floatBuffer;
        //float minvElement;
	int Msub_MAFge_minMAFtoConstructGRM = geno.getMsub_MAFge_minMAFtoConstructGRM_in();
        diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno_LOCO()) /(Msub_MAFge_minMAFtoConstructGRM) + tauVec(0)*diagofWminusUinv;
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
	
	//for(int j = 0; j < 10; j++){
        //        std::cout << "bVec(j): " << bVec(j) << std::endl;
        //        std::cout << "crossProd1(j): " << crossProd1(j) << std::endl;

        //}


        crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;

	//for(int j = 0; j < 10; j++){
        //        std::cout << "crossProdVec(j): " << j << " " << crossProdVec(j) << std::endl;

        //}A



        return(crossProdVec);
}

// [[Rcpp::export]]
arma::fvec extractVecatTimek(unsigned int ktime , arma::fvec & rvecIndex, arma::fvec & winvn) {
	arma::fvec kthVec;
	unsigned int n_kthVec = winvn.n_elem;
	kthVec.zeros(n_kthVec);
	for(unsigned int i=0; i< n_kthVec; i++){
		if(rvecIndex(i) >= ktime){
			kthVec(i) = winvn(i);
		} 	
        }
	return(kthVec);
}



//extractUvecforkthTime(i, n_RvecIndex, n_NVec, n_sqrtDVec, vec);
// [[Rcpp::export]]
void extractUvecforkthTime(unsigned int kthtime, arma::fvec & RvecIndex,  arma::fvec& NVec,  arma::fvec & sqrtDVec, arma::fvec & kthVec){
        //unsigned int ktime=RvecIndex(nthsample);
        unsigned int nsample = RvecIndex.n_elem;
        kthVec.zeros();
        float sqrtDKth = sqrtDVec(kthtime);
        for(unsigned int j = 0; j < nsample; j++){
                if((RvecIndex(j)-1) >= kthtime){
                        kthVec(j) = NVec(j);
                }
        }
        kthVec = kthVec * sqrtDKth;
}



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_UandbVec_surv : public Worker
{
        // source vectors
        arma::fcolvec n_bVec;
        arma::fcolvec n_RvecIndex;
        arma::fcolvec n_NVec;
        arma::fcolvec n_sqrtDVec;
        //unsigned int k_uniqTime;
        // product that I have accumulated
        arma::fvec m_bout;
        unsigned int m_N;

        // constructors
        CorssProd_UandbVec_surv(arma::fcolvec & x, arma::fvec & y,  arma::fvec & z,  arma::fvec & q)
                : n_bVec(x),n_RvecIndex(y),n_NVec(z),n_sqrtDVec(q) {
                  m_N = geno.getNnomissing();
                  m_bout.zeros(m_N);
        }
        CorssProd_UandbVec_surv(const CorssProd_UandbVec_surv& CorssProd_UandbVec_surv, Split)
                : n_bVec(CorssProd_UandbVec_surv.n_bVec),n_RvecIndex(CorssProd_UandbVec_surv.n_RvecIndex),n_NVec(CorssProd_UandbVec_surv.n_NVec),n_sqrtDVec(CorssProd_UandbVec_surv.n_sqrtDVec)
        {
                m_N = CorssProd_UandbVec_surv.m_N;
                m_bout.zeros(m_N);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                vec.zeros(m_N);
                //int nthsample;
                //int ktime;
                //arma::fvec vec.zeros(m_N);
                for(unsigned int i = begin; i < end; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractUvecforkthTime(i, n_RvecIndex, n_NVec, n_sqrtDVec, vec);
                        //for(unsigned int j = 0; j < ktime; j++){
                        //        vec(j) = n_Dvec(j)*n_sqrtWinvNVec(i);
                        //}
                        float val1 = dot(vec,  n_bVec);
                        m_bout += val1 * (vec);
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_UandbVec_surv & rhs) {
                m_bout += rhs.m_bout;
        }
};




// [[Rcpp::export]]
arma::fvec parallelCrossProd_UandbVec_surv(arma::fcolvec & bVec, arma::fvec & RvecIndex, arma::fvec& NVec,  arma::fvec & sqrtDVec) {

//  // declare the InnerProduct instance that takes a pointer to the vector data
        unsigned int ktime = sqrtDVec.n_elem;
        CorssProd_UandbVec_surv  CorssProd_UandbVec_surv(bVec, RvecIndex, NVec, sqrtDVec);
        //int m_N = geno.getNnomissing();

//  // call paralleReduce to start the work
        parallelReduce(0, ktime, CorssProd_UandbVec_surv);

        return CorssProd_UandbVec_surv.m_bout;
}



// [[Rcpp::export]]
arma::fcolvec getProdWminusUb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex, arma::fvec& NVec, arma::fvec& sqrtDVec, arma::fvec& wVec){
        //unsigned int nsample = geno.getNnomissing();
        //unsigned int kuniqtime = Dvec.n_elem;

        arma::fcolvec Ub = parallelCrossProd_UandbVec_surv(bVec, RvecIndex, NVec, sqrtDVec);
        arma::fcolvec WminusUb = wVec % bVec - Ub;
        return WminusUb;
}

/*
arma::fcolvec dofWminusU;
// [[Rcpp::export]]
void getDiagofWminusU(arma::fcolvec& diagVec){
	arma::fcolvec dofWminusU = diagVec;
	cout << "dofWminusU" << endl;
	dofWminusU.print();
}
*/
/*
// [[Rcpp::export]]
arma::fvec extractVecatSamplen(unsigned int nthsample , unsigned int ktime, arma::fvec & rvecIndex, arma::fvec & winvn) {
	arma::fvec nthVec;
	nthVec.zeros(ktime);
	
	for(unsigned int i=0; i< ktime; i++){
		if(rvecIndex(nthsample) >= i){
			nthVec(i) = winvn(rvecIndex(nthsample)-1);
		} 	
        }
	return(kthVec);
}
*/

// [[Rcpp::export]]
arma::fcolvec getCrossprod_Surv(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat & WinvNRt, arma::fmat & ACinv){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fmat WinvNRtG;
        arma::fmat ACivWinvNRtG;
 	//cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
 	//cout << "OKKKKK4" << endl;
        WinvNRtG = (WinvNRt.t()) * bVec;
 	//cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
 	//cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
 	//cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);

        crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

        return(crossProdVec);
}



// [[Rcpp::export]]
arma::fcolvec getCrossprod_Surv_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat & WinvNRt, arma::fmat & ACinv){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fmat WinvNRtG;
        arma::fmat ACivWinvNRtG;
 	//cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
 	//cout << "OKKKKK4" << endl;
        WinvNRtG = (WinvNRt.t()) * bVec;
 	//cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
 	//cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
 	//cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);

        crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

        return(crossProdVec);
}

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_WinvNRttandVec : public Worker
{
        // source vectors
        arma::fcolvec & n_bVec;
	arma::fcolvec & n_RvecIndex;
	arma::fcolvec & n_WinvN;

        unsigned int k_uniTime;

        // product that I have accumulated
        arma::fvec m_bout;


        // constructors
        CorssProd_WinvNRttandVec(arma::fcolvec & x, arma::fvec & y, arma::fvec & z, unsigned int k)
                : n_bVec(x),n_RvecIndex(y),n_WinvN(z),k_uniTime(k) {
                m_bout.zeros(k_uniTime);
        }
        CorssProd_WinvNRttandVec(const CorssProd_WinvNRttandVec& CorssProd_WinvNRttandVec, Split)
                : n_bVec(CorssProd_WinvNRttandVec.n_bVec),n_RvecIndex(CorssProd_WinvNRttandVec.n_RvecIndex),n_WinvN(CorssProd_WinvNRttandVec.n_WinvN),k_uniTime(CorssProd_WinvNRttandVec.k_uniTime)
        {

                m_bout.zeros(k_uniTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                float val1;
                int ktime;
                for(unsigned int i = begin; i < end; i++){
                        ktime = i;
			vec=extractVecatTimek(ktime, n_RvecIndex, n_WinvN);
//			std::cout << "j: " << j << std::endl;	
                        val1 = dot(vec,  n_bVec);
                        m_bout[i] += m_bout[i] + val1;
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_WinvNRttandVec & rhs) {
        m_bout += rhs.m_bout;
        }
};





/*
// [[Rcpp::export]]
arma::fvec parallelCrossProd_WinvNRttandVec(arma::fcolvec & bVec, arma::fvec & RvecIndex, arma::fvec & WinvN, unsigned int kuniqtime) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        CorssProd_WinvNRttandVec  CorssProd_WinvNRttandVec(bVec, RvecIndex, WinvN, kuniqtime);

  // call paralleReduce to start the work
        CorssProd_WinvNRttandVec(0, kuniqtime, CorssProd_WinvNRttandVec);

        //cout << "print test; M: " << M << endl;
        //for(int i=0; i<10; ++i)
        //{
        //        cout << (CorssProd.m_bout)[i] << ' ' << endl;
        //        cout << bVec[i] << ' ' << endl;
        //        cout << (CorssProd.m_bout/M)[i] << ' ' << endl;
        //}
        ////cout << endl;
  // return the computed product
        //std::cout << "number of markers with maf ge " << minMAFtoConstructGRM << " is " << CorssProd.Msub_mafge1perc << std::endl;
        return CorssProd_WinvNRttandVec.m_bout;
        //return CorssProd.m_bout;
}



// [[Rcpp::export]]
arma::fcolvec getprodWinvNRttandVec(arma::fcolvec& bVec, arma::fvec & RvecIndex, arma::fvec & WinvN, unsigned int kuniqtime){
	arma::fcolvec WinvNRttandVec = parallelCrossProd_WinvNRttandVec(bVec, RvecIndex, WinvN,kuniqtime);	
	return(WinvNRttandVec);	
}

// [[Rcpp::export]]
arma::fcolvec getprodDRNWinvbvec(arma::fcolvec & bVec, arma::fcolvec & DVec, arma::fvec & RvecIndex, arma::fvec & NWinv, unsigned int kuniqtime){
	
	arma::fcolvec DRNWinv.zeros(kuniqtime);
		


	arma::fcolvec cumuSumVec.zeros(bVec.n_elem);
	cumuSumVec(0) = NWinv(0);
	for(unsigned int j = 1; j < kuniqtime; i++){
		cumuSumVec(j) = NWinv(j) + cumuSumVec(j-1); 	
	}
	
	for(unsigned int i = 0; i < kuniqtime; i++){
                        
                        vec=extractVecatSamplen(ktime,nthsample,n_RvecIndex, n_sqrtWinvNVec )
//                      std::cout << "j: " << j << std::endl;
                        val1 = dot(vec,  n_bVec);
                        m_bout += val1 * (vec);
	}	
	
	WinvNRttandVec(arma::fcolvec& bVec, arma::fvec & RvecIndex, arma::fvec & WinvN, unsigned int kuniqtime){
        arma::fcolvec WinvNRttandVec = parallelCrossProd_WinvNRttandVec(bVec, RvecIndex, WinvN,kuniqtime);
        return(WinvNRttandVec);
}


*/

// [[Rcpp::export]]
void extractVecfornthSample(unsigned int nthsample, unsigned int k_uniqTime, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & nthVec) {
	unsigned int ktime=RvecIndex(nthsample);
	nthVec.zeros(k_uniqTime);
        for(unsigned int j = 0; j < ktime; j++){
        	nthVec(j) = sqrtWinvNVec(nthsample);
        }
}


// [[Rcpp::export]]
void extractVecfornthSample_double(unsigned int nthsample, unsigned int k_uniqTime, arma::vec & RvecIndex, arma::vec & sqrtWinvNVec, arma::vec & nthVec) {
	unsigned int ktime=RvecIndex(nthsample);
	nthVec.zeros(k_uniqTime);
        for(unsigned int j = 0; j < ktime; j++){
        	nthVec(j) = sqrtWinvNVec(nthsample);
        }
}



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_RandbVec_surv : public Worker
{
        // source vectors
        arma::fcolvec & n_bVec;
	arma::fcolvec & n_RvecIndex;
        unsigned int k_uniqTime;
        // product that I have accumulated
        arma::fvec m_bout;
	unsigned int m_N;

        // constructors
        CorssProd_RandbVec_surv(arma::fcolvec & x, arma::fvec & y, unsigned int k)
                : n_bVec(x),n_RvecIndex(y),k_uniqTime(k) {
		  m_N = geno.getNnomissing();
                  m_bout.zeros(k_uniqTime);
        }
        CorssProd_RandbVec_surv(const CorssProd_RandbVec_surv& CorssProd_RandbVec_surv, Split)
                : n_bVec(CorssProd_RandbVec_surv.n_bVec),n_RvecIndex(CorssProd_RandbVec_surv.n_RvecIndex),k_uniqTime(CorssProd_RandbVec_surv.k_uniqTime)
        {
		m_N = CorssProd_RandbVec_surv.m_N;
                m_bout.zeros(k_uniqTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
		vec.zeros(k_uniqTime);
                float val1;
                int nthsample;
		int ktime;
                for(unsigned int i = begin; i < end; i++){
                        nthsample = i;
			ktime=n_RvecIndex(i);
			for(unsigned int j = 0; j < ktime; j++){
				m_bout(j) += n_bVec(i); 
			}
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_RandbVec_surv & rhs) {
        m_bout += rhs.m_bout;
        }
};




//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_AandbVec_surv : public Worker
{
        // source vectors
        arma::fcolvec n_bVec;
        arma::fcolvec n_RvecIndex;
	arma::fcolvec n_sqrtWinvNVec;
        unsigned int k_uniqTime;
        // product that I have accumulated
        arma::fvec m_bout;
        unsigned int m_N;

        // constructors
        CorssProd_AandbVec_surv(arma::fcolvec & x, arma::fvec & y,  arma::fvec & z,  unsigned int k)
                : n_bVec(x),n_RvecIndex(y),n_sqrtWinvNVec(z),k_uniqTime(k) {
                  m_N = geno.getNnomissing();
                  m_bout.zeros(k_uniqTime);
        }
        CorssProd_AandbVec_surv(const CorssProd_AandbVec_surv& CorssProd_AandbVec_surv, Split)
                : n_bVec(CorssProd_AandbVec_surv.n_bVec),n_RvecIndex(CorssProd_AandbVec_surv.n_RvecIndex),n_sqrtWinvNVec(CorssProd_AandbVec_surv.n_sqrtWinvNVec),k_uniqTime(CorssProd_AandbVec_surv.k_uniqTime)
        {
                m_N = CorssProd_AandbVec_surv.m_N;
                m_bout.zeros(k_uniqTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
		vec.zeros(k_uniqTime);
                //int nthsample;
                //int ktime;
		//arma::fvec vec.zeros(m_N);
                for(unsigned int i = begin; i < end; i++){
			//nthsample = i;
			//ktime=n_RvecIndex(i);
			//vec.zeros(k_uniqTime);
			extractVecfornthSample(i, k_uniqTime, n_RvecIndex, n_sqrtWinvNVec, vec);	
			//for(unsigned int j = 0; j < ktime; j++){
                        //        vec(j) = n_Dvec(j)*n_sqrtWinvNVec(i);
                        //}
                        float val1 = dot(vec,  n_bVec);
                        m_bout += val1 * (vec);
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_AandbVec_surv & rhs) {
        	m_bout += rhs.m_bout;
        }
};



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_AandbVec_surv_double : public Worker
{
        // source vectors
        arma::colvec n_bVec;
        arma::colvec n_RvecIndex;
	arma::colvec n_sqrtWinvNVec;
        unsigned int k_uniqTime;
        // product that I have accumulated
        arma::vec m_bout;
        unsigned int m_N;

        // constructors
        CorssProd_AandbVec_surv_double(arma::colvec & x, arma::vec & y,  arma::vec & z,  unsigned int k)
                : n_bVec(x),n_RvecIndex(y),n_sqrtWinvNVec(z),k_uniqTime(k) {
                  m_N = geno.getNnomissing();
                  m_bout.zeros(k_uniqTime);
        }
        CorssProd_AandbVec_surv_double(const CorssProd_AandbVec_surv_double& CorssProd_AandbVec_surv_double, Split)
                : n_bVec(CorssProd_AandbVec_surv_double.n_bVec),n_RvecIndex(CorssProd_AandbVec_surv_double.n_RvecIndex),n_sqrtWinvNVec(CorssProd_AandbVec_surv_double.n_sqrtWinvNVec),k_uniqTime(CorssProd_AandbVec_surv_double.k_uniqTime)
        {
                m_N = CorssProd_AandbVec_surv_double.m_N;
                m_bout.zeros(k_uniqTime);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::vec vec;
		vec.zeros(k_uniqTime);
                //int nthsample;
                //int ktime;
		//arma::fvec vec.zeros(m_N);
                for(unsigned int i = begin; i < end; i++){
			//nthsample = i;
			//ktime=n_RvecIndex(i);
			//vec.zeros(k_uniqTime);
			extractVecfornthSample_double(i, k_uniqTime, n_RvecIndex, n_sqrtWinvNVec, vec);	
			//for(unsigned int j = 0; j < ktime; j++){
                        //        vec(j) = n_Dvec(j)*n_sqrtWinvNVec(i);
                        //}
                        double val1 = dot(vec,  n_bVec);
                        m_bout += val1 * (vec);
                }
        }
        // join my value with that of another InnerProduct
        void join(const  CorssProd_AandbVec_surv_double & rhs) {
        	m_bout += rhs.m_bout;
        }
};


// [[Rcpp::export]]
arma::fvec parallelCrossProd_AandbVec_surv(arma::fcolvec & bVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, unsigned int kuniqtime) {

//  // declare the InnerProduct instance that takes a pointer to the vector data
//	unsigned int ktime = sqrtWinvNVec.n_elem;
        CorssProd_AandbVec_surv  CorssProd_AandbVec_surv(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
	int m_N = geno.getNnomissing();

//  // call paralleReduce to start the work
        parallelReduce(0, m_N, CorssProd_AandbVec_surv);

        return CorssProd_AandbVec_surv.m_bout;
}


// [[Rcpp::export]]
arma::vec parallelCrossProd_AandbVec_surv_double(arma::colvec & bVec, arma::vec & RvecIndex, arma::vec & sqrtWinvNVec, unsigned int kuniqtime) {

//  // declare the InnerProduct instance that takes a pointer to the vector data
//	unsigned int ktime = sqrtWinvNVec.n_elem;
        CorssProd_AandbVec_surv_double  CorssProd_AandbVec_surv_double(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
	int m_N = geno.getNnomissing();

//  // call paralleReduce to start the work
        parallelReduce(0, m_N, CorssProd_AandbVec_surv_double);

        return CorssProd_AandbVec_surv_double.m_bout;
}







// [[Rcpp::export]]
arma::fvec parallelCrossProd_RandbVec_surv(arma::fcolvec & bVec, arma::fvec & RvecIndex, unsigned int kuniqtime) {
	
  // declare the InnerProduct instance that takes a pointer to the vector data
        CorssProd_RandbVec_surv  CorssProd_RandbVec_surv(bVec, RvecIndex, kuniqtime);
	int m_N = geno.getNnomissing();
  // call paralleReduce to start the work
        parallelReduce(0, m_N, CorssProd_RandbVec_surv);

        return CorssProd_RandbVec_surv.m_bout;
}

// [[Rcpp::export]]
arma::fcolvec getProdRb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex, unsigned int kuniqtime){
        //unsigned int nsample = geno.getNnomissing();
        //unsigned int kuniqtime = Dvec.n_elem;

        arma::fcolvec Rb = parallelCrossProd_RandbVec_surv(bVec, RvecIndex, kuniqtime);
        return Rb;
}




// [[Rcpp::export]]
arma::fcolvec getProdAb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex, arma::fvec& sqrtWinvNVec,arma::fvec& Dvec){
	//unsigned int nsample = geno.getNnomissing();
	unsigned int kuniqtime = Dvec.n_elem;
	
	arma::fcolvec Ab = parallelCrossProd_AandbVec_surv(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
	//cout << "Ab 1st part " << endl;
	//Ab.print();
	Ab = Ab +  (-1/Dvec) % bVec;
	return Ab;
}




// [[Rcpp::export]]
arma::colvec getProdAb_Surv_double(arma::colvec& bVec, arma::vec & RvecIndex, arma::vec& sqrtWinvNVec, arma::vec& Dvec){
        //unsigned int nsample = geno.getNnomissing();
        unsigned int kuniqtime = Dvec.n_elem;

        arma::colvec Ab = parallelCrossProd_AandbVec_surv_double(bVec, RvecIndex, sqrtWinvNVec, kuniqtime);
        Ab = Ab +  (-1/Dvec) % bVec;
        return Ab;
}







// [[Rcpp::export]]
arma::fvec getDiagofA( arma::fvec& RvecIndex, arma::fvec& sqrtWinvNVec,arma::fvec& Dvec){
	arma::fvec diagA;
	diagA = (-1/Dvec);
	unsigned int nsample = sqrtWinvNVec.n_elem;
	arma::fvec vec;	
	unsigned int k_uniqTime = Dvec.n_elem;

	for(unsigned int i = 0; i < nsample; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractVecfornthSample(i, k_uniqTime, RvecIndex, sqrtWinvNVec, vec);
			//cout << "vec.n_elem: " << vec.n_elem << endl;
			diagA = diagA + vec % vec;
        }


	for(unsigned int i = 0; i < k_uniqTime; i++){
		if(diagA(i) == 0.0){
			diagA(i) = 0.0001;
		}
        }

	return(diagA);
}


// [[Rcpp::export]]
arma::vec getDiagofA_double( arma::vec& RvecIndex, arma::vec& sqrtWinvNVec,arma::vec& Dvec){
	arma::vec diagA;
	diagA = (-1/Dvec);
	unsigned int nsample = sqrtWinvNVec.n_elem;
	arma::vec vec;	
	unsigned int k_uniqTime = Dvec.n_elem;

	for(unsigned int i = 0; i < nsample; i++){
                        //nthsample = i;
                        //ktime=n_RvecIndex(i);
                        //vec.zeros(k_uniqTime);
                        extractVecfornthSample_double(i, k_uniqTime, RvecIndex, sqrtWinvNVec, vec);
			//cout << "vec.n_elem: " << vec.n_elem << endl;
			diagA = diagA + vec % vec;
        }


	for(unsigned int i = 0; i < k_uniqTime; i++){
		if(diagA(i) == 0.0){
			diagA(i) = 0.0001;
		}
        }

	return(diagA);

}


// [[Rcpp::export]]
arma::vec getPCG1ofACinvAndVector_test(arma::vec& bVec,  arma::vec& RvecIndex, arma::vec& sqrtWinvNVec,arma::vec& Dvec, int maxiterPCG, float tolPCG, arma::vec & wVec, arma::vec & tauVec, arma::mat & Rmat){
    unsigned int kuniqtime = Dvec.n_elem;
    //cout << "kuniqtime is " << kuniqtime << endl;
    arma::vec xVec(kuniqtime);
    xVec.zeros();
	//bVec = bVec/(1e+3);
        arma::vec rVec = bVec;
        arma::vec r1Vec;
        arma::vec zVec(kuniqtime);
        arma::vec minvVec(kuniqtime);

        minvVec = 1/getDiagofA_double(RvecIndex,sqrtWinvNVec,Dvec);/////To update
        zVec = minvVec % rVec;
	cout << "minvVec(10): " << minvVec(10) << endl;
	cout << "minvVec(20): " << minvVec(20) << endl;
	//zVec = rVec;
        double sumr2 = sum(rVec % rVec);
        arma::vec z1Vec(kuniqtime);
        arma::vec pVec = zVec;
	cout << "Rmat.n_cols " << Rmat.n_cols << endl;
        cout << "Rmat.n_rows " << Rmat.n_rows << endl;
        cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
	//arma::fmat sqrtWinvNmat = arma::diagmat(sqrtWinvNVec);
	// cout << "OK" << endl;
	//arma::fmat ApVectemp = (Rmat.t()) * sqrtWinvNmat;
        //arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
	arma::colvec ApVec = getProdAb_Surv_double(pVec,RvecIndex,sqrtWinvNVec,Dvec);
        //cout << "ApVec(10): " << ApVec(10) << endl;
        //cout << "ApVec0(10): " << ApVec0(10) << endl;
        int iter = 0;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                arma::colvec ApVec = getProdAb_Surv_double(pVec,RvecIndex,sqrtWinvNVec,Dvec);
		cout << "iter: " << iter << endl;	
		//arma::fcolvec ApVectemp = (Rmat.t()) * sqrtWinvNVec;
		//arma::fcolvec ApVec = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
		//cout << "Rmat.n_cols " << Rmat.n_cols << endl;
		//cout << "Rmat.n_rows " << Rmat.n_rows << endl;
		//cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
		//arma::fmat ApVectemp = (Rmat.t()) * sqrtWinvNmat; 	
		//arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
		cout << "ApVec(10): " << ApVec(10) << endl;
		cout << "pVec(10): " << pVec(10) << endl;
		//cout << "ApVec0(10): " << ApVec0(10) << endl;
                arma::vec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);
		cout << "rVec.t() * zVec " << rVec.t() * zVec << endl;
		cout << "pVec.t() * ApVec " << pVec.t() * ApVec << endl;
                float a = preA(0);
		cout << "a: " << a << endl;
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
		arma::vec z1Vec = minvVec % r1Vec;
		//arma::fvec z1Vec = r1Vec;
		arma::vec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                double bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
		cout << "bet: " << bet << endl;	
		cout << "Prebet.n_elem: " << Prebet.n_elem << endl;	
		cout << "z1Vec(10): " << z1Vec(10) << endl;
		cout << "r1Vec(10): " << r1Vec(10) << endl;


		zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
		cout << "sumr2 is " << sumr2 << endl;
        }
        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;
        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
	//xVec = xVec *(1e+3);
        return(xVec);
}




// [[Rcpp::export]]
arma::fvec getPCG1ofACinvAndVector(arma::fvec& bVec,  arma::fvec& RvecIndex, arma::fvec& sqrtWinvNVec,arma::fvec& Dvec, int maxiterPCG, float tolPCG, arma::fvec & wVec, arma::fvec & tauVec){
    maxiterPCG = 200;	
    unsigned int kuniqtime = Dvec.n_elem;
    //cout << "kuniqtime is " << kuniqtime << endl;
    arma::fvec xVec(kuniqtime);
    xVec.zeros();

	//bVec = bVec /(1e+3);	
        arma::fvec rVec = bVec;
        arma::fvec r1Vec;
        arma::fvec zVec(kuniqtime);
        arma::fvec minvVec(kuniqtime);

        minvVec = 1/getDiagofA(RvecIndex,sqrtWinvNVec,Dvec);/////To update
        zVec = minvVec % rVec;
	//cout << "minvVec(10): " << minvVec(10) << endl;
	//cout << "minvVec(20): " << minvVec(20) << endl;
	//zVec = rVec;
        float sumr2 = sum(rVec % rVec);
        arma::fvec z1Vec(kuniqtime);
        arma::fvec pVec = zVec;

	//arma::fcolvec ApVec = getProdAb_Surv(pVec,RvecIndex,sqrtWinvNVec,Dvec);


	//cout << "RmatIndex.n_cols " << RmatIndex.n_cols << endl;
        //cout << "RmatIndex.n_rows " << RmatIndex.n_rows << endl;
        //cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
        //arma::fcolvec ApVectemp = (Rmat.t()) * sqrtWinvNVec;
        //arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
        //cout << "ApVec(10): " << ApVec(10) << endl;
        //cout << "ApVec0(10): " << ApVec0(10) << endl;



        int iter = 0;
	arma::fcolvec ApVec;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                ApVec = getProdAb_Surv(pVec,RvecIndex,sqrtWinvNVec,Dvec);
		//cout << "iter: " << iter << endl;	
		//for(size_t j=0; j< 10; j++){
		//	cout << "j: " << j << " ApVec(j) " << ApVec(j) << endl;
		//}



		//arma::fcolvec ApVectemp = (Rmat.t()) * sqrtWinvNVec;
		//arma::fcolvec ApVec = ApVectemp * (ApVectemp.t()) * pVec - (1/Dvec) % pVec;
		//cout << "RmatIndex.n_cols " << RmatIndex.n_cols << endl;
		//cout << "RmatIndex.n_rows " << RmatIndex.n_rows << endl;
		//cout << "sqrtWinvNVec.n_elem " << sqrtWinvNVec.n_elem << endl;
		//arma::fcolvec ApVectemp = (RmatIndex.t()) * sqrtWinvNVec; 	
		//arma::fcolvec ApVec0 = ApVectemp * (ApVectemp.t()) * pVec + Dvec % pVec;
		//cout << "ApVec(10): " << ApVec(10) << endl;
		//cout << "pVec(10): " << pVec(10) << endl;
		//cout << "ApVec0(0): " << ApVec0(10) << endl;
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);
                float a = preA(0);
		//cout << "a: " << a << endl;
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
		z1Vec = minvVec % r1Vec;
		//arma::fvec z1Vec = r1Vec;
		arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
		//cout << "bet: " << bet << endl;	
		//cout << "Prebet.n_elem: " << Prebet.n_elem << endl;	
		//cout << "z1Vec(10): " << z1Vec(10) << endl;
		//cout << "r1Vec(10): " << r1Vec(10) << endl;
		zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
        }
        //if (iter >= maxiterPCG){
        //        cout << "pcg did not converge. You may increase maxiter number." << endl;
        //}
	//cout << "sumr2 is " << sumr2 << endl;
        //cout << "iter from getPCG1ofACinvAndVector " << iter << endl;
	//xVec = xVec *(1e+3);
        return(xVec);
}

// [[Rcpp::export]]
arma::fcolvec getProdRtb_Surv(arma::fcolvec& bVec, arma::fvec & RvecIndex){
	unsigned int kuniqtime = bVec.n_elem;
	arma::fcolvec bsumVec;
	arma::fcolvec Rtbvec;
	unsigned int m_N = geno.getNnomissing();
	Rtbvec.zeros(m_N);	
	bsumVec.zeros(kuniqtime);
	bsumVec(0) = bVec(0);
	for(unsigned int i = 1; i < kuniqtime; i++){
		bsumVec(i) = bsumVec(i-1) + bVec(i);	
	}
	int ktime;
        for(unsigned int j = 0; j < m_N; j++){
        	ktime = RvecIndex(j);
		Rtbvec(j) = bsumVec(ktime-1);
        }
	return(Rtbvec);
}

// [[Rcpp::export]]
arma::fcolvec getCrossprod_Surv_new(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & NWinv, arma::fvec & Dvec, unsigned int kuniqtime, int maxiterPCG, float tolPCG){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
 	//cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
 	//cout << "OKKKKK4" << endl;
 	//cout << "crossProdVec0(0) " << crossProdVec0(0) << endl;
        //WinvNRtG = (WinvNRt.t()) * bVec;
	//cout << "NWinv: " << endl;
        //for(size_t i=0; i< 10; i++){

          //      cout << NWinv(i) << " " << endl;
        //}

	//cout << "bVec: " << endl;
	//for(size_t i=0; i< 10; i++){

          //      cout << bVec(i) << " " << endl;
        //}


	arma::fcolvec NWinvbVec =  NWinv % bVec;
        //cout << "OKKKKK5" << endl;
	


	/*
	cout << NWinvbVec.n_elem << endl;
	cout << NWinvbVec(0) << endl;
	cout << Rmat.n_cols << endl;
	cout << Rmat.n_rows << endl;
	
	arma::fmat Rmatt = Rmat.t();
	cout << Rmatt.n_rows << endl;
	cout << Rmatt.n_cols << endl;

	cout << "NWinvbVec: " << endl;
	for(size_t i=0; i< 10; i++){

		cout << NWinvbVec(i) << " " << endl;
  	}



// 	arma::fcolvec RNWinvb0 = Rmatt * NWinvbVec;


//	cout << "RNWinvb(0) " << RNWinvb(0) << endl;
*/
	RNWinvb = getProdRb_Surv(NWinvbVec, RvecIndex, kuniqtime);
 //	RNWinvb0 = (Rmat.t()) * NWinvbVec;
//	arma::fcolvec RNWinvb1 =   Rmat.t() * (NWinv % bVec);

//	cout << "RNWinvb(0) " << RNWinvb(0) << endl;
//	//cout << "RNWinvb0(0) " << RNWinvb0(0) << endl;
//	cout << "RNWinvb1(0) " << RNWinvb1(0) << endl;
		
//	cout << "OKKKKK5" << endl;
	//arma::fcolvec RNWinvb;
	// cout << "OKKKKK6" << endl;
	//cout << "RNWinvb(0) is " << RNWinvb(0) << endl;
	//arma::fcolvec DRNWinvb = Dvec % RNWinvb;
	//cout << "DRNWinvb(0) is " << DRNWinvb(0) << endl;
	// cout << "OKKKKK7" << endl;
	arma::fcolvec AinvRNWinvb;
	  //for(size_t i=0; i< 5; i++){
	  //                cout << "i: " << i << " RNWinvb(i) " << RNWinvb(i) << endl;
	//		     }
/*
	  for(size_t i=0; i< 10; i++){
	                  cout << "i: " << i << " RvecIndex(i) " << RvecIndex(i) << endl;
			          }

	  for(size_t i=0; i< 10; i++){
	                  cout << "i: " << i << " sqrtWinvNVec(i) " << sqrtWinvNVec(i) << endl;
			          }

	
	  for(size_t i=0; i< 10; i++){
	                  cout << "i: " << i << " Dvec(i) " << Dvec(i) << endl;
			          }
	*/
	float pxnorm = arma::norm(RNWinvb);

	RNWinvb = RNWinvb/pxnorm;
	AinvRNWinvb = getPCG1ofACinvAndVector(RNWinvb, RvecIndex, sqrtWinvNVec, Dvec, maxiterPCG, tolPCG, wVec, tauVec);
	AinvRNWinvb = AinvRNWinvb * pxnorm;
	//for(size_t i=0; i< 5; i++){
	//	cout << "i: " << i << " AinvRNWinvb(i) " << AinvRNWinvb(i) << endl;
	//}
	//arma::fmat sqrtWinvNRtDt(wVec.n_elem, Dvec.n_elem);
	//arma::fmat Dmat = diagmat(-1/Dvec);
	//arma::fmat sqrtWinvNmat = diagmat(sqrtWinvNVec);
	//cout << "OKKKKK7c" << endl;
	//arma::fmat sqrtWinvNRt0 = sqrtWinvNmat * Rmat;
	//cout << "OKKKKK7d" << endl;
	//cout << "sqrtWinvNRt0.n_cols: " << sqrtWinvNRt0.n_cols << endl;
	//cout << "sqrtWinvNVec.n_elem: " << sqrtWinvNVec.n_elem << endl;
	//arma::fmat sqrtWinvNRt2 = (sqrtWinvNRt0.t()) * sqrtWinvNRt0;
	//cout << "OKKKKK7e" << endl;
	//arma::fmat A = sqrtWinvNRt2 + Dmat;
	//cout << "A(0,0) " << A(0,0) << endl;	
	//cout << "A(1,1) " << A(1,1) << endl;	
	//cout << "A(2,2) " << A(2,2) << endl;	

	//arma::fcolvec Adiag = getDiagofA(RvecIndex,sqrtWinvNVec,Dvec);
	//cout << "Adiag(0) " << Adiag(0)<< endl;
	//cout << "Adiag(1) " << Adiag(1)<< endl;
	//cout << "Adiag(1) " << Adiag(1)<< endl;

	//arma::fvec AinvRNWinvb0 = solve(A, RNWinvb);
	//cout << "AinvRNWinvb0(0) is " << AinvRNWinvb0(0) << endl;
	//cout << "AinvRNWinvb(0) is " << AinvRNWinvb(0) << endl;
	//cout << "AinvDRNWinvb(0) is " << AinvDRNWinvb(0) << endl;
	//cout << "AinvDRNWinvb0(0) is " << AinvDRNWinvb0(0) << endl;
	//cout << "OKKKKK7b" << endl;
		
	//arma::fcolvec DAinvDRNWinvb = Dvec % AinvDRNWinvb0;
	// cout << "OKKKKK8" << endl;
	// cout << "DAinvDRNWinvb(0) is " << DAinvDRNWinvb(0) << endl;
	//arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb0, RvecIndex);
	arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb, RvecIndex);
	//cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
	crossProdVec1 = NWinv % RtAinvDRNWinvb;
	//cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
 	//cout << "OKKKKK9" << endl;
        //RtAinvDRNWinvb = Rmat/(1e+10)  * AinvRNWinvb;
	//cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
        //crossProdVec1 = NWinv % (Rmat  * AinvRNWinvb);
	//cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
	//crossProdVec1 = getprodWinvNRttandVec(ACivWinvNRtG, RmatIndex, WinvN, kuniqtime);
 	//cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }else{
        	arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);

        	crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

        	return(crossProdVec);
	}
}



// [[Rcpp::export]]
arma::fvec getPCG1ofWminusUAndVector(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & NVec, arma::fvec & sqrtDVec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){

                   //  Start Timers
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = sqrtDVec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;
    //cout << "xVec: " << endl;	
    //xVec.print();
   // arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
	//cout << "rVec: " << endl;
	//rVec.print();
   arma::fvec rVec = bVec;

        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        //double wall1 = get_wall_time();
        //double cpu1  = get_cpu_time();
	//minvVec = diagofWminusUinv;
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
	//cout << "rVec1: " << endl;
	//dofWminusU.print();
	minvVec = 1/dofWminusU; 
	//cout << "minvVec: " << endl;
	//minvVec.print();
        zVec = minvVec % rVec;
	//cout << "rVec: " << endl;
	//rVec.print();
	//cout << "rVec2: " << endl;
		//zVec = rVec;
	//double wall2 = get_wall_time();
	//double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
  //      cout << "sumr2: " << sumr2 << endl;
	//cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
		//cout << "OKKKKKK" << endl;

		//arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
		//arma::fcolvec RWinN =  Rmat.t() * WinvN;
		//cout << "RWinN(0) is " << RWinN(0) << endl;


		//cout << "RWinNpVec(0) is " << RWinNpVec(0) << endl;	
                arma::fcolvec ApVec = getProdWminusUb_Surv(pVec, RvecIndex, NVec, sqrtDVec, wVec);
		//arma::fcolvec ApVec = getCrossprod_Surv_new(pVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
		//cout << "ApVec is " << ApVec(0) << endl;
		//cout << "OKKKKKK2" << endl;
		/*
		arma::fcolvec ApVec0;
 		arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
        	WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
	*/




                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        //if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        //}else{
        //        z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        //}

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
//        cout << "iter from getPCG1ofSigmaAndVector_WminusU " << iter << endl;
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}


// [[Rcpp::export]]
arma::fcolvec getCrossprod_Surv_new2(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex,  arma::fvec & NVec, arma::fvec & sqrtDVec, arma::fcolvec & diagofWminusUinv, unsigned int kuniqtime, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
	//arma::fcolvec sqrtDVec = arma::sqrt(Dvec);
	arma::fcolvec x0Vec(bVec.n_elem);
	x0Vec.zeros();
 	//cout << "OKKKKK2" << endl;
        crossProdVec0 = tauVec(0)*getPCG1ofWminusUAndVector(wVec,tauVec,bVec,RvecIndex,NVec,sqrtDVec,diagofWminusUinv,x0Vec,maxiterPCG, tolPCG, dofWminusU);

 	//cout << "OKKKKK3" << endl;
        if(tauVec(1) == 0){

                return(crossProdVec0);
        }else{
        	arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);

        	crossProdVec = crossProdVec0 +  tauVec(1)*crossProd1;

        	return(crossProdVec);
	}	
}



// [[Rcpp::export]]
arma::fcolvec getCrossprod_Surv_new2_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex,  arma::fvec & NVec, arma::fvec & sqrtDVec,  arma::fcolvec & diagofWminusUinv, unsigned int kuniqtime, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        //arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
 	//cout << "OKKKKK3" << endl;
	//arma::fcolvec sqrtDVec = arma::sqrt(Dvec);
	arma::fcolvec x0Vec(bVec.n_elem);
	x0Vec.zeros();
        crossProdVec0 = tauVec(0)*getPCG1ofWminusUAndVector(wVec,tauVec,bVec,RvecIndex,NVec,sqrtDVec,diagofWminusUinv,x0Vec,maxiterPCG, tolPCG, dofWminusU);

        if(tauVec(1) == 0){

                return(crossProdVec0);
        }else{
        	arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);

        	crossProdVec = crossProdVec0 +  tauVec(1)*crossProd1;

        	return(crossProdVec);
	}	
}


// [[Rcpp::export]]
arma::fcolvec getCrossprod_Surv_new_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & NWinv, arma::fvec & Dvec, unsigned int kuniqtime, int maxiterPCG, float tolPCG){
        arma::fcolvec crossProdVec;
        arma::fcolvec crossProdVec0;
        arma::fcolvec crossProdVec1;
        arma::fcolvec RNWinvb;
        //arma::fmat ACivWinvNRtG;
 	//cout << "OKKKKK3" << endl;
        crossProdVec0 = tauVec(0)*(bVec % (1/wVec));
 	//cout << "OKKKKK4" << endl;
 	//cout << "crossProdVec0(0) " << crossProdVec0(0) << endl;
        //WinvNRtG = (WinvNRt.t()) * bVec;
	//cout << "NWinv: " << endl;
        //for(size_t i=0; i< 10; i++){

          //      cout << NWinv(i) << " " << endl;
        //}

	//cout << "bVec: " << endl;
	//for(size_t i=0; i< 10; i++){

          //      cout << bVec(i) << " " << endl;
        //}


	arma::fcolvec NWinvbVec =  NWinv % bVec;
        //cout << "OKKKKK5" << endl;
	


	/*
	cout << NWinvbVec.n_elem << endl;
	cout << NWinvbVec(0) << endl;
	cout << Rmat.n_cols << endl;
	cout << Rmat.n_rows << endl;
	
	arma::fmat Rmatt = Rmat.t();
	cout << Rmatt.n_rows << endl;
	cout << Rmatt.n_cols << endl;

	cout << "NWinvbVec: " << endl;
	for(size_t i=0; i< 10; i++){

		cout << NWinvbVec(i) << " " << endl;
  	}



// 	arma::fcolvec RNWinvb0 = Rmatt * NWinvbVec;


//	cout << "RNWinvb(0) " << RNWinvb(0) << endl;
*/
	RNWinvb = getProdRb_Surv(NWinvbVec, RvecIndex, kuniqtime);
 //	RNWinvb0 = (Rmat.t()) * NWinvbVec;
//	arma::fcolvec RNWinvb1 =   Rmat.t() * (NWinv % bVec);

//	cout << "RNWinvb(0) " << RNWinvb(0) << endl;
//	//cout << "RNWinvb0(0) " << RNWinvb0(0) << endl;
//	cout << "RNWinvb1(0) " << RNWinvb1(0) << endl;
		
//	cout << "OKKKKK5" << endl;
	//arma::fcolvec RNWinvb;
	// cout << "OKKKKK6" << endl;
	//cout << "RNWinvb(0) is " << RNWinvb(0) << endl;
	//arma::fcolvec DRNWinvb = Dvec % RNWinvb;
	//cout << "DRNWinvb(0) is " << DRNWinvb(0) << endl;
	// cout << "OKKKKK7" << endl;
	arma::fcolvec AinvRNWinvb;
	/*  for(size_t i=0; i< 10; i++){
	                  cout << "i: " << i << " RNWinvb(i) " << RNWinvb(i) << endl;
			          }

	  for(size_t i=0; i< 10; i++){
	                  cout << "i: " << i << " RvecIndex(i) " << RvecIndex(i) << endl;
			          }

	  for(size_t i=0; i< 10; i++){
	                  cout << "i: " << i << " sqrtWinvNVec(i) " << sqrtWinvNVec(i) << endl;
			          }

	
	  for(size_t i=0; i< 10; i++){
	                  cout << "i: " << i << " Dvec(i) " << Dvec(i) << endl;
			          }
	*/
	float pxnorm = arma::norm(RNWinvb);
	RNWinvb = RNWinvb/pxnorm;
	AinvRNWinvb = getPCG1ofACinvAndVector(RNWinvb, RvecIndex, sqrtWinvNVec, Dvec, maxiterPCG, tolPCG, wVec, tauVec);
	AinvRNWinvb = AinvRNWinvb * pxnorm;
	//for(size_t i=0; i< 10; i++){
	//	cout << "i: " << i << " AinvRNWinvb(i) " << AinvRNWinvb(i) << endl;
	//}
	//arma::fmat sqrtWinvNRtDt(wVec.n_elem, Dvec.n_elem);
	//arma::fmat Dmat = diagmat(-1/Dvec);
	//arma::fmat sqrtWinvNmat = diagmat(sqrtWinvNVec);
	//cout << "OKKKKK7c" << endl;
	//arma::fmat sqrtWinvNRt0 = sqrtWinvNmat * Rmat;
	//cout << "OKKKKK7d" << endl;
	//cout << "sqrtWinvNRt0.n_cols: " << sqrtWinvNRt0.n_cols << endl;
	//cout << "sqrtWinvNVec.n_elem: " << sqrtWinvNVec.n_elem << endl;
	//arma::fmat sqrtWinvNRt2 = (sqrtWinvNRt0.t()) * sqrtWinvNRt0;
	//cout << "OKKKKK7e" << endl;
	//arma::fmat A = sqrtWinvNRt2 + Dmat;
	//cout << "A(0,0) " << A(0,0) << endl;	
	//cout << "A(1,1) " << A(1,1) << endl;	
	//cout << "A(2,2) " << A(2,2) << endl;	

	//arma::fcolvec Adiag = getDiagofA(RvecIndex,sqrtWinvNVec,Dvec);
	//cout << "Adiag(0) " << Adiag(0)<< endl;
	//cout << "Adiag(1) " << Adiag(1)<< endl;
	//cout << "Adiag(1) " << Adiag(1)<< endl;

	//arma::fvec AinvRNWinvb0 = solve(A, RNWinvb);
	//cout << "AinvRNWinvb0(0) is " << AinvRNWinvb0(0) << endl;
	//cout << "AinvRNWinvb(0) is " << AinvRNWinvb(0) << endl;
	//cout << "AinvDRNWinvb(0) is " << AinvDRNWinvb(0) << endl;
	//cout << "AinvDRNWinvb0(0) is " << AinvDRNWinvb0(0) << endl;
	//cout << "OKKKKK7b" << endl;
		
	//arma::fcolvec DAinvDRNWinvb = Dvec % AinvDRNWinvb0;
	// cout << "OKKKKK8" << endl;
	// cout << "DAinvDRNWinvb(0) is " << DAinvDRNWinvb(0) << endl;
	//arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb0, RvecIndex);
	arma::fcolvec RtAinvDRNWinvb = getProdRtb_Surv(AinvRNWinvb, RvecIndex);
	//cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
	crossProdVec1 = NWinv % RtAinvDRNWinvb;
	//cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
 	//cout << "OKKKKK9" << endl;
        //RtAinvDRNWinvb = Rmat/(1e+10)  * AinvRNWinvb;
	//cout << "RtAinvDRNWinvb is " << RtAinvDRNWinvb(0) << endl;
        //crossProdVec1 = NWinv % (Rmat  * AinvRNWinvb);
	//cout << "crossProdVec1(0) is " << crossProdVec1(0) << endl;
	//crossProdVec1 = getprodWinvNRttandVec(ACivWinvNRtG, RmatIndex, WinvN, kuniqtime);
 	//cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);

        crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;

        return(crossProdVec);
}


// [[Rcpp::export]]
arma::fcolvec getCrossprod_LOCO(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1  = getCrossprodMatAndKin_LOCO(bVec);
        crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;

        return(crossProdVec);
}




double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}



// [[Rcpp::export]]
arma::sp_mat gen_sp_GRM() {
    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}


// [[Rcpp::export]]
arma::sp_mat gen_sp_Sigma(arma::fvec& wVec,  arma::fvec& tauVec){
   arma::fvec dtVec = (1/wVec) * (tauVec(0));
//   dtVec.print();
   arma::vec valueVecNew = valueVec * tauVec(1);

   int nnonzero = valueVec.n_elem;
   for(size_t i=0; i< nnonzero; i++){
     if(locationMat(0,i) == locationMat(1,i)){
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
       valueVecNew(i) = valueVecNew(i) + dtVec(locationMat(0,i));
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
	if(valueVecNew(i) < 1e-4){
  			valueVecNew(i) = 1e-4 ;
  		}


     }
   }

    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVecNew, dimNum, dimNum);
//    std::cout << "result.n_rows " << result.n_rows << std::endl;
//    std::cout << "result.n_cols " << result.n_cols << std::endl;
    //result.print();
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}



// [[Rcpp::export]]
arma::vec gen_spsolve_v3(arma::vec & yvec){
    // sparse x sparse -> sparse
    //arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    //arma::vec y = arma::linspace<arma::vec>(0, 5, dimNum);
    arma::sp_mat result = gen_sp_GRM();

    std::cout << "yvec.n_elem: " << yvec.n_elem << std::endl;
    std::cout << "result.n_rows: " << result.n_rows << std::endl;
    std::cout << "result.n_cols: " << result.n_cols << std::endl;
    arma::vec x = arma::spsolve(result, yvec);

    return x;
}

// [[Rcpp::export]]
arma::fvec gen_spsolve_v4(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec & yvec){

//    double wall0in = get_wall_time();
// double cpu0in  = get_cpu_time();

    arma::vec yvec2 = arma::conv_to<arma::vec>::from(yvec);
//double wall1in = get_wall_time();
// double cpu1in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall1in - wall0in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu1in - cpu0in  << endl;

    arma::sp_mat result = gen_sp_Sigma(wVec, tauVec);

//double wall2in = get_wall_time();
// double cpu2in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall2in - wall1in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu2in - cpu1in  << endl;

//    std::cout << "yvec.n_elem: " << yvec.n_elem << std::endl;
//    std::cout << "yvec2.n_elem: " << yvec2.n_elem << std::endl;
//    std::cout << "result.n_rows: " << result.n_rows << std::endl;
//    std::cout << "result.n_cols: " << result.n_cols << std::endl;
    arma::vec x = arma::spsolve(result, yvec2);

//double wall3in = get_wall_time();
// double cpu3in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;


    arma::fvec z = arma::conv_to<arma::fvec>::from(x);

//double wall4in = get_wall_time();
// double cpu4in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall4in - wall3in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu4in - cpu3in  << endl;


    return z;
}


//bool isUsePrecondM = false;
//bool isUseSparseSigmaforInitTau = false;



// [[Rcpp::export]]
void setisUsePrecondM(bool isUseSparseSigmaforPCG){
	isUsePrecondM = isUseSparseSigmaforPCG;
}

// [[Rcpp::export]]
void setisUseSparseSigmaforInitTau(bool isUseSparseSigmaforInitTau0){
	isUseSparseSigmaforInitTau = isUseSparseSigmaforInitTau0;
}

// [[Rcpp::export]]
void setisUseSparseSigmaforNullModelFitting(bool isUseSparseSigmaforModelFitting0){
	isUseSparseSigmaforModelFitting = isUseSparseSigmaforModelFitting0;
}
//Modified on 11-28-2018 to allow for a preconditioner for CG (the sparse Sigma)                                                                                                                                     //Sigma = tau[1] * diag(1/W) + tau[2] * kins
//This function needs the function getDiagOfSigma and function getCrossprod



// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG){
       
                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
	int Nnomissing = geno.getNnomissing();
	arma::fvec xVec(Nnomissing);
        xVec.zeros();

if(isUseSparseSigmaforInitTau){
	cout << "use sparse kinship to estimate initial tau " <<  endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else if(isUseSparseSigmaforModelFitting){
	cout << "use sparse kinship to fit the model " << endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else{
        arma::fvec rVec = bVec;
        //cout << "HELLOa: "  << endl;
        arma::fvec r1Vec;
        //cout << "HELLOb: "  << endl;
        //int Nnomissing = geno.getNnomissing();
        //cout << "HELL1: "  << endl;

        arma::fvec crossProdVec(Nnomissing);
        //cout << "HELL2: "  << endl;

        //arma::SpMat<float> precondM = sparseGRMinC;
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);

     double wall1 = get_wall_time();
       double cpu1  = get_cpu_time();

//    cout << "Wall Time 1= " << wall1 - wall0 << endl;
//    cout << "CPU Time 1 = " << cpu1  - cpu0  << endl;
        if (!isUsePrecondM){
	//arma::fvec diagVec = getDiagOfSigma(wVec, tauVec);
        //cout<<"diagVec:" << endl;
        //diagVec.print();

                minvVec = 1/getDiagOfSigma(wVec, tauVec);
	//cout<<"minvVec:" << endl;
	//minvVec.print();
                zVec = minvVec % rVec;
        }else{


//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;




		zVec = gen_spsolve_v4(wVec, tauVec, rVec);
                //sparseGRMinC = (sparseGRMinC) * (tauVec(1));
                //arma::fvec dtVec = (1/wVec) * (tauVec(0));
                //(sparseGRMinC).diag() = (sparseGRMinC).diag() + dtVec;
                //zVec = gen_spsolve_v4(sparseGRMinC, rVec) ;
        }
 double wall2 = get_wall_time();
 double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);

/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
	//std::cout << "pVec: " << std::endl;
	//pVec.print();
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;



        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
		//cout << "ApVec: " << endl;
		//ApVec.print();
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
        }else{
		z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                        //std::cout << "sumr2: " << sumr2 << std::endl;
                        //std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}


// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_Surv(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
        int Nnomissing = geno.getNnomissing();
        arma::fvec xVec(Nnomissing);
        //xVec.zeros();
	xVec = x0Vec;
if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else if(isUseSparseSigmaforModelFitting){
	cout << "use sparse kinship to fit the model " << endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else{
        //arma::fvec rVec = bVec;
	arma::fvec rVec = bVec - getCrossprod_Surv(xVec, wVec, tauVec, WinvNRt, ACinv); 

        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);

       double wall1 = get_wall_time();
       double cpu1  = get_cpu_time();

        if (!isUsePrecondM){
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
		minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
                zVec = minvVec % rVec;
		//zVec = rVec;
	}else{


                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }


        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
	//cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
	//cout << "OKKKKKK2" << endl;
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
		//z1Vec = r1Vec;
	}else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "sumr2: " << sumr2 << std::endl;
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}


// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_Surv_LOCO(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec,int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
        int Nnomissing = geno.getNnomissing();
        arma::fvec xVec(Nnomissing);
        //xVec.zeros();
	xVec = x0Vec;

if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else if(isUseSparseSigmaforModelFitting){
	cout << "use sparse kinship to fit the model " << endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec);
}else{
        arma::fvec rVec = bVec -  getCrossprod_Surv_LOCO(xVec, wVec, tauVec, WinvNRt, ACinv);
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);

       double wall1 = get_wall_time();
       double cpu1  = get_cpu_time();

        if (!isUsePrecondM){
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
		minvVec = 1/getDiagOfSigma_surv_LOCO(diagofWminusUinv, tauVec);
                //zVec = minvVec % rVec;
		zVec = rVec;
	}else{


                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }


        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
	//cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                arma::fcolvec ApVec = getCrossprod_Surv_LOCO(pVec, wVec, tauVec, WinvNRt, ACinv);
	//cout << "OKKKKKK2" << endl;
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                //z1Vec = minvVec % r1Vec;
		z1Vec = r1Vec;
	}else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "sumr2: " << sumr2 << std::endl;
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}




// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_Surv_new(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = Dvec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;
    //cout << "xVec: " << endl;	
    //xVec.print();


    if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else if(isUseSparseSigmaforModelFitting){
	cout << "use sparse kinship to fit the model " << endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else{
        arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
	//cout << "rVec: " << rVec[0] << endl;
	//cout << "bVec: " << bVec[0] << endl;
	//arma::fvec csn = getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
	//cout << "csn: " << csn[0] << endl;
	//rVec.print();


        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
		minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                zVec = minvVec % rVec;
		//zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
	double wall2 = get_wall_time();
	double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
	//cout << "sumr2: " << sumr2 << endl;
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
	//cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
		//cout << "OKKKKKK" << endl;

		//arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
		//arma::fcolvec RWinN =  Rmat.t() * WinvN;
		//cout << "RWinN(0) is " << RWinN(0) << endl;


		//cout << "RWinNpVec(0) is " << RWinNpVec(0) << endl;

		arma::fcolvec ApVec = getCrossprod_Surv_new(pVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
		//cout << "ApVec is " << ApVec(0) << endl;
		//cout << "OKKKKKK2" << endl;
		/*
		arma::fcolvec ApVec0;
 		arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
        	WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
	*/




                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}



// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_Surv_new2(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & NVec, arma::fvec & sqrtDvec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = sqrtDvec.n_elem;
    arma::fvec xVec(Nnomissing);
    xVec.zeros();
    //xVec = x0Vec;
    cout << "xVec: " << endl;	
    //xVec.print();


    if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else if(isUseSparseSigmaforModelFitting){
	cout << "use sparse kinship to fit the model " << endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else{
        //arma::fvec rVec = bVec; 
        arma::fvec rVec = bVec - getCrossprod_Surv_new2(xVec,  wVec, tauVec, RvecIndex,NVec, sqrtDvec,diagofWminusUinv, kuniqtime, maxiterPCG, tolPCG, dofWminusU);
        //arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
	cout << "rVec: " << endl;
	//rVec.print();


        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
		minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                zVec = minvVec % rVec;
		//zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
	double wall2 = get_wall_time();
	double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
	//cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
		//cout << "OKKKKKK" << endl;

		//arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
		//arma::fcolvec RWinN =  Rmat.t() * WinvN;
		//cout << "RWinN(0) is " << RWinN(0) << endl;


		cout << "pVec(0) is " << pVec(0) << endl;
		arma::fcolvec ApVec = getCrossprod_Surv_new2(pVec, wVec, tauVec, RvecIndex,NVec, sqrtDvec, diagofWminusUinv, kuniqtime, maxiterPCG, tolPCG, dofWminusU);
		cout << "ApVec is " << ApVec(0) << endl;


		//cout << "OKKKKKK2" << endl;
		/*
		arma::fcolvec ApVec0;
 		arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
        	WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
	*/
		arma::fvec pAp = pVec.t() * ApVec;
		if(pAp(0) == 0){
			return(xVec);
		}
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/

                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
       std::cout << "sumr2: " << sumr2 << std::endl;
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}


// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_Surv_LOCO_new2(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & NVec, arma::fvec & sqrtDvec, arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = sqrtDvec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;
    //cout << "xVec: " << endl;	
    //xVec.print();


    if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else if(isUseSparseSigmaforModelFitting){
	cout << "use sparse kinship to fit the model " << endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else{
        //arma::fvec rVec = bVec - getCrossprod_Surv_new2(xVec,  wVec, tauVec, RvecIndex,NVec, Dvec, kuniqtime, maxiterPCG, tolPCG);
        arma::fvec rVec = bVec - getCrossprod_Surv_new2_LOCO(xVec,  wVec, tauVec, RvecIndex,NVec, sqrtDvec, diagofWminusUinv, kuniqtime, maxiterPCG, tolPCG, dofWminusU);
        //arma::fvec rVec = bVec - getCrossprod_Surv_new(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
	//cout << "rVec: " << endl;
	//rVec.print();


        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
		minvVec = 1/getDiagOfSigma_surv_LOCO(diagofWminusUinv, tauVec);
                //minvVec = 1/getDiagOfSigma(wVec, tauVec);
                zVec = minvVec % rVec;
		//zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
	double wall2 = get_wall_time();
	double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
	//cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
		//cout << "OKKKKKK" << endl;

		//arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
		//arma::fcolvec RWinN =  Rmat.t() * WinvN;
		//cout << "RWinN(0) is " << RWinN(0) << endl;


		//cout << "RWinNpVec(0) is " << RWinNpVec(0) << endl;
		arma::fcolvec ApVec = getCrossprod_Surv_new2_LOCO(pVec, wVec, tauVec, RvecIndex,NVec, sqrtDvec, diagofWminusUinv,  kuniqtime, maxiterPCG, tolPCG, dofWminusU);
		//cout << "ApVec is " << ApVec(0) << endl;
		//cout << "OKKKKKK2" << endl;
		/*
		arma::fcolvec ApVec0;
 		arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
        	WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
	*/




                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}




// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_Surv_new_LOCO(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec,  arma::fvec & diagofWminusUinv, arma::fvec & x0Vec, int maxiterPCG, float tolPCG){

                   //  Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    unsigned int kuniqtime = Dvec.n_elem;
    arma::fvec xVec(Nnomissing);
    //xVec.zeros();
    xVec = x0Vec;

    if(isUseSparseSigmaforInitTau){
        cout << "use sparse kinship to estimate initial tau " <<  endl;
        xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else if(isUseSparseSigmaforModelFitting){
	cout << "use sparse kinship to fit the model " << endl;
	xVec = gen_spsolve_v4(wVec, tauVec, bVec); //to update
    }else{
        arma::fvec rVec = bVec - getCrossprod_Surv_new_LOCO(xVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        double wall1 = get_wall_time();
        double cpu1  = get_cpu_time();
        if (!isUsePrecondM){
                //minvVec = 1/getDiagOfSigma_surv(diagofWminusUinv, tauVec);
		minvVec = 1/getDiagOfSigma_surv_LOCO(diagofWminusUinv, tauVec);
                zVec = minvVec % rVec;
		//zVec = rVec;
        }else{
                zVec = gen_spsolve_v4(wVec, tauVec, rVec);
        }
	double wall2 = get_wall_time();
	double cpu2  = get_cpu_time();
// cout << "Wall Time 2 = " << wall2 - wall1 << endl;
// cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;


//      cout << "HELL3: "  << endl;
//      for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
        float sumr2 = sum(rVec % rVec);
/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
	//cout << "OKKKKKK" << endl;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                //arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
                //arma::fcolvec ApVec = getCrossprod_Surv(pVec, wVec, tauVec, WinvNRt, ACinv);
		//cout << "OKKKKKK" << endl;

		//arma::fcolvec RWinNpVec =  Rmat.t() * (WinvN % pVec);
		//arma::fcolvec RWinN =  Rmat.t() * WinvN;
		//cout << "RWinN(0) is " << RWinN(0) << endl;


		//cout << "RWinNpVec(0) is " << RWinNpVec(0) << endl;

		arma::fcolvec ApVec = getCrossprod_Surv_new_LOCO(pVec, wVec, tauVec, RvecIndex, sqrtWinvNVec,WinvN,Dvec, kuniqtime, maxiterPCG, tolPCG);
		//cout << "ApVec is " << ApVec(0) << endl;
		//cout << "OKKKKKK2" << endl;
		/*
		arma::fcolvec ApVec0;
 		arma::fcolvec crossProdVec0 = tauVec(0)*(pVec % (1/wVec));
        	WinvNRtG = (WinvNRt.t()) * bVec;
        //cout << "OKKKKK5" << endl;
        ACivWinvNRtG = ACinv * WinvNRtG;
        //cout << "OKKKKK6" << endl;
        crossProdVec1 = WinvNRt * ACivWinvNRtG;
        //cout << "OKKKKK7" << endl;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = crossProdVec0 - tauVec(0)*crossProdVec1;

                return(crossProdVec);
        }
        arma::fvec crossProd1  = getCrossprodMatAndKin(bVec);
        crossProdVec = crossProdVec0 + tauVec(0)*crossProdVec1 + tauVec(1)*crossProd1;
	*/




                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
                //z1Vec = r1Vec;
        }else{
                z1Vec = gen_spsolve_v4(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }
       //std::cout << "sumr2: " << sumr2 << std::endl;

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}



//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
//This function needs the function getDiagOfSigma and function getCrossprod
// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_old(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG){
	           //  Start Timers
//    double wall0 = get_wall_time();
//    double cpu0  = get_cpu_time();

//	 std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	//cout << "HELLO: "  << endl;
	//cout << "HELL2: "  << endl;
  	arma::fvec rVec = bVec;
	//cout << "HELLOa: "  << endl;
  	arma::fvec r1Vec;
	//cout << "HELLOb: "  << endl;
  	int Nnomissing = geno.getNnomissing();
	//cout << "HELL1: "  << endl;

  	arma::fvec crossProdVec(Nnomissing);
	//cout << "HELL2: "  << endl;

  	arma::fvec minvVec = 1/getDiagOfSigma(wVec, tauVec);
//	cout << "HELL3: "  << endl;
//	for(int i = 0; i < 10; i++){
//                cout << "full set minvVec[i]: " << minvVec[i] << endl;
//        }
  	float sumr2 = sum(rVec % rVec);
/*
	if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
  	arma::fvec zVec = minvVec % rVec;
  	arma::fvec z1Vec;
 	arma::fvec pVec = zVec;
	/*
        if(bVec[0] == 1 && bVec[2] == 1){
	for(int i = 0; i < 10; i++){ 
		cout << "pVec[i]: " << i << " " << pVec[i] << endl;
 	}
	}
*/
  	arma::fvec xVec(Nnomissing);
  	xVec.zeros();
  
  	int iter = 0;
  	while (sumr2 > tolPCG && iter < maxiterPCG) {
    		iter = iter + 1;
    		arma::fcolvec ApVec = getCrossprod(pVec, wVec, tauVec);
    		arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

    		float a = preA(0);
	
/*	     if(bVec[0] == 1 && bVec[2] == 1){
			cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
        		for(int i = 0; i < 10; i++){
			
                		cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                		cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                		cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                		cout << "rVec[i]: " << i << " " << rVec[i] << endl;
        		}
        	    }   
*/	
 
    		xVec = xVec + a * pVec;
/*
		if(bVec[0] == 1 && bVec[2] == 1){
        		for(int i = 0; i < 10; i++){
                		cout << "xVec[i]: " << i << " " << xVec[i] << endl;
        		}
        	}   

*/

 
    		r1Vec = rVec - a * ApVec;
/*
		if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
    		z1Vec = minvVec % r1Vec;
    		arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
    		float bet = Prebet(0);
    		pVec = z1Vec+ bet*pVec;
    		zVec = z1Vec;
    		rVec = r1Vec;
    
    		sumr2 = sum(rVec % rVec);
/*
		if(bVec[0] == 1 && bVec[2] == 1){
			std::cout << "sumr2: " << sumr2 << std::endl;
			std::cout << "tolPCG: " << tolPCG << std::endl;
		}
*/
  	}
  
  	if (iter >= maxiterPCG){
    		cout << "pcg did not converge. You may increase maxiter number." << endl;
     
  	}
  	cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;

//        double wall1 = get_wall_time();
//    double cpu1  = get_cpu_time();

//    cout << "Wall Time = " << wall1 - wall0 << endl;
//    cout << "CPU Time  = " << cpu1  - cpu0  << endl;
	
//	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
  	return(xVec);
}



//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
//This function needs the function getDiagOfSigma and function getCrossprod
// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_LOCO(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG){
  	arma::fvec rVec = bVec;
  	arma::fvec r1Vec;
  	int Nnomissing = geno.getNnomissing();

  	arma::fvec crossProdVec(Nnomissing);
  	arma::fvec minvVec = 1/getDiagOfSigma_LOCO(wVec, tauVec);
        //for(int i = 0; i < 10; i++){
	//	cout << "minvVec[i]: " << minvVec[i] << endl;
        //} 
 
  	float sumr2 = sum(rVec % rVec);

  	arma::fvec zVec = minvVec % rVec;
  	arma::fvec z1Vec;
 	arma::fvec pVec = zVec;
  
  	arma::fvec xVec(Nnomissing);
  	xVec.zeros();
  
  	int iter = 0;
  	while (sumr2 > tolPCG && iter < maxiterPCG) {
    		iter = iter + 1;
    		arma::fcolvec ApVec = getCrossprod_LOCO(pVec, wVec, tauVec);
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

/*
// [[Rcpp::export]] 
void setChromosomeIndicesforLOCO(vector<int> chromosomeStartIndexVec, vector<int> chromosomeEndIndexVec, vector<int> chromosomeVecVec){
  LOCO = true;
  chromosomeStartIndex = chromosomeStartIndexVec;
  chromosomeEndIndexV = chromosomeEndIndexVec;
  chromosomeVec = chromosomeVecVec;
}
*/

// [[Rcpp::export]]
void setStartEndIndex(int startIndex, int endIndex){
  geno.startIndex = startIndex;
  geno.endIndex = endIndex;
  geno.Msub = geno.M - (endIndex - startIndex + 1);
}

//This function calculates the coefficients of variation for mean of a vector
// [[Rcpp::export]]
float calCV(arma::fvec& xVec){
  int veclen = xVec.n_elem;
  float vecMean = arma::mean(xVec);
  float vecSd = arma::stddev(xVec);
  float vecCV = (vecSd/vecMean)/veclen;
  return(vecCV);
}

// [[Rcpp::export]]
float GetTrace(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& cov1, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){
  set_seed(200);
  int Nnomissing = geno.getNnomissing();
  arma::fmat Sigma_iXt = Sigma_iX.t();
  arma::fvec Sigma_iu;  
  arma::fcolvec Pu;
  arma::fvec Au;
  arma::fvec uVec;

  int nrunStart = 0;
  int nrunEnd = nrun;
  float traceCV = traceCVcutoff + 0.1;
  arma::fvec tempVec(nrun);
  tempVec.zeros();

  while(traceCV > traceCVcutoff){     
    //arma::fvec tempVec(nrun);
    //tempVec.zeros();
    for(int i = nrunStart; i < nrunEnd; i++){
      Rcpp::NumericVector uVec0;
      uVec0 = nb(Nnomissing);
      uVec = as<arma::fvec>(uVec0);
      uVec = uVec*2 - 1;
      Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG);
      Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
      Au = getCrossprodMatAndKin(uVec);
      tempVec(i) = dot(Au, Pu);
      Au.clear();
      Pu.clear();
      Sigma_iu.clear();
      uVec.clear();
    }
    traceCV = calCV(tempVec);
    if(traceCV > traceCVcutoff){
      nrunStart = nrunEnd;
      nrunEnd = nrunEnd + 10;
      tempVec.resize(nrunEnd); 
      cout << "CV for trace random estimator using "<< nrun << " runs is " << traceCV <<  " > " << traceCVcutoff << endl;
      cout << "try " << nrunEnd << " runs" << endl;      
    }
  }

  float tra = arma::mean(tempVec);
  tempVec.clear();
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
  	//arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
	//arma::fmat XSX = Xmatt * Sigma_iX;
	//std::cout << "Sigma_iX: " << std::endl;
	//Sigma_iX.print();
	
	//std::cout << "Xmatt * Sigma_iX: " << std::endl;
	//XSX.print();

	arma::fmat cov;
	try {
	  //cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
	  cov = arma::inv_sympd(Xmatt * Sigma_iX);
	} catch (const std::exception& e) {
	  //cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
	  cov = arma::pinv(Xmatt * Sigma_iX);
	  cout << "inv_sympd failed, inverted with pinv" << endl;
	}

 	arma::fmat Sigma_iXt = Sigma_iX.t();
  	arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
  	arma::fvec alpha = cov * SigmaiXtY;

  	arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
  	return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
}


// [[Rcpp::export]]
Rcpp::List getCoefficients_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG){

        int Nnomissing = geno.getNnomissing();
        arma::fvec Sigma_iY;

        Sigma_iY = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Yvec, maxiterPCG, tolPCG);
        int colNumX = Xmat.n_cols;
        arma::fmat Sigma_iX(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);

                Sigma_iX.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

        }

        arma::fmat Xmatt = Xmat.t();
	arma::fmat cov;
        try {
          cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
        } catch (const std::exception& e) {
          cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
          cout << "inv_sympd failed, inverted with pinv" << endl;
        }
        //arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
        arma::fmat Sigma_iXt = Sigma_iX.t();
        arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
        arma::fvec alpha = cov * SigmaiXtY;

        arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
        return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
}


// [[Rcpp::export]]
Rcpp::List getCoefficients_q_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG){

        int Nnomissing = geno.getNnomissing();
        arma::fvec Sigma_iY;

        Sigma_iY = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Yvec, maxiterPCG, tolPCG);
        int colNumX = Xmat.n_cols;
        arma::fmat Sigma_iX(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);

                Sigma_iX.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

        }

        arma::fmat Xmatt = Xmat.t();
	arma::fmat cov;
        try {
          cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
        } catch (const std::exception& e) {
          cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
          cout << "inv_sympd failed, inverted with pinv" << endl;
        }

        //arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
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
int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){
        //std::cout << "getAIScore starts" << std::endl;
	arma::fmat Sigma_iXt = Sigma_iX.t();

  	arma::fvec PY1 = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));
  	arma::fvec APY = getCrossprodMatAndKin(PY1);
  	float YPAPY = dot(PY1, APY);

  	float Trace = GetTrace(Sigma_iX, Xmat, wVec, tauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);
  	arma::fvec PAPY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG, tolPCG);
  	arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
  	float AI = dot(APY, PAPY);
        //std::cout << "getAIScore ends" << std::endl;

  	return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY1, Named("AI") = AI);
}

//Rcpp::List fitglmmaiRPCG(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec,
// Modified by SLEE, 04/16/2017
// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore
// [[Rcpp::export]]
Rcpp::List fitglmmaiRPCG(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec &tauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){
        //tauVec(1) = 0.6577277;
  	//std::cout << "fitglmmaiRPCG starts" << std::endl;
        Rcpp::List re = getAIScore(Yvec, Xmat,wVec,  tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);
  	float YPAPY = re["YPAPY"];
  	float Trace = re["Trace"];
  	float score1 = YPAPY - Trace;
  	float AI1 = re["AI"];
  	float Dtau = score1/AI1;
  	arma::fvec tau0 = tauVec;
  	tauVec(1) = tau0(1) + Dtau;
  	//std::cout << "fitglmmaiRPCG ends" << std::endl;

	//std::cout << "AI1: " << AI1 << std::endl;
	//std::cout << "score1: " << score1 << std::endl;
	//std::cout << "Dtau: " << Dtau << std::endl;


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
        //tauVec(1) = 0.6577277;
  	return List::create(Named("tau") = tauVec);
}



/*add for SPA by Wei 04222017*/
// [[Rcpp::export]]
arma::fmat getSigma_X(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG){


  	int Nnomissing = Xmat.n_rows;
  	int colNumX = Xmat.n_cols;

  	//cout << colNumX << endl;
  	//cout << size(wVec) << endl;
  	//cout << size(tauVec) << endl;


  	arma::fmat Sigma_iX1(Nnomissing,colNumX);
  	arma::fvec XmatVecTemp;

  	for(int i = 0; i < colNumX; i++){
    		XmatVecTemp = Xmat.col(i);
    		Sigma_iX1.col(i) = getPCG1ofSigmaAndVector(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);
  	}
  	return(Sigma_iX1);
}



// [[Rcpp::export]]
arma::fmat getSigma_X_Surv(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv,  arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;

        //cout << colNumX << endl;
        //cout << size(wVec) << endl;
        //cout << size(tauVec) << endl;
        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
        arma::fvec x0Vec;


        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
		//x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
		x0Vec = wVec % XmatVecTemp - sqrtDRN.t() * sqrtDRN * XmatVecTemp;
		//x0Vec.print();

		if(tauVec(1) != 0){

			Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv(wVec, tauVec, XmatVecTemp, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
		}else{
			Sigma_iX1.col(i) = x0Vec;
		}

        }
        return(Sigma_iX1);
}


// [[Rcpp::export]]
arma::fmat getSigma_X_Surv_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv,  arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;

        //cout << colNumX << endl;
        //cout << size(wVec) << endl;
        //cout << size(tauVec) << endl;
        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
	arma::fvec x0Vec;

	for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
                //x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
                x0Vec = wVec % XmatVecTemp - sqrtDRN.t() * sqrtDRN * XmatVecTemp;

                if(tauVec(1) != 0){

                        Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_LOCO(wVec, tauVec, XmatVecTemp, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
                }else{
                        Sigma_iX1.col(i) = x0Vec;
                }

        }


        return(Sigma_iX1);
}


// [[Rcpp::export]]
arma::fmat getSigma_X_Surv_new(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

	//arma::fvec XmatVecTemp1 = Xmat.col(0);
	//arma::fvec Sigma_iX1temp = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, XmatVecTemp1, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, maxiterPCG, tolPCG);
	//cout << "XmatVecTemp1, 1st column in Xmat" << endl;	

	arma::fvec x0Vec;
	arma::fvec sqrtDVec = arma::sqrt(Dvec);
        for(int i = 0; i < colNumX; i++){
		//if (i == 0){
		//	cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
		//}
                XmatVecTemp = Xmat.col(i);
	//	cout << "i is " << i << endl;
	//	cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
		x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
	
		if(tauVec(1) != 0){
                	Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, XmatVecTemp, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
		}else{
			Sigma_iX1.col(i) = x0Vec;
	//		Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
		}
	//	cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;
		
        }
        return(Sigma_iX1);
}


// [[Rcpp::export]]
arma::fmat getSigma_X_Surv_new_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
	  arma::fvec x0Vec;
        //arma::fvec XmatVecTemp1 = Xmat.col(0);
        //arma::fvec Sigma_iX1temp = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, XmatVecTemp1, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG);
        //cout << "XmatVecTemp1, 1st column in Xmat" << endl;
	arma::fvec sqrtDVec = arma::sqrt(Dvec);

        for(int i = 0; i < colNumX; i++){
                //if (i == 0){
                //      cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
                //}
                XmatVecTemp = Xmat.col(i);
		x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
        //      cout << "i is " << i << endl;
        //      cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
		if(tauVec(1)!=0){
                	Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new_LOCO(wVec, tauVec, XmatVecTemp, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
		}else{
			//Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
			Sigma_iX1.col(i) = x0Vec;
		}
        //      cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;

        }
        return(Sigma_iX1);
}



// [[Rcpp::export]]
arma::fmat  getSigma_X_Surv_new2(arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& Xmat,  arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

	arma::fvec x0Vec;
	arma::fvec sqrtDVec = arma::sqrt(Dvec);
        for(int i = 0; i < colNumX; i++){
		if (i == 0){
			cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
		}
                XmatVecTemp = Xmat.col(i);
	//	cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
		x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
		cout << "i is " << i << endl;
	
		if(tauVec(1) != 0){
                	Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new2(wVec, tauVec, XmatVecTemp, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
		}else{
                	Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_new2(wVec, tauVec, XmatVecTemp, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
			//Sigma_iX1.col(i) = x0Vec;
	//		Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
		}
	//	cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;
		
        }
        return(Sigma_iX1);
}


// [[Rcpp::export]]
arma::fmat  getSigma_X_Surv_new2_LOCO(arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& Xmat,  arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG,  arma::fvec & dofWminusU){

        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
	  arma::fvec x0Vec;
	arma::fvec sqrtDVec = arma::sqrt(Dvec);

        for(int i = 0; i < colNumX; i++){
                //if (i == 0){
                //      cout << "XmatVecTemp1, 1st column in Xmat b" << endl;
                //}
                XmatVecTemp = Xmat.col(i);
		x0Vec = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
        //      cout << "i is " << i << endl;
        //      cout << "XmatVecTemp(0) " << XmatVecTemp(0) << endl;
		if(tauVec(1)!=0){
                	Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_Surv_LOCO_new2(wVec, tauVec, XmatVecTemp, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);

		}else{
			//Sigma_iX1.col(i) = getProdWminusUb_Surv(XmatVecTemp, RvecIndex, Nvec, sqrtDVec, wVec);
			Sigma_iX1.col(i) = x0Vec;
		}
        //      cout << "Sigma_iX1(0,i) " << Sigma_iX1(0,i) << endl;

        }
        return(Sigma_iX1);
}


// [[Rcpp::export]]
arma::fmat getSigma_X_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;

        //cout << colNumX << endl;
        //cout << size(wVec) << endl;
        //cout << size(tauVec) << endl;


        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
                Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);
        }
        return(Sigma_iX1);
}



// [[Rcpp::export]]
arma::fvec  getSigma_G(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG){
  	arma::fvec Sigma_iG;
  	Sigma_iG = getPCG1ofSigmaAndVector(wVec, tauVec, Gvec, maxiterPCG, tolPCG);
  	return(Sigma_iG);
}


// [[Rcpp::export]]
arma::fvec  getSigma_G_Surv(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec,  arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;

	arma::fvec x0Vec = wVec % Gvec - sqrtDRN.t() * sqrtDRN * Gvec;

	if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv(wVec, tauVec, Gvec, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
        }else{
                //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }

        return(Sigma_iG);
}

// [[Rcpp::export]]
arma::fvec  getSigma_G_Surv_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec,  arma::fmat & WinvNRt, arma::fmat & ACinv, arma::fvec & diagofWminusUinv, arma::fmat & sqrtDRN, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;

	arma::fvec x0Vec = wVec % Gvec - sqrtDRN.t() * sqrtDRN * Gvec;

        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv_LOCO(wVec, tauVec, Gvec, WinvNRt, ACinv, diagofWminusUinv, x0Vec, maxiterPCG, tolPCG);
        }else{
                //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }


        return(Sigma_iG);
}

// [[Rcpp::export]]
arma::fvec  getSigma_G_Surv_new(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;
	arma::fvec sqrtDVec = arma::sqrt(Dvec);
	arma::fvec x0Vec;
	x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
	

	if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
        	Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG);
	}else{
		//Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);	
		Sigma_iG = x0Vec;
	}
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}


// [[Rcpp::export]]
arma::fvec  getSigma_G_Surv_new2(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
        arma::fvec Sigma_iG;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        arma::fvec x0Vec;
        x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);


        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
		Sigma_iG = getPCG1ofSigmaAndVector_Surv_new2(wVec, tauVec, Gvec, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
        }else{
        //      //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}



// [[Rcpp::export]]
arma::fvec  getSigma_G_Surv_new2_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG, arma::fvec & dofWminusU){
        arma::fvec Sigma_iG;
        arma::fvec sqrtDVec = arma::sqrt(Dvec);
        arma::fvec x0Vec;
        x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);


        if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
                Sigma_iG = getPCG1ofSigmaAndVector_Surv_LOCO_new2(wVec, tauVec, Gvec, RvecIndex, Nvec, sqrtDVec, diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG, dofWminusU);
        }else{
        //      //Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
                Sigma_iG = x0Vec;
        }
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}



// [[Rcpp::export]]
arma::fvec  getSigma_G_Surv_new_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, arma::fvec & WinvN, arma::fvec & Dvec, arma::fvec & diagofWminusUinv, arma::fvec & Nvec, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;
	arma::fvec x0Vec;
	arma::fvec sqrtDVec = arma::sqrt(Dvec);
	x0Vec = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
		

	if(tauVec(1) != 0){
        //Sigma_iG = getPCG1ofSigmaAndVector_Surv_new(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec, maxiterPCG, tolPCG, Rmat);
        	Sigma_iG = getPCG1ofSigmaAndVector_Surv_new_LOCO(wVec, tauVec, Gvec, RvecIndex, sqrtWinvNVec, WinvN, Dvec,diagofWminusUinv, x0Vec,  maxiterPCG, tolPCG);
	}else{
		Sigma_iG = x0Vec;
        //        Sigma_iG = getProdWminusUb_Surv(Gvec, RvecIndex, Nvec, sqrtDVec, wVec);
	}
        //cout << "Sigma_iG: " << Sigma_iG << endl;
        return(Sigma_iG);
}




// [[Rcpp::export]]
arma::fvec  getSigma_G_LOCO(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;
        Sigma_iG = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Gvec, maxiterPCG, tolPCG);
        return(Sigma_iG);
}

//This function needs the function getPCG1ofSigmaAndVector and function getCrossprodMatAndKin
// [[Rcpp::export]]
arma::fvec GetTrace_q(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::fmat& cov1,  int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){
  	set_seed(200);
  	arma::fmat Sigma_iXt = Sigma_iX.t();
  	int Nnomissing = geno.getNnomissing();
  	arma::fvec tempVec(nrun);
  	tempVec.zeros();
  	arma::fvec tempVec0(nrun);
  	tempVec0.zeros();

        arma::fvec Sigma_iu;
        arma::fcolvec Pu;
        arma::fvec Au;
        arma::fvec uVec;

        int nrunStart = 0;
        int nrunEnd = nrun;
        float traceCV = traceCVcutoff + 0.1;
        float traceCV0 = traceCVcutoff + 0.1;

        while((traceCV > traceCVcutoff) | (traceCV0 > traceCVcutoff)){


  	for(int i = nrunStart; i < nrunEnd; i++){

    		Rcpp::NumericVector uVec0;
    		uVec0 = nb(Nnomissing);
    		uVec = as<arma::fvec>(uVec0);
    		uVec = uVec*2 - 1;
  //  		arma::fvec Sigma_iu;
    		Sigma_iu = getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG);
  //  		arma::fcolvec Pu;
    		Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
  //  		arma::fvec Au;
    		Au = getCrossprodMatAndKin(uVec);
    		tempVec(i) = dot(Au, Pu);
    		tempVec0(i) = dot(uVec, Pu);
                Au.clear();
      		Pu.clear();
      		Sigma_iu.clear();
      		uVec.clear();
  	}
	traceCV = calCV(tempVec);
	traceCV0 = calCV(tempVec0);
	
	if((traceCV > traceCVcutoff) | (traceCV0 > traceCVcutoff)){
          nrunStart = nrunEnd;
          nrunEnd = nrunEnd + 10;
          tempVec.resize(nrunEnd);
          tempVec0.resize(nrunEnd);

	  //std::cout << "arma::mean(tempVec0): " << arma::mean(tempVec0) << std::endl;	
          cout << "CV for trace random estimator using "<< nrunStart << " runs is " << traceCV <<  "(> " << traceCVcutoff << endl;
          cout << "try " << nrunEnd << "runs" << endl;
        }

    }   

  	arma::fvec traVec(2);
  	traVec(1) = arma::mean(tempVec);
  	traVec(0) = arma::mean(tempVec0);
	tempVec.clear();
	tempVec0.clear();
  	return(traVec);
}

//Rcpp::List getAIScore_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){


//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace
// [[Rcpp::export]]
Rcpp::List getAIScore_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){


  	arma::fmat Sigma_iXt = Sigma_iX.t();
  	arma::fmat Xmatt = Xmat.t();

	//arma::fmat cov1;
        //try {
        //  cov1 = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
        //} catch (const std::exception& e) {
        //  cov1 = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
        //  cout << "inv_sympd failed, inverted with pinv" << endl;
        //}


  	//arma::fmat cov1 = inv_sympd(Xmatt * Sigma_iX);
  	arma::fvec PY1 = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));
  	arma::fvec APY = getCrossprodMatAndKin(PY1);

  	float YPAPY = dot(PY1, APY);

  	arma::fvec A0PY = PY1; ////Quantitative


  	float YPA0PY = dot(PY1, A0PY); ////Quantitative


  	arma::fvec Trace = GetTrace_q(Sigma_iX, Xmat, wVec, tauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);

  	arma::fmat AI(2,2);
  	arma::fvec PA0PY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, A0PY, maxiterPCG, tolPCG);
  	arma::fvec PA0PY = PA0PY_1 - Sigma_iX * (cov * (Sigma_iXt * PA0PY_1));

  	AI(0,0) =  dot(A0PY, PA0PY);

  	//cout << "A1(0,0) " << AI(0,0)  << endl;
  	arma::fvec PAPY_1 = getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG, tolPCG);
  	arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
  	AI(1,1) = dot(APY, PAPY);

  	AI(0,1) = dot(A0PY, PAPY);

  	AI(1,0) = AI(0,1);

  	//cout << "AI " << AI << endl;
  	//cout << "Trace " << Trace << endl;
  	//cout << "YPAPY " << YPAPY << endl;
  	//cout << "cov " << cov1 << endl;
	return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("YPA0PY") = YPA0PY,Named("Trace") = Trace,Named("PY") = PY1,Named("AI") = AI);

}






//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace
// [[Rcpp::export]]
Rcpp::List getAIScore_q_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){

        int Nnomissing = geno.getNnomissing();
        arma::fvec Sigma_iY1;
        Sigma_iY1 = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, Yvec, maxiterPCG, tolPCG);

//	for(int j = 0; j < 10; j++){
//                std::cout << "Sigma_iY1(j): " << Sigma_iY1(j) << std::endl;
//        }


        int colNumX = Xmat.n_cols;
        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);

                Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_LOCO(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);

        }


        arma::fmat Sigma_iX1t = Sigma_iX1.t();
        arma::fmat Xmatt = Xmat.t();

        //arma::fmat cov1 = inv_sympd(Xmatt * Sigma_iX1);
	arma::fmat cov1;
        try {
          cov1 = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX1));
        } catch (const std::exception& e) {
          cov1 = arma::pinv(arma::symmatu(Xmatt * Sigma_iX1));
          cout << "inv_sympd failed, inverted with pinv" << endl;
        }


        //cout << "cov " << cov1 << endl;


	return Rcpp::List::create(Named("cov") = cov1, Named("Sigma_iX") = Sigma_iX1, Named("Sigma_iY") = Sigma_iY1);
        //return Rcpp::List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace,Named("Sigma_iY") = Sigma_iY1, Named("Sigma_iX") = Sigma_iX1, Named("PY") = PY1, Named("AI") = AI, Named("cov") = cov1);
}



//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore_q
// [[Rcpp::export]]
Rcpp::List fitglmmaiRPCG_q_LOCO(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){

        arma::uvec zeroVec = (tauVec < tol); //for Quantitative, GMMAT
        Rcpp::List re = getAIScore_q_LOCO(Yvec, Xmat, wVec, tauVec, nrun, maxiterPCG, tolPCG, traceCVcutoff);
//return Rcpp::List::create(Named("cov") = cov1, Named("Sigma_iX") = Sigma_iX1, Named("Sigma_iY") = Sigma_iY1);
        arma::fmat cov = re["cov"];
        arma::fmat Sigma_iX = re["Sigma_iX"];
        arma::fmat Sigma_iXt = Sigma_iX.t();

        arma::fvec alpha1 = cov * (Sigma_iXt * Yvec);
        arma::fvec Sigma_iY = re["Sigma_iY"];
        arma::fvec eta1 = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha1) / wVec;
	return List::create(Named("tau") = tauVec, Named("cov") = cov, Named("alpha") = alpha1, Named("eta") = eta1);
}



//Rcpp::List fitglmmaiRPCG_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){



//This function needs the function getPCG1ofSigmaAndVector and function getCrossprod, getAIScore_q
// [[Rcpp::export]]
Rcpp::List fitglmmaiRPCG_q(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec &tauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff){

  	arma::uvec zeroVec = (tauVec < tol); //for Quantitative, GMMAT
	Rcpp::List re = getAIScore_q(Yvec, Xmat,wVec,  tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);

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
  	//cout << "0,0" << AI(0,0) << endl;
  	//cout << "0,1" << AI(0,1) << endl;
  	//cout << "1,0" << AI(1,0) << endl;
  	//cout << "1,1" << AI(1,1) << endl;

  	arma::fvec Dtau = solve(AI, scoreVec);


  	arma::fvec tau0 = tauVec;
  	tauVec = tau0 + Dtau;


  	tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros(); //for Quantitative Copied from GMMAT  

  	float step = 1.0;


  	//cout << "tau2 " << tauVec(0) << " " << tauVec(1) << endl;
  	while (tauVec(0) < 0.0 || tauVec(1)  < 0.0){ //for Quantitative
     	//	cout << "tauVec Here: " << tauVec << endl;
    		step = step*0.5;
    		tauVec = tau0 + step * Dtau; //for Quantitative
    	//	cout << "tau_4: " << tauVec << endl;
    		tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros(); //for Quantitative Copied from GMMAT
    	//	cout << "tau_5: " << tauVec << endl;
 	}



  	tauVec.elem( find(tauVec < tol) ).zeros();
	return List::create(Named("tau") = tauVec);

  	//return List::create(Named("tau") = tauVec, Named("cov") = cov, Named("alpha") = alpha1, Named("eta") = eta1);
}



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_usingSubMarker : public Worker
{
        // source vectors
        arma::fcolvec & m_bVec;
        unsigned int m_N;
        unsigned int m_M_Submarker;
        unsigned int m_M;
        arma::ivec subMarkerIndex ;

        // product that I have accumulated
        arma::fvec m_bout;


        // constructors
        CorssProd_usingSubMarker(arma::fcolvec & y)
                : m_bVec(y) {

                //m_Msub = geno.getMsub();
                subMarkerIndex = getSubMarkerIndex();
                m_M_Submarker = subMarkerIndex.n_elem;
                m_N = geno.getNnomissing();
                m_bout.zeros(m_N);
        }
        CorssProd_usingSubMarker(const CorssProd_usingSubMarker& CorssProd_usingSubMarker, Split)
                : m_bVec(CorssProd_usingSubMarker.m_bVec)
        {

                m_N = CorssProd_usingSubMarker.m_N;
                //m_M = CorssProd_usingSubMarker.m_M;
                m_M_Submarker = CorssProd_usingSubMarker.m_M_Submarker;
                subMarkerIndex = CorssProd_usingSubMarker.subMarkerIndex;
                m_bout.zeros(m_N);

        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                float val1;
                int j;
                for(unsigned int i = begin; i < end; i++){
                        j = subMarkerIndex[i];
//			std::cout << "j: " << j << std::endl;	
                        geno.Get_OneSNP_StdGeno(j, &vec);
                        val1 = dot(vec,  m_bVec);
                        m_bout += val1 * (vec);
                }
        }

        // join my value with that of another InnerProduct
        void join(const  CorssProd_usingSubMarker & rhs) {
        m_bout += rhs.m_bout;
        }
};


// [[Rcpp::export]]
arma::fvec parallelCrossProd_usingSubMarker(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        int m_M_Submarker = getSubMarkerNum();

//	std::cout << "m_M_Submarker: " << m_M_Submarker << std::endl;
        CorssProd_usingSubMarker CorssProd_usingSubMarker(bVec);
//	std::cout << "m_M_Submarker: 2 " << m_M_Submarker << std::endl;
  // call paralleReduce to start the work
        parallelReduce(0, m_M_Submarker, CorssProd_usingSubMarker);
//	std::cout << "m_M_Submarker: 3 " << m_M_Submarker << std::endl;
//	std::cout << "CorssProd_usingSubMarker.m_bout " << CorssProd_usingSubMarker.m_bout << std::endl;
  // return the computed product
        //cout << "Msub: " << Msub << endl;
        //for(int i=0; i<100; ++i)
        //{
        //      cout << (CorssProd_usingSubMarker.m_bout/m_M_Submarker)[i] << ' ';
        //}
//        cout << endl;

//	cout << (CorssProd_usingSubMarker.m_bout).n_elem << endl;
        return CorssProd_usingSubMarker.m_bout/m_M_Submarker;
}



// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin_usingSubMarker(arma::fcolvec& bVec){

        arma::fvec crossProdVec = parallelCrossProd_usingSubMarker(bVec) ;

        return(crossProdVec);
}









//std::vector<int> calGRMvalueUsingSubMarker_forOneInv(int sampleIndex, float relatednessCutoff){
//        //sampleIndex starts with 0
//        std::vector<int> relatedIndex;
//        int Ntotal = geno.getNnomissing();
//        arma::fcolvec bindexvec(Ntotal);
//        bindexvec.zeros();
//        bindexvec(sampleIndex) = 1;
//        arma::fvec crossProdVec = getCrossprodMatAndKin_usingSubMarker(bindexvec);
//        for(int i=sampleIndex; i< Ntotal; i++){
//                if(crossProdVec(i) >= relatednessCutoff){
//                        relatedIndex.push_back(i);
//                }
//        }
//        return(relatedIndex);
//}




//The code below is from http://gallery.rcpp.org/articles/parallel-inner-product/
struct InnerProduct : public Worker
{
   // source vectors
   std::vector<float> x;
   std::vector<float> y;

   // product that I have accumulated
   float product;

   // constructors
   InnerProduct(const std::vector<float> x, const std::vector<float> y)
      : x(x), y(y), product(0) {}
   InnerProduct(const InnerProduct& innerProduct, Split)
      : x(innerProduct.x), y(innerProduct.y), product(0) {}

   // process just the elements of the range I've been asked to
   void operator()(std::size_t begin, std::size_t end) {
      product += std::inner_product(x.begin() + begin,
                                    x.begin() + end,
                                    y.begin() + begin,
                                    0.0);
   }

   // join my value with that of another InnerProduct
   void join(const InnerProduct& rhs) {
     product += rhs.product;
   }
};


// [[Rcpp::export]]
float parallelInnerProduct(std::vector<float> &x, std::vector<float> &y) {

   int xsize = x.size();
   // declare the InnerProduct instance that takes a pointer to the vector data
   InnerProduct innerProduct(x, y);

   // call paralleReduce to start the work
   parallelReduce(0, x.size(), innerProduct);

   // return the computed product
   return innerProduct.product/xsize;
}


// [[Rcpp::export]]
float calGRMValueforSamplePair(arma::ivec &sampleidsVec){
        //std::vector<float> stdGenoforSamples = geno.Get_Samples_StdGeno(sampleidsVec);
        geno.Get_Samples_StdGeno(sampleidsVec);
	//std::cout << "here5" << std::endl;
	//for(int i = 0; i < 10; i++){
	//	std::cout << geno.stdGenoforSamples[i] << " ";
	//}
	//std::cout << std::endl;
	//std::cout << geno.stdGenoforSamples.size() << std::endl;
        int Ntotal = geno.getNnomissing();
        float grmValue;
	std::vector<float> stdGenoforSamples2;
	//std::cout << "here5b" << std::endl;
	//std::cout << sampleidsVec.n_elem << std::endl;
	//std::cout << "here5c" << std::endl;
        if(sampleidsVec.n_elem == 2){
                std::vector<float> s1Vec;
                //s1Vec.zeros(Ntotal);

                std::vector<float> s2Vec;
                //arma::fvec s2Vec;
                //s2Vec.zeros(Ntotal);

                for(int i = 0; i < Ntotal; i++){
                        //s1Vec[i] = stdGenoforSamples[i*2+0];
                        s1Vec.push_back(geno.stdGenoforSamples[i*2]);
                        s2Vec.push_back(geno.stdGenoforSamples[i*2+1]);
                }
                grmValue = parallelInnerProduct(s1Vec, s2Vec);
                //grmValue = innerProductFun(s1Vec, s2Vec);
        }else{
	//	std::cout << "here5d" << std::endl;
	//	std::cout << "geno.stdGenoforSamples.size() " << geno.stdGenoforSamples.size() << std::endl;
		stdGenoforSamples2.clear();
		for (int i=0; i< geno.stdGenoforSamples.size(); i++){
			//std::cout << i << " " << geno.stdGenoforSamples[i] << " ";
        		stdGenoforSamples2.push_back(geno.stdGenoforSamples[i]);
		}
	//	std::cout << std::endl;
	//	std::cout << "here6" << std::endl;
                grmValue = parallelInnerProduct(stdGenoforSamples2, geno.stdGenoforSamples);
                //grmValue = innerProductFun(stdGenoforSamples2, geno.stdGenoforSamples);
	//	std::cout << "here7" << std::endl;
        }
        return(grmValue);
}


//Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
//arma::sp_fmat createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){



// [[Rcpp::export]]
Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){

        int nSubMarker = markerIndexVec.n_elem;
        int Ntotal = geno.getNnomissing();
        std::vector<unsigned int>     iIndexVec;
        std::vector<unsigned int>     iIndexVec2;
        std::vector<unsigned int>     jIndexVec;
        std::vector<unsigned int>     jIndexVec2;
        std::vector<unsigned int>     allIndexVec;
        std::vector<float>     kinValueVec;
        std::vector<float>     kinValueVec2;
	std::vector<float> stdGenoMultiMarkers;	
	stdGenoMultiMarkers.resize(Ntotal*nSubMarker);

	//std::cout << "createSparseKin1" << std::endl;
	size_t sizeTemp;
	float kinValue;
	float kinValueTemp;
	//std::cout << "createSparseKin1b" << std::endl;

	Get_MultiMarkersBySample_StdGeno(markerIndexVec, stdGenoMultiMarkers);
	std::cout << "createSparseKin2" << std::endl;
	//arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), Ntotal, nSubMarker);
	arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), nSubMarker, Ntotal);
	//std::cout << "createSparseKin3" << std::endl;
	//std::cout << "stdGenoMultiMarkersMat.n_rows: " << stdGenoMultiMarkersMat.n_rows << std::endl;
	//std::cout << "stdGenoMultiMarkersMat.n_cols: " << stdGenoMultiMarkersMat.n_cols << std::endl;



        for(unsigned int i=0; i< Ntotal; i++){
              for(unsigned int j = i; j < Ntotal; j++){
                        //kinValueTemp = arma::dot(stdGenoMultiMarkersMat.row(i), stdGenoMultiMarkersMat.row(j));
			if(j > i){
                		kinValueTemp = arma::dot(stdGenoMultiMarkersMat.col(i), stdGenoMultiMarkersMat.col(j));
                		kinValueTemp = kinValueTemp/nSubMarker;
                		if(kinValueTemp >= relatednessCutoff){
//                              if(i == 0){
                                //std::cout << "kinValueTemp: " << kinValueTemp << std::endl;
                                //std::cout << "relatednessCutoff: " << relatednessCutoff << std::endl;
                                //std::cout << "i: " << i << std::endl;
//                              std::cout << "j: " << j;
//                              }
                        		iIndexVec.push_back(i);
					jIndexVec.push_back(j);

                		}
			}else{
				iIndexVec.push_back(i);
				jIndexVec.push_back(j);
			}
        	}
	}
	
	arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
        size_t ni = iIndexVec.size();
        kinValueVec.resize(ni);
        std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

        int Mmarker = geno.getM();
        for(size_t i=0; i< Mmarker; i++){
                geno.Get_OneSNP_StdGeno(i, temp);
                for(size_t j=0; j < ni; j++){
                        kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
                }

        }






/*
//	(stdGenoMultiMarkersMat.row(0)).print("stdGenoMultiMarkersMat.row(0):");
	//std::cout << stdGenoMultiMarkersMat << std::endl;
	//std::cout << stdGenoMultiMarkersMat.row(487) << std::endl;	
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(16); // Use 16 threads for all consecutive parallel regions
	int totalCombination = Ntotal*(Ntotal-1)/2 - 1;

	#pragma omp parallel
	{
	std::vector<unsigned int> vec_privatei;	
	std::vector<unsigned int> vec_privatej;	
	#pragma omp for nowait //fill vec_private in parallel
	for(int k = 0; k < totalCombination; k++){
        	int i = k / Ntotal;
        	int j = k % Ntotal;
        	if((j <= i)){
            		i = Ntotal - i - 2;
            		j = Ntotal - j - 1;
        	}

//        for(i=0; i< Ntotal; i++){
//		for(j = i; j < Ntotal; j++){
			//kinValueTemp = arma::dot(stdGenoMultiMarkersMat.row(i), stdGenoMultiMarkersMat.row(j));
		kinValueTemp = arma::dot(stdGenoMultiMarkersMat.col(i), stdGenoMultiMarkersMat.col(j));
		kinValueTemp = kinValueTemp/nSubMarker;
		if(kinValueTemp >= relatednessCutoff){
//				if(i == 0){
				//std::cout << "kinValueTemp: " << kinValueTemp << std::endl;
				//std::cout << "relatednessCutoff: " << relatednessCutoff << std::endl;
				//std::cout << "i: " << i << std::endl;
//				std::cout << "j: " << j;
//				}
			vec_privatei.push_back((unsigned int)i);
			//allIndexVec.push_back(i);
			//iIndexVec.push_back(i);
			//iIndexVec.push_back(i);
			//allIndexVec.push_back(j);
			vec_privatej.push_back((unsigned int)j);
				
								
		}
	}
//	#pragma omp critical
	#pragma omp for schedule(static) ordered
    	for(int i=0; i<omp_get_num_threads(); i++) {
        	#pragma omp ordered
        	iIndexVec.insert(iIndexVec.end(), vec_privatei.begin(), vec_privatei.end());  
        	jIndexVec.insert(jIndexVec.end(), vec_privatej.begin(), vec_privatej.end());  
    	}
//		}
	}
//	int nall = allIndexVec.size();
//	 std::cout << "nall: " << nall << std::endl;
//	int k = 0;
//	while(k < nall){
	//	std::cout << "allIndexVec[k]: " << k << " " << allIndexVec[k] << std::endl;
	//	std::cout << "allIndexVec[k+1]: " << k+1 << " " << allIndexVec[k+1] << std::endl;
//        	iIndexVec.push_back(allIndexVec[k]);
//                jIndexVec.push_back(allIndexVec[k+1]);
//		k = k + 2;
//        }
//	allIndexVec.clear();

	for(int k = 0; k < Ntotal; k++){
		iIndexVec.push_back((unsigned int)k);
		jIndexVec.push_back((unsigned int)k);
	}

        arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
        size_t ni = iIndexVec.size();
        //size_t ni = nall/2 + Ntotal;
        kinValueVec.resize(ni);
        std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

        int Mmarker = geno.getM();
        for(size_t i=0; i< Mmarker; i++){
                geno.Get_OneSNP_StdGeno(i, temp);
                for(size_t j=0; j < ni; j++){
//                for(size_t k=0; k < nall/2; k++){
                        kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
//                        kinValueVec[j] = kinValueVec[j] + (((*temp)[allIndexVec[k*2]])*((*temp)[allIndexVec[k*2+1]]))/Mmarker;
                }
//		for(size_t k=nall/2; k < ni; k++){
			
//			kinValueVec[j] = kinValueVec[j] + (((*temp)[allIndexVec[k*2]])*((*temp)[allIndexVec[k*2+1]]))/Mmarker;

//		}
        }	


*/   // end of the openMP version 

	std::cout << "ni: " << ni << std::endl;
/*	for(size_t j=0; j < 10; j++){
		std::cout << "iIndexVec[j]: " << iIndexVec[j] << std::endl;
		std::cout << "jIndexVec[j]: " << jIndexVec[j] << std::endl;
		std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
	}
*/
	for(size_t j=0; j < ni; j++){
		if(kinValueVec[j] >= relatednessCutoff){
	//	std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
			kinValueVec[j] = tauVec(1)*kinValueVec[j];
			iIndexVec2.push_back(iIndexVec[j]+1);
			jIndexVec2.push_back(jIndexVec[j]+1);
			if(iIndexVec[j] == jIndexVec[j]){
				kinValueVec[j] = kinValueVec[j] + tauVec(0)/(wVec(iIndexVec[j]));	
			}
			kinValueVec2.push_back(kinValueVec[j]);
		}

	}

//	std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;

	//arma::fvec x(kinValueVec2);
	//arma::umat locations(iIndexVec2);
	//arma::uvec jIndexVec2_b(jIndexVec2);
	//locations.insert_cols(locations.n_cols, jIndexVec2_b); 
	//arma::umat locationst = locations.t();
	//locations.clear();
	
	//create a sparse Sigma
//	arma::sp_fmat sparseSigma(locationst, x);
//	arma::sp_fmat sparseSigmab  = arma::symmatu(sparseSigma);
	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);
//	return sparseSigmab;
}



// [[Rcpp::export]]
arma::fmat getColfromStdGenoMultiMarkersMat(arma::uvec & a){
	return((geno.stdGenoMultiMarkersMat).cols(a));
}

// [[Rcpp::export]]
int getNColStdGenoMultiMarkersMat(){
	return((geno.stdGenoMultiMarkersMat).n_cols);
}

// [[Rcpp::export]]
int getNRowStdGenoMultiMarkersMat(){
        return((geno.stdGenoMultiMarkersMat).n_rows);
}


// [[Rcpp::export]]
void setSubMarkerIndex(arma::ivec &subMarkerIndexRandom){
	geno.subMarkerIndex = subMarkerIndexRandom;
//	std::cout << "(geno.subMarkerIndex).n_elem: " << (geno.subMarkerIndex).n_elem << std::endl;
	int Nnomissing = geno.getNnomissing();
	(geno.stdGenoMultiMarkersMat).set_size(subMarkerIndexRandom.n_elem, Nnomissing);
}

// [[Rcpp::export]]
void setRelatednessCutoff(float a){
	geno.relatednessCutoff = a;
}


// [[Rcpp::export]]
double innerProduct(NumericVector x, NumericVector y) {
   return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}


//Rcpp::List refineKin(std::vector<unsigned int> &iIndexVec, std::vector<unsigned int> & jIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
//Rcpp::List refineKin(arma::imat &iMat, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){

// [[Rcpp::export]]
Rcpp::List refineKin(float relatednessCutoff){
        std::vector<unsigned int>     iIndexVec2;
        std::vector<unsigned int>     jIndexVec2;
//	std::vector<float>     kinValueVec;
        std::vector<float>     kinValueVec2;
 //       std::vector<float>     kinValueVec_orig; //for test original kinship

	arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
	(*temp).clear();
        //size_t ni = iIndexVec.size();
        //size_t ni = iMat.n_rows;
        size_t ni = geno.indiceVec.size();
	std::cout << "ni: " << ni << std::endl;
 
	initKinValueVecFinal(ni);

//	std::cout << "OKK: "  << std::endl;
//        kinValueVec.resize(ni);
//        std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

        int Mmarker = geno.getM();
        
        //for(size_t i=0; i< Mmarker; i++){
        //        geno.Get_OneSNP_StdGeno(i, temp);
        //        for(size_t j=0; j < ni; j++){
        //                kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
        //        }
        //}
	//arma::fvec kinValueVecTemp;
	arma::fvec kinValueVecTemp2;
	arma::fvec GRMvec;
	GRMvec.set_size(ni);
	int Mmarker_mafgr1perc = 0;
  	for(size_t i=0; i< Mmarker; i++){
//		std::cout << "OKKK: "  << std::endl;
//		std::cout << "Mmarker: " << std::endl;

//                geno.Get_OneSNP_StdGeno(i, temp);
		float freqv = geno.alleleFreqVec[i];
		if(freqv >= minMAFtoConstructGRM && freqv <= 1-minMAFtoConstructGRM){
		Mmarker_mafgr1perc = Mmarker_mafgr1perc + 1;

                geno.Get_OneSNP_Geno(i);
		float invstdv = geno.invstdvVec[i];
		geno.setSparseKinLookUpArr(freqv, invstdv);			

		//std::cout << "freqv: " << freqv << std::endl;
		//std::cout << "invstdv: " << invstdv << std::endl;
		//for (int j = 0; j < 3; j++){
		//	std::cout << geno.sKinLookUpArr[j][0] << std::endl;	
		//	std::cout << geno.sKinLookUpArr[j][1] << std::endl;	
		//	std::cout << geno.sKinLookUpArr[j][2] << std::endl;	

		//}
		//std::cout << "geno.m_OneSNP_StdGeno(i) " << geno.m_OneSNP_StdGeno(i) <<  std::endl;	
		//kinValueVecTemp = parallelcalsparseGRM(iMat);
//		parallelcalsparseGRM(iMat, GRMvec);

		parallelcalsparseGRM(GRMvec);
		//std::cout << "kinValueVecTemp.n_elem: " << kinValueVecTemp.n_elem << std::endl;
//		std::cout << "OKKK2: "  << std::endl;
		parallelsumTwoVec(GRMvec);
//		for(size_t j=0; j< ni; j++){
//			(geno.kinValueVecFinal)[j] = (geno.kinValueVecFinal)[j] + GRMvec(j);
//		}
		(*temp).clear();
	   }//if(freqv >= 0.01 && freqv <= 0.99){
		//kinValueVecTemp.clear();
        }



       // for(size_t j=0; j < 100; j++){
       //         std::cout << "iIndexVec[j]: " << iIndexVec[j] << std::endl;
       //         std::cout << "jIndexVec[j]: " << jIndexVec[j] << std::endl;
       //         std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
       // }

	int a1;
	int a2;
        for(size_t j=0; j < ni; j++){
		geno.kinValueVecFinal[j] = (geno.kinValueVecFinal[j]) /(Mmarker_mafgr1perc);

//		std::cout << "j: " << j << " geno.kinValueVecFinal[j]: " << geno.kinValueVecFinal[j] << std::endl;
            //    if(geno.kinValueVecFinal[j] >= relatednessCutoff){
                if((geno.kinValueVecFinal[j]) >= relatednessCutoff){
        //      std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
			//kinValueVec_orig.push_back((geno.kinValueVecFinal)[j]); //for test	
                        //(geno.kinValueVecFinal)[j] = tauVec(1)*(geno.kinValueVecFinal)[j];
                        //(geno.kinValueVecFinal)[j] = tauVec(1)*(geno.kinValueVecFinal)[j];
 				 a1 = (geno.indiceVec)[j].first + 1;
				 a2 = (geno.indiceVec)[j].second + 1;
				 iIndexVec2.push_back(a1);
				 jIndexVec2.push_back(a2);

                        kinValueVec2.push_back((geno.kinValueVecFinal)[j]);
                }

        }



	std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;
	//return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2,  Named("kinValue_orig") = kinValueVec_orig);	
	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);	
}


// [[Rcpp::export]]
Rcpp::List shortenList(arma::imat &iMat, arma::fvec &kinValueVecTemp, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
	        std::vector<unsigned int>     iIndexVec2;
        std::vector<unsigned int>     jIndexVec2;
	std::vector<float>     kinValueVec2;
	size_t ni = iMat.n_rows;

	for(size_t j=0; j < ni; j++){
                if(kinValueVecTemp(j) >= relatednessCutoff){
        //      std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
                        kinValueVecTemp(j) = tauVec(1)*(kinValueVecTemp(j));
                        iIndexVec2.push_back(iMat(j,1)+1);
                        //iIndexVec2.push_back(iIndexVec[j]+1);
                        jIndexVec2.push_back(iMat(j,2)+1);
                        //jIndexVec2.push_back(jIndexVec[j]+1);
        //                if(iIndexVec[j] == jIndexVec[j]){
        //                        kinValueVec[j] = kinValueVec[j] + tauVec(0)/(wVec(iIndexVec[j]));
        //                }

                        if(iMat(j,1) == iMat(j,2)){
                                kinValueVecTemp(j) = kinValueVecTemp(j) + tauVec(0)/(wVec(iMat(j,1)));
                        }

                        kinValueVec2.push_back(kinValueVecTemp(j));
                }

        }

        std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;
	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);

}

// [[Rcpp::export]]
arma::fvec testTime(int i, arma::fcolvec & m_bVec){
	arma::fvec vec;
	arma::fvec mvec;
	std::cout << "i is " << i << std::endl;
	clock_t t_0;
	t_0 = clock();
        geno.Get_OneSNP_StdGeno(i, &vec);
	clock_t t_1;
	t_1 = clock();
	std::cout << "t_1-t_0 is " << t_1-t_0 << std::endl;
        float val1 = dot(vec,  m_bVec);
	clock_t t_2;
	t_2 = clock();
	std::cout << "t_2-t_1 is " << t_2-t_1 << std::endl;
        mvec = val1 * (vec);
	clock_t t_3;
	t_3 = clock();
	std::cout << "t_3-t_2 is " << t_3-t_2 << std::endl;
	return(mvec);
}


// [[Rcpp::export]]
arma::sp_mat gen_sp_v2(const arma::sp_mat& a) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;

    return result;
}


// [[Rcpp::export]]
arma::vec gen_spsolve_v2(const arma::sp_mat& a) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a);
    int r = result.n_rows;
    arma::vec y = arma::linspace<arma::vec>(0, 5, r);	
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    arma::vec x = arma::spsolve( result, y ); 
    	
    return x;
}

// [[Rcpp::export]]
arma::vec gen_spsolve_inR(const arma::sp_mat& a, arma::vec & y) {
    // sparse x sparse -> sparse
    //arma::sp_mat result1 = result * A;
    arma::vec x = arma::spsolve( a, y );

    return x;
}

// [[Rcpp::export]]
arma::fvec get_DiagofKin(){
    int M = geno.getM();
    int Nnomissing = geno.getNnomissing();
    int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();

    arma::fvec x(Nnomissing);

    if(!(geno.setKinDiagtoOne)){
           x  = (*geno.Get_Diagof_StdGeno()) /MminMAF; 
    }else{
	   x  = arma::ones<arma::fvec>(Nnomissing);	
    }	
    return(x);
}





//The code below is modified from http://gallery.rcpp.org/articles/parallel-inner-product/
struct stdgenoVectorScalorProduct : public Worker
{
   // source vectors
   arma::fvec & m_bout;
   float  y;
   //unsigned int m_N;
   int jthMarker;

   // constructors
   stdgenoVectorScalorProduct(const int jth, const float y, arma::fvec & prodVec)
      : jthMarker(jth), y(y), m_bout(prodVec) {
        //m_N = geno.getNnomissing();
//      m_bout.zeros(m_N);

  }


   // process just the elements of the range I've been asked to

        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                geno.Get_OneSNP_StdGeno(jthMarker, &vec);
                for(unsigned int i = begin; i < end; i++){
                        m_bout[i] = m_bout[i]+vec[i] * y;
                }
        }



};


// [[Rcpp::export]]
void getstdgenoVectorScalorProduct(int jth, float y, arma::fvec & prodVec) {


   stdgenoVectorScalorProduct stdgenoVectorScalorProduct(jth, y, prodVec);

   unsigned int m_N = geno.getNnomissing();

   parallelFor(0, m_N, stdgenoVectorScalorProduct);

   // return the computed product
}





struct getP_mailman : public Worker
{
        // source vectors
        unsigned int ithMarker;
	unsigned int powVal;
	arma::ivec ithGeno;
        // destination vector
        arma::ivec Psubvec;


        // constructors
        getP_mailman(unsigned int ith, unsigned int mmchunksize)
                : ithMarker(ith){
		ithGeno = Get_OneSNP_Geno(ith);			
		//unsigned int k  = pow(3, mmchunksize);
		unsigned m_N = geno.getNnomissing();
		Psubvec.zeros(m_N);
		unsigned int powNumber = mmchunksize - 1 - ith % mmchunksize; 		
		powVal = pow(3, powNumber);
        }


	// take the square root of the range of elements requested
     void operator()(std::size_t begin, std::size_t end) {

	for(unsigned int j = begin; j < end; j++){
		Psubvec[j] = ithGeno[j] * powVal;
        }		 
     }

};


int computePindex(arma::ivec &ithGeno){
	int a = ithGeno.n_elem;
	int q = 0;
	int baseNum;
	for(unsigned int i = 0; i < a; i++){
		baseNum = pow(3, a - i - 1);
		q = q + ithGeno[i] * baseNum;
	}
	return(q);
}


struct getP_mailman_NbyM : public Worker
{
        // source vectors
        unsigned int jthChunk;
        unsigned int mmchunksize;
        //arma::ivec ithGeno;
        // destination vector
        arma::ivec Psubvec;


        // constructors
        getP_mailman_NbyM(unsigned int jthChunk,unsigned int mmchunksize)
                : jthChunk(jthChunk), mmchunksize(mmchunksize){
                //ithGeno = Get_OneSNP_Geno(ith);
                //unsigned int k  = pow(3, mmchunksize);
                unsigned m_M = geno.getM();
                Psubvec.zeros(m_M);
                //powNumber = mmchunksize - 1 - ith % mmchunksize;
        }


        // take the square root of the range of elements requested
     void operator()(std::size_t begin, std::size_t end) {
	arma::ivec ithGeno;	
	arma::ivec ithGenosub;
	unsigned int jthIndvStart = jthChunk * mmchunksize;
	unsigned int jthIndvEnd = (jthChunk+1) * mmchunksize - 1;
	arma::uvec indvIndex = arma::linspace<arma::uvec>(jthIndvStart, jthIndvEnd);
        for(unsigned int i = begin; i < end; i++){
		ithGeno = Get_OneSNP_Geno(i);		
		ithGenosub = ithGeno.elem(indvIndex);
		Psubvec[i] = computePindex(ithGenosub);
        }
     }

};



// // [[Rcpp::export]]
//arma::ivec parallelmmGetP(unsigned int ith, unsigned int mmchunksize) {
  
//  	int M = geno.getM();
//	int N = geno.getNnomissing();	
//  	Pvec.zeros(N);

//  	getP_mailman getP_mailman(ith, mmchunksize);
  
//  	parallelFor(0, N, getP_mailman);
 	
//  	return getP_mailman.Psubvec;
//}


// [[Rcpp::export]]
void sumPz(arma::fvec & Pbvec, arma::fvec & Ubvec, unsigned int mmchunksize){

        for (int i = 0; i < Pbvec.n_elem; i++){
                std::cout << "i: " << i << " " << Pbvec[i] << std::endl;
        }

        unsigned int d = Pbvec.n_elem;;
        Ubvec.zeros(mmchunksize);
        unsigned int i = 0;
        arma::fvec z0;
        arma::fvec z1;
        arma::fvec z2;
        z0.zeros(d/3);
        z1.zeros(d/3);
        z2.zeros(d/3);

        while(i < mmchunksize){
                d = d / 3;
//              std::cout << "d: " << d << std::endl;
                z0.resize(d);
                z1.resize(d);
                z2.resize(d);

//              arma::uvec indexvec = arma::linspace<arma::uvec>(0, d-1);
                z0 = Pbvec.subvec(0, d-1);
/*
                 for (int j = 0; j < z0.n_elem; j++){
                std::cout << "j: " << j << " " << z0[j] << std::endl;
        }
*/
                //indexvec = arma::linspace<arma::uvec>(d, 2*d-1);
                //z1 = Pbvec.elem(indexvec);
                z1 = Pbvec.subvec(d, 2*d-1);
                //indexvec = arma::linspace<arma::uvec>(2*d, 3*d-1);
                //z2 = Pbvec.elem(indexvec);
                z2 = Pbvec.subvec(2*d, 3*d-1);

                Pbvec.resize(d);
                Pbvec = z0 + z1 + z2;
                Ubvec[i] = sum(z1) + 2*sum(z2);
                i = i + 1;
              std::cout << "i: " << i << std::endl;
              std::cout << "Ubvec[i]: " << Ubvec[i] << std::endl;

        }
}



// [[Rcpp::export]]
void mmGetPb_MbyN(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec, arma::fvec & kinbvec) {
	std::cout << "OKKK" << std::endl;
        int M = geno.getM();
        int N = geno.getNnomissing();
	int k = pow(3,mmchunksize);
        arma::ivec Pvec;
	Pvec.zeros(N);
	Pbvec.zeros(k);
	arma::ivec ithGeno;
	ithGeno.ones(N);
	unsigned int Ptemp;
	Ptemp = 1;
	int indL = cthchunk*mmchunksize;
	int indH = (cthchunk+1)*mmchunksize - 1;
	unsigned int j0 = 0;
	//arma::fmat stdGenoMat(mmchunksize, N);
	float ithfreq; 
	float ithinvstd; 
	arma::fvec chunkfreq = geno.alleleFreqVec.subvec(indL, indH);
	arma::fvec chunkinvstd = geno.invstdvVec.subvec(indL, indH); 
	arma::fvec chunkbvec = bvec.subvec(indL, indH); 

	for (int i = indH; i >= indL; i--){
		ithGeno = Get_OneSNP_Geno(i);
		cout << "Ptemp: " << Ptemp << endl;
		//ithfreq = geno.alleleFreqVec(i);
		//ithinvstd = geno.invstdvVec(i);
		Pvec = Pvec + Ptemp * ithGeno; 
		Ptemp = Ptemp * 3;
		//stdGenoMat.row(j) = ithGeno*ithinvstd - 2*ithfreq*ithinvstd;
		//j0 = j0 + 1;

                //unsigned int k  = pow(3, mmchunksize);
                //unsigned m_N = geno.getNnomissing();
                //Psubvec.zeros(m_N);
                //unsigned int powNumber = mmchunksize - 1 - ith % mmchunksize;

	
	//	getP_mailman getP_mailman(i, mmchunksize);
	//	parallelFor(0, N, getP_mailman);
	//	Pvec = Pvec + getP_mailman.Psubvec;
	//	getP_mailman.Psubvec.clear();
  	}
	

	for (int i = 0; i < N; i++){	
//		std::cout << "i: " << i << " " << Pvec[i] << std::endl;	
		Pbvec[Pvec[i]] = Pbvec[Pvec[i]] + bvec[i];
//		std::cout << "Pbvec[Pvec[i]] " << Pbvec[Pvec[i]] << std::endl;
	}

	arma::fvec Gbvectemp;
	sumPz(Pbvec, Gbvectemp, mmchunksize);
	arma::fvec crossKinVec;
	arma::fvec GbvecInvStd = Gbvectemp % chunkinvstd;
        arma::fvec secondTerm = 2*chunkfreq % chunkinvstd * sum(chunkbvec);
        crossKinVec  = GbvecInvStd - secondTerm;

	//getstdgenoVectorScalorProduct(j, crossKinVec[j], kinbvec);
	j0 = 0;
	arma::fvec stdvec;
	for (int i = indL; i <= indH; i++){
		geno.Get_OneSNP_StdGeno(i, &stdvec);
                kinbvec = kinbvec + crossKinVec[j0]*(stdvec);
		j0 = j0 + 1;
	}

//	for (int i = 0; i < k; i++){
//                std::cout << "Pbvec[i]: " << i << " " << Pbvec[i] << std::endl;
//        }

        //return Pbvec;
}

// [[Rcpp::export]]
void mmGetPb_NbyM(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec) {

        int M = geno.getM();
        int N = geno.getNnomissing();
        int k = pow(3,mmchunksize);
        arma::ivec Pvec;
        Pvec.zeros(M);
        Pbvec.zeros(k);
	getP_mailman_NbyM getP_mailman_NbyM(cthchunk,mmchunksize);
	parallelFor(0, M, getP_mailman_NbyM);
	Pvec = getP_mailman_NbyM.Psubvec;
	for (int i = 0; i < M; i++){
		Pbvec[Pvec[i]] = Pbvec[Pvec[i]] + bvec[i];
	}
}



// [[Rcpp::export]]
void muliplyMailman(arma::fvec & bvec, arma::fvec & Gbvec, arma::fvec & kinbvec){
	int M = geno.getM();
        int N = geno.getNnomissing();

        Gbvec.zeros(M);
	std::cout << "Gbvec.n_elem " << Gbvec.n_elem << std::endl;
        unsigned int mmchunksize = ceil(log(N)/log(3));
	std::cout << "mmchunksize " << mmchunksize << std::endl;

        int numchunk = M / mmchunksize; 
	std::cout << "numchunk " << numchunk << std::endl;
        int reschunk = M % mmchunksize;
	std::cout << "reschunk " << reschunk << std::endl;
	//unsigned int indL;
	//unsigned int indH;
	//mmGetPb(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec)
	arma::fvec Pbvec;
	//arma::fvec Gbvectemp;	


	
	//for (unsigned int j = 0; j < 1; j++){
	for (unsigned int j = 0; j < numchunk; j++){
//		std::cout << "j: " << j << std::endl;
		//Pbvec.zeros(M);
		//indL = j*mmchunksize;
		//indH = (j+1)*mmchunksize-1;
//		if(j == 0){
		double wall0ain = get_wall_time();
 		double cpu0ain  = get_cpu_time();
//		}
//		mmGetPb_MbyN(j, mmchunksize, bvec, Pbvec);

//		if(j == 0){

		mmGetPb_MbyN(j, mmchunksize, bvec, Pbvec, kinbvec);



	double wall1ain = get_wall_time();
 double cpu1ain  = get_cpu_time();
 cout << "Wall Time in mmGetPb_MbyN = " << wall1ain - wall0ain << endl;
 cout << "CPU Time  in mmGetPb_MbyN = " << cpu1ain - cpu0ain  << endl;


//}

//		sumPz(Pbvec, Gbvectemp, mmchunksize);

//if(j == 0){
cout << "ith chunk " << j << endl;
//}
		//getstdgenoVectorScalorProduct(int jth, float y, arma::fvec & prodVec)

//		Gbvec.subvec(j*mmchunksize, (j+1)*mmchunksize-1) = Gbvectemp;
  	}

        if(reschunk > 0){
			arma::fvec vec;
		//arma::uvec indexvec = arma::linspace<arma::uvec>(M-reschunk, M-1);
		for (unsigned int j = M-reschunk; j < M; j++){
     		           geno.Get_OneSNP_StdGeno(j, &vec);
			kinbvec = kinbvec + arma::dot(vec, bvec) * vec;
		}	
        }

	kinbvec = kinbvec / M;
}


// [[Rcpp::export]]
void muliplyMailman_NbyM(arma::fvec & bvec, arma::fvec & tGbvec){
        int M = geno.getM();
        int N = geno.getNnomissing();

        tGbvec.zeros(N);

        unsigned int mmchunksize = ceil(log(M)/log(3));

        int numchunk = N / mmchunksize;
        int reschunk = N % mmchunksize;
        unsigned int indL;
        unsigned int indH;
        //mmGetPb(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec)
        arma::fvec Pbvec;
        Pbvec.zeros(M);
	arma::fvec tGbvectemp;

        for (unsigned int j = 0; j < numchunk; j++){
                indL = j*mmchunksize;
                indH = (j+1)*mmchunksize-1;
		mmGetPb_NbyM(j, mmchunksize, bvec, Pbvec);
           	sumPz(Pbvec, tGbvectemp, mmchunksize);
                tGbvec.subvec(j*mmchunksize, (j+1)*mmchunksize-1) = tGbvectemp;     
        }

        if(reschunk > 0){
		arma::imat A(reschunk,M);
		A.zeros();
		arma::ivec Gtemp(N);
		Gtemp.zeros();
		arma::ivec Gtemp2(reschunk);
		Gtemp2.zeros();
		arma::uvec indexvec = arma::linspace<arma::uvec>(M - reschunk -1, M);
                for (unsigned int j = 0; j < M; j++){
			Gtemp = Get_OneSNP_Geno(j);
			Gtemp2 = Gtemp.elem(indexvec);
			A.col(j) = Gtemp2;
                }

		Pbvec.elem(indexvec) = Gtemp2 * (bvec.elem(indexvec));
        }
}

// [[Rcpp::export]]
void freqOverStd(arma::fcolvec& freqOverStdVec){
	freqOverStdVec = 2 * (geno.alleleFreqVec) % (geno.invstdvVec);

	 int M = geno.getM();
/*
	for (unsigned int j = 0; j < M; j++){
		std::cout << "geno.alleleFreqVec " << j << " " << geno.alleleFreqVec[j] << std::endl; 
		std::cout << "geno.invstdvVec " << j << " " << geno.invstdvVec[j] << std::endl; 
		std::cout << "freqOverStdVec " << j << " " << freqOverStdVec[j] << std::endl; 
               }
*/

}


// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin_mailman(arma::fcolvec& bVec){
	std::cout << "b0: " << std::endl;
	int M = geno.getM();
        int N = geno.getNnomissing();
	arma::fvec Gbvec;


	double wall0in = get_wall_time();
 	double cpu0in  = get_cpu_time();
 	arma::fvec kinbvec;
        kinbvec.zeros(N);

	muliplyMailman(bVec, Gbvec, kinbvec);


double wall1in = get_wall_time();
 double cpu1in  = get_cpu_time();
 cout << "Wall Time in muliplyMailman = " << wall1in - wall0in << endl;
 cout << "CPU Time  in muliplyMailman = " << cpu1in - cpu0in  << endl;



//	for (unsigned int j = 0; j < M; j++){
//                std::cout << "Gbvec " << j << " " << Gbvec[j] << std::endl;
//               }
/*
//	std::cout << "b: " << std::endl;
	arma::fvec freqOverStdVec;
//	std::cout << "a: " << std::endl;
	freqOverStd(freqOverStdVec);
//	std::cout << "c: " << std::endl;
	arma::fvec crossKinVec;
	arma::fvec GbvecInvStd = Gbvec % (geno.invstdvVec);
	arma::fvec secondTerm = freqOverStdVec * sum(bVec);
	crossKinVec  = GbvecInvStd - secondTerm;

double wall2in = get_wall_time();
 double cpu2in  = get_cpu_time();
 cout << "Wall Time in Gtb = " << wall2in - wall1in << endl;
 cout << "CPU Time  in Gtb = " << cpu2in - cpu1in  << endl;


	 for (unsigned int j = 0; j < M; j++){
                std::cout << "GbvecInvStd " << j << " " << GbvecInvStd[j] << std::endl;
                std::cout << "secondTerm " << j << " " << secondTerm[j] << std::endl;
		std::cout << "crossKinVec " << j << " " << crossKinVec[j] << std::endl;
               }
*/
/*	
	arma::fvec kinbvec;
	kinbvec.zeros(N);

	for (unsigned int j = 0; j < M; j++){
		getstdgenoVectorScalorProduct(j, crossKinVec[j], kinbvec);
	}


double wall3in = get_wall_time();
 double cpu3in  = get_cpu_time();
 cout << "Wall Time in getstdgenoVectorScalorProduct = " << wall3in - wall2in << endl;
 cout << "CPU Time  in getstdgenoVectorScalorProduct = " << cpu3in - cpu2in  << endl;



	kinbvec = kinbvec / M;
*/	
        return(kinbvec);

}


// [[Rcpp::export]]
arma::fvec get_GRMdiagVec(){
  int mMarker = gettotalMarker();
  int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        //cout << "MminMAF=" << MminMAF << endl;

  arma::fvec diagGRMVec = (*geno.Get_Diagof_StdGeno())/MminMAF;
  return(diagGRMVec);
}


// [[Rcpp::export]]
void setminMAFforGRM(float minMAFforGRM){
  minMAFtoConstructGRM = minMAFforGRM;
}

// // [[Rcpp::export]] 
//int getNumofMarkersforGRM(){
//  int a = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
//  return(a);
//}
