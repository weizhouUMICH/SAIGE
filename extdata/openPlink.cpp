// reading a text file
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main () {
  string line;
  std::ifstream m_ibedFile;
  std::string m_bedFile = "/net/dumbo/home/zhowei/projects/Dec2021/SAIGE/extdata/input/genotype_100markers.bed"
  m_ibedFile.open(m_bedFile.c_str(), std::ios::binary);

  if (m_ibedFile.is_open())
  {
    	  
    m_ibedFile.close();
  }

  else cout << "Unable to open file"; 

  return 0;
}
