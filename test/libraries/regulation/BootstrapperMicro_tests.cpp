#include <gtest/gtest.h>

#include <genetrail2/core/Matrix.h>
#include <config.h>

#include <genetrail2/core/NameDatabases.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <Eigen/Core>

#include <vector>
#include <genetrail2/regulation/RegulationFile.h>
#include <genetrail2/regulation/RegulationFileParser.h>

#include <genetrail2/regulation/RegulationBootstrapperMicro.h>
#include <genetrail2/regulation/RegulatorAssociationScore.h>


using namespace GeneTrail;

TEST(BootstrapperMicro, perform_bootstrapping_run) {
  unsigned int opts = DenseMatrixReader::NO_OPTIONS;
	opts |= DenseMatrixReader::READ_ROW_NAMES;
    opts |= DenseMatrixReader::READ_COL_NAMES;
    DenseMatrixReader reader;
    std::ifstream strm(TEST_DATA_PATH("RegulationFile_Matrix.txt"), std::ios::binary);
    std::ifstream strm_micro(TEST_DATA_PATH("RegulationFile_MatrixMicro.txt"), std::ios::binary);

    DenseMatrix matrix = reader.read(strm, opts);
    DenseMatrix micro = reader.read(strm_micro,opts);
  
    MatrixNameDatabase name_database_matrix(&matrix);
    MatrixNameDatabase name_database_micro(&micro);  
    
    std::unordered_set<size_t> tset;
      
      
    RegulationFileParser<MatrixNameDatabase, double> parser(name_database_matrix, name_database_micro, tset, TEST_DATA_PATH("RegulationFileMicro.txt"), 0.0);
    RegulationFile<double>& rFile = parser.getRegulationFile();
    std::vector<std::tuple<size_t,size_t,double>> regulations;
    std::tuple<size_t,size_t,double> t1(0,0,2.5);
    regulations.push_back(t1);
    std::tuple<size_t,size_t,double> t2(0,1,2.0);
    regulations.push_back(t2);
    std::tuple<size_t,size_t,double> t3(1,0,1.7);
    regulations.push_back(t3);
    
    for(auto& reg :regulations){
	EXPECT_TRUE( std::get<2>(reg) == 2.5 || std::get<2>(reg) == 2.0 || std::get<2>(reg) == 1.7);
    }
      
      
    RegulationBootstrapperMicro<double> bootstrapper(&matrix,&micro,5,0);
    bootstrapper.perform_bootstrapping_run(rFile,regulations,false,false,false,PearsonCorrelation());

    for(auto& reg :regulations){      
      if(std::get<0>(reg) == name_database_micro("miRNA2") && std::get<1>(reg) == name_database_matrix("GeneA")){
	EXPECT_TRUE( std::get<2>(reg) == (0.55)/sqrt(0.91));
      }
    }
    std::vector<size_t> vec{0,0,1};
    bootstrapper.set_bootstrap_sample(vec);
    bootstrapper.perform_bootstrapping_run(rFile,regulations,false,false,false,PearsonCorrelation());
    for(auto& reg :regulations){      
      if(std::get<0>(reg) == name_database_micro("miRNA2") && std::get<1>(reg) == name_database_matrix("GeneA")){
	EXPECT_EQ(1.0000000000000001, std::get<2>(reg));
      }
    }

}