#include <gtest/gtest.h>

#include <genetrail2/core/NameDatabases.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>

#include <genetrail2/regulation/RegulationFile.h>
#include <genetrail2/regulation/RegulationFileParser.h>

#include <config.h>
#include <unordered_set>
#include <tuple>

using namespace GeneTrail;

TEST(RegulationFile, MapNameDatabase) {
    MapNameDatabase name_database(TEST_DATA_PATH("RegulationFile.txt"));
    EXPECT_EQ(name_database.size(), 5);
    
    EXPECT_EQ(name_database(0), "GeneA");
    EXPECT_EQ(name_database("GeneA"), 0);
    
    EXPECT_EQ(name_database(1), "GeneB");
    EXPECT_EQ(name_database("GeneB"), 1);
    
    EXPECT_EQ(name_database(2), "GeneC");
    EXPECT_EQ(name_database("GeneC"), 2);
    
    EXPECT_EQ(name_database(3), "GeneD");
    EXPECT_EQ(name_database("GeneD"), 3);
    
    EXPECT_EQ(name_database(4), "GeneE");
    EXPECT_EQ(name_database("GeneE"), 4);
    
    std::unordered_set<size_t> tset;
    for(size_t i=0; i<5; ++i){
        tset.emplace(i);
    }
    
    RegulationFileParser<MapNameDatabase, double> parser(name_database, tset, TEST_DATA_PATH("RegulationFile.txt"), 0.0);
	RegulationFile<double>& regulationFile = parser.getRegulationFile();
    
    for(const auto& i : regulationFile.regulators()){
        EXPECT_TRUE(name_database(i) == "GeneA" || name_database(i) == "GeneC" || name_database(i) == "GeneD");
    }
    
    EXPECT_TRUE(regulationFile.checkRegulator(name_database("GeneA")) || regulationFile.checkRegulator(name_database("GeneC")) || regulationFile.checkRegulator(name_database("GeneD")));
    EXPECT_FALSE(regulationFile.checkRegulator(name_database("GeneB")));
    EXPECT_FALSE(regulationFile.checkRegulator(name_database("GeneE")));
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database("GeneA")).size(), 3);
    for(const auto& reg : regulationFile.regulator2regulations(name_database("GeneA"))){
        EXPECT_TRUE(std::get<0>(reg) == name_database("GeneA"));
        EXPECT_TRUE(name_database(std::get<1>(reg)) == "GeneB" || name_database(std::get<1>(reg)) == "GeneC" || name_database(std::get<1>(reg)) == "GeneE");
        if(name_database(std::get<1>(reg)) == "GeneB"){
            EXPECT_EQ(std::get<2>(reg), 1.0);
        } else if(name_database(std::get<1>(reg)) == "GeneC"){
            EXPECT_EQ(std::get<2>(reg), 2.0);
        } else if(name_database(std::get<1>(reg)) == "GeneE"){
            EXPECT_EQ(std::get<2>(reg), 5.0);
        } else {
            // If this happens we are screwed!!
            EXPECT_TRUE(false);
        }
    }
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database("GeneC")).size(), 1);
    auto reg = regulationFile.regulator2regulations(name_database("GeneC"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneD");
    EXPECT_EQ(std::get<2>(reg), 3.0);
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database("GeneD")).size(), 1);
    reg = regulationFile.regulator2regulations(name_database("GeneD"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneB");
    EXPECT_EQ(std::get<2>(reg), 4.0);
    
    EXPECT_TRUE(regulationFile.checkTarget(name_database("GeneB")) || regulationFile.checkTarget(name_database("GeneC")) || regulationFile.checkTarget(name_database("GeneD")) || regulationFile.checkTarget(name_database("GeneE")));
    EXPECT_FALSE(regulationFile.checkTarget(name_database("GeneA")));
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneB")).size(), 2);
    for(const auto& reg : regulationFile.target2regulations(name_database("GeneB"))){
        EXPECT_TRUE(std::get<1>(reg) == name_database("GeneB"));
        EXPECT_TRUE(name_database(std::get<0>(reg)) == "GeneA" || name_database(std::get<0>(reg)) == "GeneD");
        if(name_database(std::get<0>(reg)) == "GeneA"){
            EXPECT_EQ(std::get<2>(reg), 1.0);
        } else if(name_database(std::get<0>(reg)) == "GeneD"){
            EXPECT_EQ(std::get<2>(reg), 4.0);
        } else {
            // If this happens we are screwed!!
            EXPECT_TRUE(false);
        }
    }
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneC")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneC"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneC");
    EXPECT_EQ(name_database(std::get<0>(reg)), "GeneA");
    EXPECT_EQ(std::get<2>(reg), 2.0);
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneD")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneD"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneD");
    EXPECT_EQ(name_database(std::get<0>(reg)), "GeneC");
    EXPECT_EQ(std::get<2>(reg), 3.0);
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneE")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneE"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneE");
    EXPECT_EQ(name_database(std::get<0>(reg)), "GeneA");
    EXPECT_EQ(std::get<2>(reg), 5.0);
}

TEST(RegulationFile, MatrixNameDatabase) {
    unsigned int opts = DenseMatrixReader::NO_OPTIONS;
	opts |= DenseMatrixReader::READ_ROW_NAMES;
    opts |= DenseMatrixReader::READ_COL_NAMES;
    DenseMatrixReader reader;
    std::ifstream strm(TEST_DATA_PATH("RegulationFile_Matrix.txt"), std::ios::binary);
	DenseMatrix matrix = reader.read(strm, opts);
    
    MatrixNameDatabase name_database(&matrix);
    
    EXPECT_EQ(name_database.size(), 5);
    
    EXPECT_EQ(name_database(0), "GeneA");
    EXPECT_EQ(name_database("GeneA"), 0);
    
    EXPECT_EQ(name_database(1), "GeneB");
    EXPECT_EQ(name_database("GeneB"), 1);
    
    EXPECT_EQ(name_database(2), "GeneC");
    EXPECT_EQ(name_database("GeneC"), 2);
    
    EXPECT_EQ(name_database(3), "GeneD");
    EXPECT_EQ(name_database("GeneD"), 3);
    
    EXPECT_EQ(name_database(4), "GeneE");
    EXPECT_EQ(name_database("GeneE"), 4);
    
    std::unordered_set<size_t> tset;
    for(size_t i=0; i<5; ++i){
        tset.emplace(i);
    }
    
    RegulationFileParser<MatrixNameDatabase, double> parser(name_database, tset, TEST_DATA_PATH("RegulationFile.txt"), 0.0);
	RegulationFile<double>& regulationFile = parser.getRegulationFile();
    
    for(const auto& i : regulationFile.regulators()){
        EXPECT_TRUE(name_database(i) == "GeneA" || name_database(i) == "GeneC" || name_database(i) == "GeneD");
    }
    
    EXPECT_TRUE(regulationFile.checkRegulator(name_database("GeneA")) || regulationFile.checkRegulator(name_database("GeneC")) || regulationFile.checkRegulator(name_database("GeneD")));
    EXPECT_FALSE(regulationFile.checkRegulator(name_database("GeneB")));
    EXPECT_FALSE(regulationFile.checkRegulator(name_database("GeneE")));
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database("GeneA")).size(), 3);
    for(const auto& reg : regulationFile.regulator2regulations(name_database("GeneA"))){
        EXPECT_TRUE(std::get<0>(reg) == name_database("GeneA"));
        EXPECT_TRUE(name_database(std::get<1>(reg)) == "GeneB" || name_database(std::get<1>(reg)) == "GeneC" || name_database(std::get<1>(reg)) == "GeneE");
    }
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database("GeneC")).size(), 1);
    auto reg = regulationFile.regulator2regulations(name_database("GeneC"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneD");
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database("GeneD")).size(), 1);
    reg = regulationFile.regulator2regulations(name_database("GeneD"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneB");
    
    EXPECT_TRUE(regulationFile.checkTarget(name_database("GeneB")) || regulationFile.checkTarget(name_database("GeneC")) || regulationFile.checkTarget(name_database("GeneD")) || regulationFile.checkTarget(name_database("GeneE")));
    EXPECT_FALSE(regulationFile.checkTarget(name_database("GeneA")));
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneB")).size(), 2);
    for(const auto& reg : regulationFile.target2regulations(name_database("GeneB"))){
        EXPECT_TRUE(std::get<1>(reg) == name_database("GeneB"));
        EXPECT_TRUE(name_database(std::get<0>(reg)) == "GeneA" || name_database(std::get<0>(reg)) == "GeneD");
    }
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneC")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneC"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneC");
    EXPECT_EQ(name_database(std::get<0>(reg)), "GeneA");
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneD")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneD"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneD");
    EXPECT_EQ(name_database(std::get<0>(reg)), "GeneC");
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneE")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneE"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneE");
    EXPECT_EQ(name_database(std::get<0>(reg)), "GeneA");
}


TEST(RegulationFile, MatrixNameDatabaseMicro) {
    unsigned int opts = DenseMatrixReader::NO_OPTIONS;
	opts |= DenseMatrixReader::READ_ROW_NAMES;
    opts |= DenseMatrixReader::READ_COL_NAMES;
    DenseMatrixReader reader;
    std::ifstream strm(TEST_DATA_PATH("RegulationFile_Matrix.txt"), std::ios::binary);
    std::ifstream strm_micro(TEST_DATA_PATH("RegulationFile_MatrixMicro.txt"), std::ios::binary);

	DenseMatrix matrix = reader.read(strm, opts);
	DenseMatrix micro = reader.read(strm_micro,opts);
    
    MatrixNameDatabase name_database(&matrix);
    MatrixNameDatabase name_database_micro(&micro);
    
    EXPECT_EQ(name_database.size(), 5);
    
    EXPECT_EQ(name_database(0), "GeneA");
    EXPECT_EQ(name_database("GeneA"), 0);
    
    EXPECT_EQ(name_database(1), "GeneB");
    EXPECT_EQ(name_database("GeneB"), 1);
    
    EXPECT_EQ(name_database(2), "GeneC");
    EXPECT_EQ(name_database("GeneC"), 2);
    
    EXPECT_EQ(name_database(3), "GeneD");
    EXPECT_EQ(name_database("GeneD"), 3);
    
    EXPECT_EQ(name_database(4), "GeneE");
    EXPECT_EQ(name_database("GeneE"), 4);
    
    
    EXPECT_EQ(name_database_micro.size(), 4);
    
    EXPECT_EQ(name_database_micro(0), "miRNA1");
    EXPECT_EQ(name_database_micro("miRNA1"), 0);
    
    EXPECT_EQ(name_database_micro(1), "miRNA2");
    EXPECT_EQ(name_database_micro("miRNA2"), 1);
    
    EXPECT_EQ(name_database_micro(2), "miRNA3");
    EXPECT_EQ(name_database_micro("miRNA3"), 2);
    
    EXPECT_EQ(name_database_micro(3), "miRNA4");
    EXPECT_EQ(name_database_micro("miRNA4"), 3);

    
    std::unordered_set<size_t> tset;
    for(size_t i=0; i<5; ++i){
        tset.emplace(i);
    }    
    RegulationFileParser<MatrixNameDatabase, double> parser(name_database,name_database_micro, tset, TEST_DATA_PATH("RegulationFileMicro.txt"), 0.0);
	RegulationFile<double>& regulationFile = parser.getRegulationFile();
    
    for(const auto& i : regulationFile.regulators()){
        EXPECT_TRUE(name_database_micro(i) == "miRNA1" || name_database_micro(i) == "miRNA3" || name_database_micro(i) == "miRNA4");
    }
    
    EXPECT_TRUE(regulationFile.checkRegulator(name_database_micro("miRNA1")) || regulationFile.checkRegulator(name_database_micro("miRNA3")) || regulationFile.checkRegulator(name_database_micro("miRNA4")));
    EXPECT_FALSE(regulationFile.checkRegulator(name_database_micro("miRNA2")));
    EXPECT_FALSE(regulationFile.checkRegulator(name_database_micro("miRNA8")));
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database_micro("miRNA1")).size(), 2);
    for(const auto& reg : regulationFile.regulator2regulations(name_database_micro("miRNA1"))){
        EXPECT_TRUE(std::get<0>(reg) == name_database_micro("miRNA1"));
        EXPECT_TRUE(name_database(std::get<1>(reg)) == "GeneE" || name_database(std::get<1>(reg)) == "GeneB");
    }
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database_micro("miRNA3")).size(), 1);
    auto reg = regulationFile.regulator2regulations(name_database_micro("miRNA3"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneC");
    
    EXPECT_EQ(regulationFile.regulator2regulations(name_database_micro("miRNA4")).size(), 2);
    for(const auto& reg : regulationFile.regulator2regulations(name_database_micro("miRNA4"))){
        EXPECT_TRUE(std::get<0>(reg) == name_database_micro("miRNA4"));
        EXPECT_TRUE(name_database(std::get<1>(reg)) == "GeneD" || name_database(std::get<1>(reg)) == "GeneB");
    }
    
    EXPECT_TRUE(regulationFile.checkTarget(name_database("GeneB")) || regulationFile.checkTarget(name_database("GeneC")) || regulationFile.checkTarget(name_database("GeneD")) || regulationFile.checkTarget(name_database("GeneE")));
    EXPECT_FALSE(regulationFile.checkTarget(name_database("GeneA")));
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneC")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneC"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneC");
    EXPECT_EQ(name_database_micro(std::get<0>(reg)), "miRNA3");
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneD")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneD"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneD");
    EXPECT_EQ(name_database_micro(std::get<0>(reg)), "miRNA4");
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneE")).size(), 1);
    reg = regulationFile.target2regulations(name_database("GeneE"))[0];
    EXPECT_EQ(name_database(std::get<1>(reg)), "GeneE");
    EXPECT_EQ(name_database_micro(std::get<0>(reg)), "miRNA1");
    
    EXPECT_EQ(regulationFile.target2regulations(name_database("GeneB")).size(), 2);
    for (const auto& reg : regulationFile.target2regulations(name_database("GeneB"))){
      EXPECT_TRUE(name_database_micro(std::get<0>(reg)) == "miRNA1" || name_database_micro(std::get<0>(reg)) == "miRNA4");
      EXPECT_TRUE(std::get<1>(reg) == name_database("GeneB"));
      
    }
}
