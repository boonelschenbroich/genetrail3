#include <vector>

#include "Path.h"

using namespace GeneTrail;

void Path::writeToSIFFile(const std::string& name) {
    std::ofstream myfile;
    myfile.open(name);
    for (unsigned int i = 1; i < this->length(); ++i) {
        auto regulation = this->regulations_.find(this->identifier_[i-1] + this->identifier_[i]);
        if(regulation != this->regulations_.end()){
            myfile << this->identifier_[i-1] << "\tpp\t" << this->identifier_[i] << std::endl;
        }else{
            myfile << this->identifier_[i-1] << "\t" << this->getRegulation(this->identifier_[i-1],this->identifier_[i]) << "\t" << this->identifier_[i] << std::endl;
        }
    }
    myfile.close();
}

void Path::addVertex(const std::string& vertex) {
    this->identifier_.push_back(vertex);
}

std::string Path::getVertex(const int& position) {
    return this->identifier_[position];
}

void Path::addRegulation(const std::string& v1, const std::string& v2, const std::string& regulation) {
    this->regulations_[v1 + v2] = regulation;
}

std::string Path::getRegulation(const std::string& v1, const std::string& v2) {
    return this->regulations_[v1 + v2];
}

int Path::length() {
    return this->identifier_.size();
}

void Path::setRunningSum(int r){
	this->runningSum_ = r;
}

int Path::runningSum() {
    return this->runningSum_;
}

void Path::setPValue(double p){
	this->pValue_ = p;
}

double Path::pValue() {
    return this->pValue_;
}
