/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <vector>

#include "Path.h"

using namespace GeneTrail;

void Path::writeToSIFFile(const std::string& name) {
    std::ofstream myfile;
    myfile.open(name);
    for (int i = 1; i < this->length(); ++i) {
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
