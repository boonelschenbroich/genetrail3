#ifndef SCORING_FILE_PARSER_H
#define SCORING_FILE_PARSER_H

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <stdlib.h>

#include "config.h"

namespace GeneTrail {

    typedef std::tuple<std::string, double> mytuple;

    struct increasing_compare {

        bool operator() (mytuple a, mytuple b) {
            return (std::get<1>(a) < std::get<1>(b));
        }
    };

    struct decreasing_compare {

        bool operator() (mytuple a, mytuple b) {
            return (std::get<1>(a) > std::get<1>(b));
        }
    };

    struct absolute_compare {

        bool operator() (mytuple a, mytuple b) {
            return (std::abs(std::get<1>(a)) > std::abs(std::get<1>(b)));
        }
    };

    class GT2_EXPORT ScoringFileParser {
    public:

        ScoringFileParser() {
        };

        ~ScoringFileParser() {
        };

        void readScoringFile(std::string file);

        void sortScoringFile(bool descending);
        void sortScoringFileAbsolute();
        void sortScoringFileIncreasingly();
        void sortScoringFileDecreasingly();

        std::vector<std::tuple<std::string, double> > getScores();
        std::vector<std::tuple<std::string, double> > getSortedScores(bool descending);
        std::vector<std::tuple<std::string, double> > getIncreasinglySortedScores();
        std::vector<std::tuple<std::string, double> > getDecreasinglySortedScores();
        std::vector<std::tuple<std::string, double> > getAbsoluteSortedScores();

        std::vector<std::string> getFirstK(int k);
        std::vector<std::string> getFirstKInSet(int k, std::set<std::string> myset);
        std::vector<std::string> getAllInSet(std::set<std::string> myset, bool descending);
        std::vector<std::string> getAllInSet(std::set<std::string> myset);

    protected:

        std::vector<std::tuple<std::string, double> > scores_;
        std::vector<std::tuple<std::string, double> > sorted_scores_;
    };
}

#endif

