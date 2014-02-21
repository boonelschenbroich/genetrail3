#include "ScoringFileParser.h"

using namespace GeneTrail;

void ScoringFileParser::readScoringFile(std::string file) {
    std::ifstream infile(file);
    std::string line;

    while (std::getline(infile, line)) {
        if (boost::starts_with(line, "Scores")) {
            continue;
        }

        std::vector<std::string> sline;
        boost::split(sline, line, boost::is_any_of("="));

        if (sline.size() == 2) {
            boost::trim(sline[0]);
            boost::trim(sline[1]);
            std::tuple<std::string, double> tpl(sline[0], boost::lexical_cast<double>(sline[1]));
            scores_.push_back(tpl);
        }
    }
}

void ScoringFileParser::sortScoringFile(bool descending) {
    sorted_scores_ = scores_;

    if (descending) {
        std::sort(sorted_scores_.begin(), sorted_scores_.end(), descending_compare());
    } else {
        std::sort(sorted_scores_.begin(), sorted_scores_.end(), decreasing_compare());
    }
}

void ScoringFileParser::sortScoringFileAbsolute() {
    sorted_scores_ = scores_;
    std::sort(sorted_scores_.begin(), sorted_scores_.end(), absolute_compare());
}

void ScoringFileParser::sortScoringFileDescending() {
    sortScoringFile(true);
}

void ScoringFileParser::sortScoringFileDecreasing() {
    sortScoringFile(false);
}

std::vector<std::tuple<std::string, double> > ScoringFileParser::getScores() {
    return scores_;
}

std::vector<std::tuple<std::string, double> > ScoringFileParser::getSortedScores(bool descending) {
    sortScoringFile(descending);
    return sorted_scores_;
}

std::vector<std::tuple<std::string, double> > ScoringFileParser::getDescendinglySortedScores() {
    return getSortedScores(true);
}

std::vector<std::tuple<std::string, double> > ScoringFileParser::getDecreasinglySortedScores() {
    return getSortedScores(false);
}

std::vector<std::tuple<std::string, double> > ScoringFileParser::getAbsoluteSortedScores() {
    sortScoringFileAbsolute();
    return sorted_scores_;
}

std::vector<std::string> ScoringFileParser::getFirstK(int k) {
    std::vector<std::string> firstK;

    for (int i = 0; i < k; ++i) {
        firstK.push_back(std::get<0>(sorted_scores_[i]));
    }

    return firstK;
}

std::vector<std::string> ScoringFileParser::getFirstKInSet(int k, std::set<std::string> myset) {
    std::vector<std::string> firstK;

    for (int i = 0; i < k;) {
        std::set<std::string>::iterator setIt;
        setIt = myset.find(std::get<0>(sorted_scores_[i]));

        if (setIt != myset.end()) {
            ++i;
            firstK.push_back(std::get<0>(sorted_scores_[i]));
        }
    }

    return firstK;
}

std::vector<std::string> ScoringFileParser::getAllInSet(std::set<std::string> myset, bool descending) {
    if (sorted_scores_.size() == 0) {
        getSortedScores(descending);
    }

    std::vector<std::string> all;

    for (auto& elem : sorted_scores_) {
        std::set<std::string>::iterator setIt;
        setIt = myset.find(std::get<0>(elem));

        if (setIt != myset.end()) {
            all.push_back(std::get<0>(elem));
        }
    }

    return all;
}

std::vector<std::string> ScoringFileParser::getAllInSet(std::set<std::string> myset) {
    std::vector<std::string> all;

    for (auto& elem : sorted_scores_) {
        std::set<std::string>::iterator setIt;
        setIt = myset.find(std::get<0>(elem));

        if (setIt != myset.end()) {
            all.push_back(std::get<0>(elem));
        }
    }

    return all;
}

