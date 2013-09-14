#include "commandline_parser.h"

using namespace GeneTrail;

CommandLineParser::CommandLineParser(std::string description)
:description_(description.c_str())
{
}

void CommandLineParser::addOption(std::string option, std::string message)
{
	description_.add_options()
        (option.c_str(), message.c_str());
}

bool CommandLineParser::checkParameter(std::string option)
{
    return vm_.count(option.c_str());
}

void CommandLineParser::parse(int argc, char* argv[])
{
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(description_).positional(p_).run(), vm_);
    boost::program_options::notify(vm_);
}

void CommandLineParser::printHelp()
{
	std::cout << description_ << std::endl;
}

