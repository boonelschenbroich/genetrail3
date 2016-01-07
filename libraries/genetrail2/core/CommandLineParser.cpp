/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013-2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2013 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#include "CommandLineParser.h"

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

