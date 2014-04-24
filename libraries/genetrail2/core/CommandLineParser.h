/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef COMMANDLINE_PARSER_H
#define COMMANDLINE_PARSER_H

#include "macros.h"

#include <iostream>
#include <cstring>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace GeneTrail
{

	/**
     *	Wrapper for the boost program options library.
     */
	class GT2_EXPORT CommandLineParser
	{
		public:

			CommandLineParser(std::string description);

			~CommandLineParser() {};

			/**
             * This function adds a new option to the option map.
             * This option will now be used to parse the command line.
             *
             * @param option 
			 * @param message Help message for the specified option
             */
			void addOption(std::string option, std::string message);

			/**
             * This function adds a new option (with a defined type) to the option map.
             * This option will now be used to parse the command line.
             *
             * @param option 
             * @param message Help message for the specified option
             */
			template <typename T> void addTypedOption(std::string option, std::string message)
			{
				description_.add_options()
				(option.c_str(), boost::program_options::value<T>(), message.c_str());
			}

			/**
             * This function adds a new option (with a defined type and default value) to the option map.
             * This option will now be used to parse the command line.
             *
             * @param option 
             * @param default_value Default value for this option
             * @param message Help message for the specified option
             */
			template <typename T> void addDefaultOption(std::string option, T default_value, std::string message)
			{
				description_.add_options()
				(option.c_str(), boost::program_options::value<T>()->default_value(default_value), message.c_str());
			}

			/**
             * This function adds a new option (with implicit and default value) to the option map.
             * This option will now be used to parse the command line.
             *
             * @param option 
             * @param default_value Default value for this option
             * @param message Help message for the specified option
             */
			template <typename T> void addImplicitOption(std::string option, T default_value, T implicit_value, std::string message)
			{
				description_.add_options()
				(option.c_str(), boost::program_options::value<T>()->default_value(default_value)->implicit_value(implicit_value), message.c_str());
			}

			/**
             * This function checks if a given option was specified on the command line. 
             * The given value is saved as parameter. 
             *
             * @param option
             * @param parameter This variable saves the given value
             * @return Flag if the option was specified as command line argument.
             */
			template <typename T> bool getParameter(std::string option, T& parameter)
			{
				bool b;

				if((b=vm_.count(option.c_str())))
				{
					parameter = vm_[option.c_str()].as<T>();
				}

				return b;
			}

			/**
             * This function checks if a given option was specified on the command line. 
             *
             * @param option
             * @return Flag if the option was specified as command line argument.
             */
			bool checkParameter(std::string option);

			/**
             * This function parses the command line and checks for specified options. 
             *
             * @param argc Number of command line arguments
             * @param argv Command line arguments
             */
			void parse(int argc, char* argv[]);

			/**
             * Prints a message to the standard output stream containing the help messages of all options,
             */
			void printHelp();

		private:

			boost::program_options::variables_map vm_;
			boost::program_options::options_description description_;
			boost::program_options::positional_options_description p_;

	};
}

#endif

