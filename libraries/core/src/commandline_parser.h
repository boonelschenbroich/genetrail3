#ifndef COMMANDLINE_PARSER_H
#define COMMANDLINE_PARSER_H

#include <iostream>
#include <cstring>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace GeneTrail
{

	class CommandLineParser
	{
		public:

			CommandLineParser(std::string description);		

			~CommandLineParser() {};

			void addOption(std::string option, std::string message);
			
			template <typename T> void addTypedOption(std::string option, std::string message)
            {
                description_.add_options()
                    (option.c_str(), boost::program_options::value<T>(), message.c_str());
            }

			template <typename T> void addDefaultOption(std::string option, T default_value, std::string message)
			{
				description_.add_options()
        			(option.c_str(), boost::program_options::value<T>()->default_value(default_value), message.c_str());
			}

			template <typename T> void addImplicitOption(std::string option, T default_value, T implicit_value, std::string message)
            {
                description_.add_options()
                    (option.c_str(), boost::program_options::value<T>()->default_value(default_value)->implicit_value(implicit_value), message.c_str());
            }

			template <typename T> bool getParameter(std::string option, T& parameter)
			{
				bool b;
				if((b=vm_.count(option.c_str())))
				{
        			parameter = vm_[option.c_str()].as<T>();
				}
				return b;
			}
			
			bool checkParameter(std::string option);

			void parse(int argc, char* argv[]);

			void printHelp();

		private:
			
			boost::program_options::variables_map vm_;
			boost::program_options::options_description description_;
			boost::program_options::positional_options_description p_;
	
	};
}

#endif

