#include "readArgs.h"
#include <fstream>
std::vector<std::string> split(const std::string & stringToSplit, const std::string & delimiter)
{
	//return the empty vector if the string is empty
	if (stringToSplit.empty())
		return std::vector<std::string>();

	std::vector<std::string> results;
	size_t start = 0;
	size_t end = stringToSplit.find(delimiter);
	while (end != std::string::npos)
	{
		results.push_back(stringToSplit.substr(start, end - start));
		start = end + delimiter.size();
		end = stringToSplit.find(delimiter, start);
	}

	results.push_back(stringToSplit.substr(start));
	return results;
}

std::vector<std::pair<std::string, std::string>> stringToListOfKeyValuePairs(const std::string & stringToParse, const std::string & keyValueDelimiter, const std::string & elementDelimiter)
{

	std::vector<std::string> elements = split(stringToParse, elementDelimiter);

	std::vector<std::pair<std::string, std::string>> keyValuePairs;

	for (std::string element : elements)
	{
		size_t split = element.find(keyValueDelimiter);
		if (split == std::string::npos)
		{
			std::cout << "element \"" << element << "\" is  missing the key value delimiter \"" << keyValueDelimiter << "\" skipping" << std::endl;
		}
		else
		{
			keyValuePairs.push_back(std::make_pair(element.substr(0, split), element.substr(split + elementDelimiter.size())));
		}
	}
	return keyValuePairs;
}

char * getLocOfCmdLineArg(int argc, char ** argv, const std::string & command)
{
	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			//ignores leading dashes
			size_t loc = 0;
			while (argv[i][loc] == '-')
			{
				loc++;
			}
			char * arg = argv[i] + loc;

			size_t index = 0;
			bool mismatchFound = false;
			while (command[index] != '\0'&&arg[index] != '\0')
			{
				if (command[index] != arg[index])
				{
					mismatchFound = true;
					break;
				}
				index++;
			}
			//the characters don't match or the argument is shorter than the command
			if (mismatchFound || index < command.length())
				continue;
			if (arg[index] == '=')
				index++;
			else if (arg[index] != '\0')
			{
				//if this is a flag without a value (no = ) this must be the end of the argument
				//otherwise this just happens to be an argument where the first characters match the flag
				continue;
			}
			return arg + index;
		}
	}
	return nullptr;
}

bool cmdLineArgumentFound(int argc, char ** argv, const std::string & command)
{
	return getLocOfCmdLineArg(argc, argv, command) != nullptr;
}

bool setToBoolIfFlagFound(int argc, char ** argv, const std::string & command, bool setTo, bool & value)
{
	if (cmdLineArgumentFound(argc, argv, command))
	{
		value = setTo;
		return true;
	}
	return false;
}

bool getCmdLineArgument(int argc, char ** argv, const std::string & command, std::string & value)
{
	char* arg = getLocOfCmdLineArg(argc, argv, command);
	if (arg == nullptr)
		return false;
	else
	{
		value = std::string(arg);
		return true;
	}
}

std::string trim(std::string s)
{
	size_t first = s.find_first_not_of(" \n\r\t");
	if (first == std::string::npos)
		return "";
	size_t last = s.find_last_not_of(" \n\r\t");
	return s.substr(first, (last - first + 1));
}