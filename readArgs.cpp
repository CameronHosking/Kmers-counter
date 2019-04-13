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

char * getLocOfCmdLineArg(int argc, char ** argv, const char * command)
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
			if (mismatchFound || index < strlen(command))
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

bool getCmdLineArgument(int argc, char ** argv, const char * command)
{
	return getLocOfCmdLineArg(argc, argv, command) != nullptr;
}

bool setToBoolIfFlagFound(int argc, char ** argv, const char * command, bool setTo, bool & value)
{
	if (getCmdLineArgument(argc, argv, command))
	{
		value = setTo;
		return true;
	}
	return false;
}

bool getCmdLineArgument(int argc, char ** argv, const char * command, std::string & value)
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


Args::Args(const std::string &filename)
{
	std::string line;
	std::ifstream myfile(filename, std::ifstream::in);
	if (!myfile.is_open())
	{
		std::cout << "could not open " << filename << " for reading" << std::endl;
		loaded = false;
		return;
	}
	std::vector<std::string> lines;
	lines.push_back(filename);
	while (std::getline(myfile, line))
	{
		lines.push_back(trim(line));
	}
	c = (int)lines.size();
	v = new char*[c];
	for (int i = 0; i < lines.size(); ++i)
	{
		size_t size = lines.at(i).size() + 1;
		v[i] = new char[size];
		strcpy_s(v[i], size, lines.at(i).c_str());
	}
	loaded = true;
}

Args::~Args()
{
	for (int i = 0; i < c; ++i)
	{
		delete[] v[i];
	}
}
