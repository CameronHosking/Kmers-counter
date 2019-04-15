#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

std::vector<std::string> split(const std::string & stringToSplit, const std::string & delimiter);

std::vector<std::pair<std::string, std::string>> stringToListOfKeyValuePairs(const std::string &stringToParse, const std::string &keyValueDelimiter, const std::string &elementDelimiter);

char * getLocOfCmdLineArg(int argc, char ** argv, const std::string & command);

bool cmdLineArgumentFound(int argc, char ** argv, const std::string & command);

bool setToBoolIfFlagFound(int argc, char ** argv, const std::string & command, bool setTo, bool & value);

bool getCmdLineArgument(int argc, char ** argv, const std::string & command, std::string & value);

//simple getCmdLineArgument for types supported by iostream
template <typename T>
bool getCmdLineArgument(int argc, char **argv, const std::string & command, T &value)
{
	char* arg = getLocOfCmdLineArg(argc, argv, command);
	if (arg == nullptr)
		return false;
	else
	{
		std::stringstream ss;
		ss << arg;
		ss >> value;
		return true;
	}
}

