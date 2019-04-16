
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <chrono>
#include <cstring>
#include <stdio.h>

#include "readArgs.h"
class InputReader
{
	std::FILE *fp;
	bool readingFile;
public:
	InputReader()
		:fp(stdin), readingFile(false)
	{
	}

	void open(const std::string &filename)
	{
		if (readingFile)
			fclose(fp);
#ifdef _WIN32
		fopen_s(&fp, filename.c_str(), "rb");
#else // _WIN32
		fp = fopen(filename.c_str(), "rb");
#endif
		readingFile = true;
		if (fp==NULL)
		{
			std::cout << "could not open \"" << filename << "\"" << std::endl;
			exit(1);
		}
	}

	bool read(std::string &s, int64_t charsToRead)
	{
		s.resize(charsToRead);

		int charsRead = fread(&s[0] ,1 , charsToRead, fp);
		//we reached the end of the input
		if (charsRead < charsToRead)
		{
			s.resize(charsRead);
			if (charsRead == 0)
				return false;
		}
		return true;
	}
	~InputReader()
	{
		if(readingFile)
			fclose(fp);
	}
};

struct Parameters
{
	std::string inputFile;
	std::string outputFile;
	bool split;
	bool outputZeros;
	bool outputToFile;
	unsigned int k;
	Parameters(int argc, char** argv)
	{
		outputZeros = cmdLineArgumentFound(argc, argv, "outputZeros")|cmdLineArgumentFound(argc, argv, "z");
		split = cmdLineArgumentFound(argc, argv, "split")|cmdLineArgumentFound(argc, argv, "s");

		if (!(getCmdLineArgument(argc, argv, "k", k) || getCmdLineArgument(argc, argv, "K", k)) || k<1 || k>15)
		{
			std::cout << "Usage KmersCounter K=<integer> [options]" << std::endl;
			std::cout << "The integer provided for K must be between 1 and 15 inclusive." << std::endl;
			std::cout << "Input is streamed from stdin unless input file is specified via input=<filename>|i=<filename>." << std::endl;
			std::cout << "Output is to stdout unless output=<filename>|o=<filename> is set." << std::endl;
			std::cout << "-split|-s splits k-mer counts between id lines (lines starting with >)." << std::endl;
			std::cout << "-outputZeros|-z outputs 0 frequencies otherwise that k-mer is skipped in output." << std::endl;
			std::cout << "See README for more information" << std::endl;
			k = 0;
		}

		getCmdLineArgument(argc, argv, "input", inputFile);
		getCmdLineArgument(argc, argv, "i", inputFile);
		getCmdLineArgument(argc, argv, "output", outputFile);
		getCmdLineArgument(argc, argv, "o", outputFile);

		outputToFile = (outputFile != std::string());
	}
};

constexpr char chars[4] = { 'A','C','G','T' };

template <unsigned int K,bool OutputZeros>
void output(int index, std::string &kmer, const uint32_t results[], double denominator, std::ostream &os)
{
	for (int i = 0; i < 4; ++i)
	{
		kmer[kmer.length() - K] = chars[i];
		output<K - 1, OutputZeros>((index << 2) + i, kmer, results, denominator, os);
	}
}

//base case prints out the kmer and associated frequency
template<>
void output<0,true>(int index, std::string &kmer, const uint32_t results[], double denominator, std::ostream &os)
{
	os << kmer << '\t' << results[index] *denominator << '\n';
}

template<>
void output<0, false>(int index, std::string &kmer, const uint32_t results[], double denominator, std::ostream &os)
{
	int count = results[index];
	if (count != 0)
		os << kmer << '\t' << count*denominator << '\n';
}

//outputs the kmers frequency, if split is true it outputs the id as the first line,
//if output zeros is true then output kmers with zero counts
template<unsigned int K>
void output(Parameters &p, std::string id, const uint32_t results[], size_t total, int idsFound)
{
	//output string
	std::string baseKmer = std::string(K, ' ');
	std::ostream *out = &std::cout;
	if (p.outputToFile)
	{
		std::string outFileName = p.outputFile;
		if (p.split)
		{
			size_t pos = outFileName.find_first_of('%', 0);
			if (pos == std::string::npos)
			{
				outFileName += std::to_string(idsFound);
			}
			else
			{
				outFileName = outFileName.substr(0, pos) + std::to_string(idsFound) + outFileName.substr(pos + 1);
			}
		}
		out = new std::ofstream(outFileName, std::ostream::trunc| std::ostream::binary);
		std::cout << "outputting Kmers to: " << outFileName << std::endl;
	}
	if (p.split) {
		*out << id << std::endl;
	}
	if (p.outputZeros)
	{
		output<K, true>(0, baseKmer, results, 1.0 / (double)total, *out);
	}
	else
	{
		output<K, false>(0, baseKmer, results, 1.0 / (double)total, *out);
	}
	if (p.outputFile != std::string())
	{
		delete out;
	}
}

//call at the start of an id line, splits and outputs is split is true
//set the id to the contents of this line
//return the index at the end of the line
template<unsigned int K>
size_t newID(Parameters &p, size_t total, size_t i, std::string &id, uint32_t results[], std::string &s, InputReader &in, int idsFound)
{
	uint64_t kmers = 1LL << K * 2;
	//if split is true then we ouput individual counts for each id found
	if (p.split && total > 0)
	{
		output<K>(p, id, results, total, idsFound);
		//zero the results
		memset(results, 0, kmers*sizeof(int));
	}
	size_t endOfLine;
	id = "";
	//this is a while loop in case somehow the id length is longer than the read buffer we're using
	while((endOfLine = s.find_first_of('\n', i)) == std::string::npos)
	{
		id += s.substr(i);
		in.read(s, 1000000);
		i = 0;
		endOfLine = s.find_first_of('\n', i);
	}

	id += s.substr(i, endOfLine - i);
	if (p.outputToFile)
		std::cout << "counting Kmers in:" << id << std::endl;
	return endOfLine;
}

//converts Aa=0, Cc=1, Gg=2, Tt=3, some other value between 0 and 3 for all other characters
inline uint8_t charTo2bit(char c)
{
	//this v represents the character ACGT or acgt with a value of 0x00 for Aa, 0x01 for Cc, 0x11 for Gg, and 0x10 for Tt
	uint8_t v = (c & 7) >> 1;
	//xoring with itself bitshifted right by 1 doesnt effect A or C but flips the small bit of G and T
	return v ^ (v >> 1);
}

template <unsigned int K>
void countKmers(InputReader &input, Parameters &p)
{
	auto start = std::chrono::high_resolution_clock::now();
	//allocate memory for the 4^K different kmers
	uint64_t kmers = 1LL << K * 2;
	std::unique_ptr<uint32_t[]> results = std::unique_ptr<uint32_t[]>(new uint32_t[kmers]());
	//creates a bitmask that we can slap on a coded string to grab the last K characters
	uint64_t mask = kmers - 1;
	size_t count = 0;
	uint64_t currentString = 0;//encodes the current string;
	size_t total = 0;
	std::string id;
	int idsFound = 0;
	std::string s;
	size_t i = 0;
	while (input.read(s, 1000000))
	{
		for (i=0; i < s.length(); i++)
		{
			char c = s[i];

			//all alpha characters
			if (c & 64)
			{
				//if the character is any of upper or lowercase ACGT
				char cLast5bits = c & 31;
				if (cLast5bits == ('A' & 31) || cLast5bits == ('C' & 31) || cLast5bits == ('G' & 31) || cLast5bits == ('T' & 31))
				{
					count++;
					//convert the character to a 2bit number and add it to the end of the current string
					currentString += charTo2bit(c);
					if (count >= K)
					{
						//increment the number of counts for this particular string
						results[currentString&mask]++;
						total += 1;
					}
					currentString <<= 2;
				}//otherwise it is an unknown so reset the count
				else
				{
					count = 0;
				}
			}
			else if (c == '>')
			{
				i = newID<K>(p, total, i, id, results.get(), s, input, idsFound);
				idsFound++;
				if (p.split)
				{
					total = 0;
				}
				count = 0;
			}
			//other characters are ignored eg \n\r etc
		}
	}
	if (p.outputToFile)
		std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
	start = std::chrono::high_resolution_clock::now();
	
	output<K>(p, id, results.get(), total, idsFound);
	if (p.outputToFile)
		std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
}

int main(int argc, char** argv)
{
	std::string inputFile;
	InputReader s;
	Parameters p(argc, argv);
	if (p.inputFile!=std::string())
	{
		s.open(p.inputFile);
	}
	
	//this allows us to select a K at runtime even though it is a compile time constant
	switch (p.k)
	{
	case 1: countKmers<1>(s, p); break;
	case 2: countKmers<2>(s, p); break;
	case 3: countKmers<3>(s, p); break;
	case 4: countKmers<4>(s, p); break;
	case 5: countKmers<5>(s, p); break;
	case 6: countKmers<6>(s, p); break;
	case 7: countKmers<7>(s, p); break;
	case 8: countKmers<8>(s, p); break;
	case 9: countKmers<9>(s, p); break;
	case 10: countKmers<10>(s, p); break;
	case 11: countKmers<11>(s, p); break;
	case 12: countKmers<12>(s, p); break;
	case 13: countKmers<13>(s, p); break;
	case 14: countKmers<14>(s, p); break;
	case 15: countKmers<15>(s, p); break;
	}
	return 0;
}