
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <chrono>
#include <cstring>
#include <stdio.h>
#include <bitset>
#include "ryu/ryu2_header_only.hpp"
#include "ska_sort.hpp"
#include "readArgs.h"

#define SIG_FIGS 10
//#define DO_TIMING

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

	bool read(std::string &s, size_t charsToRead)
	{
		s.resize(charsToRead);

		size_t charsRead = fread(&s[0] ,1 , charsToRead, fp);
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

class OutPutter
{
	std::ostream *out;
	bool readingFile;
public:
	OutPutter()
		:out(&std::cout), readingFile(false)
	{}

	void open(const std::string &filename)
	{
		if (readingFile)
			delete out;
		out = new std::ofstream(filename, std::ostream::trunc | std::ostream::binary);
		readingFile = true;
	}

	std::ostream &getStream()
	{
		return *out;
	}

	~OutPutter()
	{
		if (readingFile)
			delete out;
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
	bool allKs;

	Parameters(int argc, char** argv)
	{
		outputZeros = cmdLineArgumentFound(argc, argv, "outputZeros")|cmdLineArgumentFound(argc, argv, "z");
		split = cmdLineArgumentFound(argc, argv, "split")|cmdLineArgumentFound(argc, argv, "s");
		allKs = cmdLineArgumentFound(argc, argv, "allUpToK") | cmdLineArgumentFound(argc, argv, "a");

		if (!(getCmdLineArgument(argc, argv, "k", k) || getCmdLineArgument(argc, argv, "K", k)) || k<1 || k>15)
		{
			std::cout << "Usage KmersCounter K=<integer> [options]" << std::endl;
			std::cout << "The integer provided for K must be between 1 and 15 inclusive." << std::endl;
			std::cout << "Input is streamed from stdin unless input file is specified via input=<filename>|i=<filename>." << std::endl;
			std::cout << "Output is to stdout unless output=<filename>|o=<filename> is set." << std::endl;
			std::cout << "-split|-s splits k-mer counts between id lines (lines starting with >)." << std::endl;
			std::cout << "-outputZeros|-z outputs 0 frequencies otherwise that k-mer is skipped in output." << std::endl;
			std::cout << "-allUpToK|-a If this flag is present all Kmers of length up to and including K are output." << std::endl;
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
void output(int index, char * kmer, uint8_t charLocation, const uint32_t results[], double denominator, std::ostream &os)
{
	for (int i = 0; i < 4; ++i)
	{
		kmer[charLocation] = chars[i];
		output<K - 1, OutputZeros>((index << 2) + i, kmer,charLocation+1, results, denominator, os);
	}
}

//base case prints out the kmer and associated frequency
template<>
void output<0,true>(int index, char * kmer, uint8_t charLocation, const uint32_t results[], double denominator, std::ostream &os)
{
	ryu::d2exp_buffered_n(results[index]*denominator, SIG_FIGS -1, &kmer[charLocation+1]);
	os.write(kmer, charLocation + SIG_FIGS + 7);
}

template<>
void output<0, false>(int index, char * kmer, uint8_t charLocation, const uint32_t results[], double denominator, std::ostream &os)
{
	int count = results[index];

	if (count != 0)
	{
		ryu::d2exp_buffered_n(count*denominator, SIG_FIGS - 1, &kmer[charLocation + 1]);
		os.write(kmer, charLocation + SIG_FIGS + 7);
	}
		
}

//prepares an outputter
void prepareOutputStream(const Parameters &p, std::string id, int idsFound, OutPutter &out)
{
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
		out.open(outFileName);
		std::cout << "outputting Kmers to: " << outFileName << std::endl;
	}
	if (p.split)
	{
		out.getStream() << id << std::endl;
	}
}

//outputs the kmers frequency, if split is true it outputs the id as the first line,
//if output zeros is true then output kmers with zero counts
template<unsigned int K>
void output(const Parameters &p, const uint32_t results[], size_t total, std::ostream &out)
{
	//output string
	//+8 is from: 1 tab, the decimal point, the e-XX, the newline character, and the null character
	char baseKmer[K+ SIG_FIGS + 8];
	baseKmer[K] = '\t';
	baseKmer[K + SIG_FIGS + 6] = '\n';
	baseKmer[K + SIG_FIGS + 7] = '\0';
	if (p.outputZeros)
	{
		output<K, true>(0, baseKmer,0, results, 1.0 / (double)total, out);
	}
	else
	{
		output<K, false>(0, baseKmer,0, results, 1.0 / (double)total, out);
	}
}

//template<unsigned int K, unsigned int S>
//class outputMultiple
//
template<unsigned int K>
class outputMultiple
{
public:

	template<unsigned int S=1>
	static void call(const Parameters &p, const std::unique_ptr<uint32_t[]> results[K], size_t totals[K], std::ostream &out)
	{
		output<S>(p, results[S-1].get(), totals[S - 1], out);
		call<S+1>(p, results, totals, out);
	}

	template<>
	static void call<K>(const Parameters &p, const std::unique_ptr<uint32_t[]> results[K], size_t totals[K], std::ostream &out)
	{
		output<K>(p, results[K - 1].get(), totals[K - 1], out);
	}


};

//call at the start of an id line
//set the id to the contents of this line
//return the index at the end of the line
size_t newID(const Parameters &p, size_t i, std::string &id, std::string &s, InputReader &in)
{

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
void prepareMultipleKResults(std::unique_ptr<uint32_t[]> *results, size_t *totals)
{
	size_t kmers = 1LL << K * 2;
	size_t lowerKindex = 0;
	for (size_t i = 0; i < kmers; i += 4)
	{
		results[K - 2][lowerKindex] += results[K - 1][i] + results[K - 1][i + 1] + results[K - 1][i + 2] + results[K - 1][i + 3];
		lowerKindex++;
	}
	totals[K - 2] += totals[K - 1];
	prepareMultipleKResults<K - 1>(results, totals);
}
template <>
void prepareMultipleKResults<1>(std::unique_ptr<uint32_t[]> *results, size_t *totals)
{}

template <unsigned int K>
inline void countSmallerKmers(std::unique_ptr<uint32_t[]> resultsForall[K], size_t totals[K], size_t count, uint64_t currentString)
{
	currentString >>= 2;
	size_t i = (count < K - 1) ? count : K - 1;
	uint64_t kmers = 1LL << i * 2;
	uint64_t mask = kmers - 1;
	for (; i > 0; i--)
	{
		resultsForall[i - 1][mask&currentString]++;
		totals[i - 1]++;
		mask >>= 2;
	}
}

template <unsigned int K, bool allUpToK>
void countKmers(InputReader &input,const Parameters &p)
{
#ifdef DO_TIMING	
	auto start = std::chrono::high_resolution_clock::now();
#endif
	//allocate memory for the 4^K different kmers
	uint64_t kmers = 1LL << K * 2;
	std::unique_ptr<uint32_t[]> results;
	std::unique_ptr<uint32_t[]> resultsForall[K];
	//creates a bitmask that we can slap on a coded string to grab the last K characters
	uint64_t mask = kmers - 1;;
	size_t totals[K];
	size_t total = 0;

	if (allUpToK)
	{
		for (int i = 0; i < K; ++i)
		{
			uint64_t kmers = 1LL << (i+1) * 2;
			totals[i] = 0;
			resultsForall[i] = std::unique_ptr<uint32_t[]>(new uint32_t[kmers]());
		}
	}
	else
	{
		results = std::unique_ptr<uint32_t[]>(new uint32_t[kmers]());
	}

	uint64_t currentString = 0;//encodes the current string;
	size_t count = 0;
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
					//convert the character to a 2bit number and add it to the end of the current string
					currentString += charTo2bit(c);
					count++;
					
					if (count >= K)
					{
						//increment the number of counts for this particular string
						(allUpToK ? resultsForall[K - 1] : results)[currentString&mask]++;
						(allUpToK ? totals[K - 1] : total) += 1;
					}
					currentString <<= 2;
				}//otherwise it is an unknown so reset the count
				else
				{
					//count the shorter kmers
					if (allUpToK && count > 0)
					{
						countSmallerKmers<K>(resultsForall, totals, count, currentString);
					}
					count = 0;
				}
			}
			else if (c == '>')
			{
				//count the shorter kmers
				if (allUpToK && count > 0)
				{
					countSmallerKmers<K>(resultsForall, totals, count, currentString);
				}
				//if split is true then we ouput individual counts for each id found
				if (p.split&&(allUpToK?totals[0]:total) > 0)
				{
					OutPutter out;
					prepareOutputStream(p, id, idsFound, out);
					if (allUpToK)
					{
						prepareMultipleKResults<K>(resultsForall, totals);
						outputMultiple<K>::call(p, resultsForall, totals, out.getStream());
						for (int Kindex = 0; Kindex < K; ++Kindex)
						{
							totals[Kindex] = 0;
							int kmers = (1LL << (Kindex + 1) * 2);
							memset(resultsForall[Kindex].get(), 0, kmers * sizeof(int));
						}
					}
					else
					{
						output<K>(p, results.get(), total, out.getStream());
						//zero the results
						memset(results.get(), 0, kmers * sizeof(int));
						total = 0;
					}
				}	

				i = newID(p,  i, id, s, input);
				idsFound++;
				count = 0;
			}
			//other characters are ignored eg \n\r etc
		}
	}
	countSmallerKmers<K>(resultsForall, totals, count, currentString);
#ifdef	DO_TIMING
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
	start = std::chrono::high_resolution_clock::now();
#endif
	OutPutter out;
	prepareOutputStream(p, id, idsFound, out);
	if (allUpToK)
	{
		prepareMultipleKResults<K>(resultsForall, totals);
		outputMultiple<K>::call(p, resultsForall, totals, out.getStream());
	}
	else{
		output<K>(p, results.get(), total, out.getStream());
	}
	
#ifdef DO_TIMING 
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl; 
#endif
}

//this wraps countKmers so that we can compile in the boolean of whether or not to count all up to kmers
template <unsigned int K>
void countMultiKmers(InputReader &input, const Parameters &p)
{
	if (p.allKs)
	{
		countKmers<K, true>(input, p);
	}
	else
	{
		countKmers<K, false>(input, p);
	}
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
	case 1: countMultiKmers<1>(s, p); break;
	case 2: countMultiKmers<2>(s, p); break;
	case 3: countMultiKmers<3>(s, p); break;
	case 4: countMultiKmers<4>(s, p); break;
	case 5: countMultiKmers<5>(s, p); break;
	case 6: countMultiKmers<6>(s, p); break;
	case 7: countMultiKmers<7>(s, p); break;
	case 8: countMultiKmers<8>(s, p); break;
	case 9: countMultiKmers<9>(s, p); break;
	case 10: countMultiKmers<10>(s, p); break;
	case 11: countMultiKmers<11>(s, p); break;
	case 12: countMultiKmers<12>(s, p); break;
	case 13: countMultiKmers<13>(s, p); break;
	case 14: countMultiKmers<14>(s, p); break;
	case 15: countMultiKmers<15>(s, p); break;
	}
	return 0;
}
