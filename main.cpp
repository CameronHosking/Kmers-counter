
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <chrono>

#include "readArgs.h"

std::string get_file_contents(const std::string &filename, bool verbose = false)
{
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	if (in)
	{
		std::string contents;
		in.seekg(0, std::ios::end);
		contents.resize(in.tellg());
		if (verbose)
		{
			std::cout << "reading in " << contents.size() << " Bytes from " << filename << std::endl;
		}
		in.seekg(0, std::ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	throw(errno);
}

constexpr char chars[4] = { 'A','C','G','T' };

template <unsigned int K,bool OutputZeros>
void output(int index, std::string &kmer, const std::unique_ptr<int[]> &results, double denominator, std::ostream &os)
{
	for (int i = 0; i < 4; ++i)
	{
		kmer[kmer.length() - K] = chars[i];
		output<K - 1, OutputZeros>((index << 2) + i, kmer, results, denominator, os);
	}
}

template<>
void output<0,true>(int index, std::string &kmer, const std::unique_ptr<int[]> &results, double denominator, std::ostream &os)
{
	os << kmer << '\t' << results[index] *denominator << '\n';
}

template<>
void output<0, false>(int index, std::string &kmer, const std::unique_ptr<int[]> &results, double denominator, std::ostream &os)
{
	int count = results[index];
	if (count != 0)
		os << kmer << '\t' << count*denominator << '\n';
}

template <unsigned int K>
void countKmers(const std::string &s, bool outputZeros)
{
	auto start = std::chrono::high_resolution_clock::now();
	//allocate memory for the 4^K different kmers
	uint64_t kmers = 1LL << K * 2;
	std::unique_ptr<int[]> results = std::unique_ptr<int[]>(new int[kmers]);
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
	start = std::chrono::high_resolution_clock::now();
	//creates a bitmask that we can slap on a coded string to grab the last K characters
	uint64_t mask = kmers - 1;
	size_t count = 0;
	uint64_t currentString = 0;//encodes the current string;
	size_t total = 0;
	std::string id;
	size_t i = 0;
	for (; i < s.length(); i++)
	{
		count++;
		switch (s[i])
		{
		case 'A':
		case 'a':
			currentString += 0;
			break;
		case 'C':
		case 'c':
			currentString += 1;
			break;
		case 'G':
		case 'g':
			currentString += 2;
			break;
		case 'T':
		case 't':
			currentString += 3;
			break;
		case '\n':
		case '\r':
			count--;
			continue;
		case '>':
		{
			if (total > 0)
			{
				std::string outFileName = "output" + std::to_string(i) + ".txt";
				std::ofstream out(outFileName, std::ostream::trunc | std::ostream::binary);
				std::cout << "outputting Kmers to: "<< outFileName << std::endl;
				out << id << std::endl;
				if (outputZeros)
				{
					output<K, true>(0, std::string(K, ' '), results, 1.0 / (double)total, out);
				}
				else
				{
					output<K, false>(0, std::string(K, ' '), results, 1.0 / (double)total, out);
				}
				out.close();
				total = 0;
				count = 0;
			}
			//zero the results
			memset(results.get(), 0, kmers*sizeof(int));
			size_t endOfLine = s.find_first_of('\n', i);
			id = s.substr(i, endOfLine - i);
			std::cout << "counting Kmers in:" << id << std::endl;
			i = endOfLine;
		}
		break;
		default:// the N character
			count = 0;
			break;
		}
		if (count >= K)
		{
			results[currentString&mask]++;
			total += 1;
		}
		currentString <<= 2;
	}
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
	start = std::chrono::high_resolution_clock::now();
	
	//output string
	std::string outFileName = "output" + std::to_string(i) + ".txt";
	std::ofstream out(outFileName, std::ostream::trunc | std::ostream::binary);
	std::cout << "outputting Kmers to: " << outFileName << std::endl;
	out << id << std::endl;
	if (outputZeros)
	{
		output<K, true>(0, std::string(K, ' '), results, 1.0 / (double)total, out);
	}
	else
	{
		output<K, false>(0, std::string(K, ' '), results, 1.0 / (double)total, out);
	}
	out.close();

	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
}

int main(int argc, char** argv)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::string inputFile;
	int k;
	if (!getCmdLineArgument(argc, argv, "input", inputFile))
	{
		std::cout << "provide input file via input=<filename>" << std::endl;
		return 0;
	}
	//accepts both lower and uppercase k between 1 and 15
	if (!(getCmdLineArgument(argc, argv, "k", k)|| getCmdLineArgument(argc, argv, "K", k))||k<1||k>15)
	{
		std::cout << "provide K via K=<integer between 1 and 15 inclusive>" << std::endl;
		return 0;
	}
	bool outputZeros = false;
	setToBoolIfFlagFound(argc, argv, "outputZeros", true, outputZeros);

	std::string s = get_file_contents(inputFile, true);
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;

	//this allows us to select a K at runtime even though it is a compile time constant
	switch (k)
	{
	case 1: countKmers<1>(s, outputZeros); break;
	case 2: countKmers<2>(s, outputZeros); break;
	case 3: countKmers<3>(s, outputZeros); break;
	case 4: countKmers<4>(s, outputZeros); break;
	case 5: countKmers<5>(s, outputZeros); break;
	case 6: countKmers<6>(s, outputZeros); break;
	case 7: countKmers<7>(s, outputZeros); break;
	case 8: countKmers<8>(s, outputZeros); break;
	case 9: countKmers<9>(s, outputZeros); break;
	case 10: countKmers<10>(s, outputZeros); break;
	case 11: countKmers<11>(s, outputZeros); break;
	case 12: countKmers<12>(s, outputZeros); break;
	case 13: countKmers<13>(s, outputZeros); break;
	case 14: countKmers<14>(s, outputZeros); break;
	case 15: countKmers<15>(s, outputZeros); break;
	}
    return 0;
}