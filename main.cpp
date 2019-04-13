
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <chrono>

#include "readArgs.h"

std::string get_file_contents(const std::string &filename)
{
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	if (in)
	{
		std::string contents;
		in.seekg(0, std::ios::end);
		contents.resize(in.tellg());
		in.seekg(0, std::ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	throw(errno);
}

constexpr char chars[4] = { 'A','C','G','T' };

template <unsigned int K>
void output(int index, std::string &kmer, const std::unique_ptr<int[]> &results, double denominator, std::stringstream &ss)
{
	for (int i = 0; i < 4; ++i)
	{
		kmer[kmer.length() - K] = chars[i];
		output<K - 1>((index << 2) + i, kmer, results, denominator, ss);
	}
}

template<>
void output<0>(int index, std::string &kmer, const std::unique_ptr<int[]> &results, double denominator, std::stringstream &ss)
{
	int count = results[index] ;
	if(count!=0)
		ss << kmer << '\t' << count*denominator << '\n';
}

template <unsigned int K>
void countKmers(const std::string &s)
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
	std::stringstream ss;
	size_t i = 0;
	for (; i < s.length(); i++)
	{
		count++;
		switch (s[i])
		{
		case 'A':
		//case 'a':
			currentString += 0;
			break;
		case 'C':
		//case 'c':
			currentString += 1;
			break;
		case 'G':
		//case 'g':
			currentString += 2;
			break;
		case 'T':
		//case 't':
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
				output<K>(0, std::string(K, ' '), results, 1.0 / (double)total, ss);
				std::ofstream out("output" + std::to_string(i) + ".txt", std::ostream::trunc);
				out << ss.rdbuf();
				ss.str(std::string());
				out.close();
				total = 0;
				count = 0;
			}
			//zero the results
			memset(results.get(), 0, kmers*sizeof(int));
			size_t endOfLine = s.find_first_of('\n', i);
			ss << s.substr(i, endOfLine - i) << '\n';
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
	output<K>(0, std::string(K, ' '), results, 1.0 / (double)total, ss);
	std::ofstream out("output" + std::to_string(i) + ".txt", std::ostream::trunc);
	out << ss.rdbuf();
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
	std::string s = get_file_contents(inputFile);
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;

	//this allows us to select a K at runtime even though it is a compile time constant
	switch (k)
	{
	case 1: countKmers<1>(s); break;
	case 2: countKmers<2>(s); break;
	case 3: countKmers<3>(s); break;
	case 4: countKmers<4>(s); break;
	case 5: countKmers<5>(s); break;
	case 6: countKmers<6>(s); break;
	case 7: countKmers<7>(s); break;
	case 8: countKmers<8>(s); break;
	case 9: countKmers<9>(s); break;
	case 10: countKmers<10>(s); break;
	case 11: countKmers<11>(s); break;
	case 12: countKmers<12>(s); break;
	case 13: countKmers<13>(s); break;
	case 14: countKmers<14>(s); break;
	case 15: countKmers<15>(s); break;		
	}
    return 0;
}