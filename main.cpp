
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
void output(size_t i, std::string id, bool split, bool outputZeros, const uint32_t results[], size_t total)
{
	//output string
	std::string outFileName = "output" + std::to_string(i) + ".txt";
	std::ofstream out(outFileName, std::ostream::trunc | std::ostream::binary);
	std::cout << "outputting Kmers to: " << outFileName << std::endl;
	if (split) {
		out << id << std::endl;
	}
	if (outputZeros)
	{
		output<K, true>(0, std::string(K, ' '), results, 1.0 / (double)total, out);
	}
	else
	{
		output<K, false>(0, std::string(K, ' '), results, 1.0 / (double)total, out);
	}
	out.close();
}

//call at the start of an id line, splits and outputs is split is true
//set the id to the contents of this line
//return the index at the end of the line
template<unsigned int K>
size_t newID(bool split, size_t total, size_t i, std::string &id, bool outputZeros, uint32_t results[], const std::string &s)
{
	uint64_t kmers = 1LL << K * 2;
	//if split is true then we ouput individual counts for each id found
	if (split && total > 0)
	{
		output<K>(i, id, split, outputZeros, results, total);
		//zero the results
		memset(results, 0, kmers*sizeof(int));
	}
	size_t endOfLine = s.find_first_of('\n', i);
	id = s.substr(i, endOfLine - i);
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
void countKmers(const std::string &s, bool outputZeros, bool split)
{
	auto start = std::chrono::high_resolution_clock::now();
	//allocate memory for the 4^K different kmers
	uint64_t kmers = 1LL << K * 2;
	std::unique_ptr<uint32_t[]> results = std::unique_ptr<uint32_t[]>(new uint32_t[kmers]());
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
		else if(c== '>')
		{
			i = newID<K>(split, total, i, id, outputZeros, results.get(), s);
			if (split)
			{
				total = 0;
			}
			count = 0;
		}
		//other characters are ignored eg \n\r etc
	}
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
	start = std::chrono::high_resolution_clock::now();
	
	output<K>(i, id, split, outputZeros, results.get(), total);

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

	bool split = false;
	setToBoolIfFlagFound(argc, argv, "split", true, split);

	std::string s = get_file_contents(inputFile, true);
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;

	//this allows us to select a K at runtime even though it is a compile time constant
	switch (k)
	{
	case 1: countKmers<1>(s, outputZeros, split); break;
	case 2: countKmers<2>(s, outputZeros, split); break;
	case 3: countKmers<3>(s, outputZeros, split); break;
	case 4: countKmers<4>(s, outputZeros, split); break;
	case 5: countKmers<5>(s, outputZeros, split); break;
	case 6: countKmers<6>(s, outputZeros, split); break;
	case 7: countKmers<7>(s, outputZeros, split); break;
	case 8: countKmers<8>(s, outputZeros, split); break;
	case 9: countKmers<9>(s, outputZeros, split); break;
	case 10: countKmers<10>(s, outputZeros, split); break;
	case 11: countKmers<11>(s, outputZeros, split); break;
	case 12: countKmers<12>(s, outputZeros, split); break;
	case 13: countKmers<13>(s, outputZeros, split); break;
	case 14: countKmers<14>(s, outputZeros, split); break;
	case 15: countKmers<15>(s, outputZeros, split); break;
	}
    return 0;
}