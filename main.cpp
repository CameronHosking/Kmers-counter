
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <chrono>

std::string get_file_contents(const char *filename)
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
	std::string s = get_file_contents("human_g1k_v37.fasta");
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;
	countKmers<4>(s);
	int itl;
	std::cin >> itl;
    return 0;
}