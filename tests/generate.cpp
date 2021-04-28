#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#include <boost/filesystem.hpp>

int main(int argc, char **argv)
{
  if (2 != argc)
  {
    std::cerr << "Usage: " << argv[0] << " path/to/fasta/reference/genome.fa" << std::endl;
    std::cerr 
      << "\ngenerate a single fastq from the first reference sequence found in the"
      << "\nfasta file. Each read uses 2 lines in the fasta file and skips one line"
      << std::endl;
    exit (1);
  }
  const boost::filesystem::path fasta = argv[1];
  std::cerr << "\ngenerating fastq from fasta input: " << fasta << std::endl;
  std::ifstream is(fasta.c_str());
  std::string line, s1, s2;
  // skip 
  while (std::getline(is, line) && line.empty())
  {
  }
  if (line.empty() || ('>' != line[0]))
  {
    std::cerr << "failed to find reference sequence: first line is: " << line << std::endl;
    exit (2);
  }
  std::cerr << "using reference sequence: " << line << std::endl;
  const auto ws = line.find(' ');
  std::cerr << "ws: " << ws << std::endl;
  const std::string name = line.substr(1, (ws == std::string::npos ? ws : ws -1));
  std::cerr << "name: " << name << std::endl;
  std::ostringstream read;
  const char Q = '@' + 32;
  std::string qualities;
  size_t position = 0;
  while(std::getline(is, line) && (!line.empty()) && ('>' != line[0]))
  {
    position += line.length();
    if (getline(is, s1) && (!s1.empty()) && ('>' != s1[0]) &&
        getline(is, s2) && (!s2.empty()) && ('>' != s2[0]))
    {
      qualities.resize(s1.length() + s2.length(), Q);
      read.str("");
      read << name << ':' << std::setw(10) << std::setfill('0') <<  position << std::setfill(' ') << '\n' << s1 << s2 << "\n+\n" << qualities;
      position += s1.length() + s2.length();
      if ((std::string::npos != s1.find('N')) || (std::string::npos != s2.find('N')))
      {
        continue;
      }
      std::cout << read.str() << std::endl;
      //break;
    }
    else
    {
      break;
    }
  }
}
