/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:01
*/

#ifndef FASTA_IO_H
#define FASTA_IO_H

#include <fstream>
#include <string>

struct FastaContainer
{
	std::string sequence;
	unsigned long lineNumber;
	int* kMer;
	unsigned int k; // Size of k
	FastaContainer() {
		sequence = "";
		lineNumber = 0;
		kMer = NULL;
		k = -1;
	}
};

class FastaIO
{
	public:
		FastaIO();
		~FastaIO();

		FastaIO(const FastaIO&) = delete;
		FastaIO& operator=(const FastaIO&) = delete;

		void closeRead();
		void closeWrite();

		int openRead(const char* filename);
		int openWrite(const char* filename);

		bool readIsOpen() const;
		bool writeIsOpen() const;

		void writeAsync();

		int getNextLine(FastaContainer& out);

	private:
		std::ifstream* m_readStream;
		std::ofstream* m_writeStream;

		unsigned int m_nextLineNumber;
};

#endif