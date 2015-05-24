/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:01
*/

#ifndef FASTA_IO_H
#define FASTA_IO_H

#include <fstream>
#include <string>

#include "KMerHash.h"

class FastaContainer
{
	public:
		FastaContainer() :
			sequence(""),
			lineNumber(0)
		{}

		FastaContainer(const FastaContainer& other) :
			sequence(other.sequence),
			lineNumber(other.lineNumber)
		{}

		FastaContainer& operator=(const FastaContainer& other) {
			if(&other != this) {
				sequence = other.sequence;
				lineNumber = other.lineNumber;
			}
			return *this;
		}

		std::string sequence;
		unsigned long lineNumber;

		KMerHashmap kMerHash;
};

class FastaIO
{
	public:
		FastaIO();
		virtual ~FastaIO();

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

        unsigned long getReadFileLength() const;
        unsigned long getReadFileRead() const;

    private:
        std::ifstream* m_readStream;
        std::ofstream* m_writeStream;

        unsigned int m_nextLineNumber;

        unsigned long m_readFileLength;
        unsigned long m_bytesRead;
};

#endif
