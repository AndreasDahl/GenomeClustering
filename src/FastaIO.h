/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:01
*/

#ifndef FASTA_IO_H
#define FASTA_IO_H

#include <fstream>
#include <string>

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

		int getNextLine(std::string& out);

	private:
		std::ifstream* m_readStream;
		std::ofstream* m_writeStream;

		std::string m_nextString;
};

#endif