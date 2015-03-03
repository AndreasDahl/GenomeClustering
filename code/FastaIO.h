/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:01
* @Last Modified time: 2015-03-03 01:02:01
* @Version: 0.0
*/

#ifndef FASTA_IO_H
#define FASTA_IO_H

#include <fstream>
#include <string>

class FastaIO
{
	public:
		FastaIO();
		FastaIO(const FastaIO& fasta);
		~FastaIO();

		void closeRead();
		void closeWrite();

		int openRead(const char* filename);
		int openWrite(const char* filename);

		int getNextLine(std::string& out);

	private:
		std::ifstream* m_readStream;
};

#endif