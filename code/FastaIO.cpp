/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:12
* @Last Modified time: 2015-03-03 01:02:12
* @Version: 0.0
*/

#include "FastaIO.h"

FastaIO::FastaIO()
{
	m_readStream = NULL;
}

FastaIO::FastaIO(const FastaIO& fasta)
{
}

FastaIO::~FastaIO()
{
	closeRead();
	//closeWrite();
}

void FastaIO::closeRead()
{
	if(m_readStream != NULL)
	{
		if(m_readStream->is_open())
			m_readStream->close();

		delete m_readStream;
		m_readStream = NULL;
	}
}

int FastaIO::openRead(const char* filename)
{
	closeRead();
	m_readStream = new std::ifstream(filename);

	if(!m_readStream->is_open()) {
		closeRead();
		return -1;
	}

	return 0;
}

int FastaIO::openWrite(const char* filename)
{
	return 0;
}

int FastaIO::getNextLine(std::string& out)
{
	while(std::getline(*m_readStream, out)) {
		if(out[0] != '>')
			return 0;
	}

	return -1;
}
