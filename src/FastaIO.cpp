/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:12
*/

#include "FastaIO.h"

FastaIO::FastaIO() :
	m_readStream(NULL),
	m_writeStream(NULL)
{}

FastaIO::~FastaIO()
{
	closeRead();
	closeWrite();
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

void FastaIO::closeWrite()
{
	if(m_writeStream != NULL)
	{
		if(m_writeStream->is_open())
			m_writeStream->close();

		delete m_writeStream;
		m_writeStream = NULL;
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
	closeWrite();
	m_writeStream = new std::ofstream(filename);

	if(!m_writeStream->is_open()) {
		closeRead();
		return -1;
	}

	return 0;
}

int FastaIO::getNextLine(std::string& out)
{
	out = "";

	while(std::getline(*m_readStream, m_nextString)) {
		if(m_nextString[0] != '>') {
			out += m_nextString;
		}
		else if(out != "") {
			return 0;
		}
	}

	return -1;
}
