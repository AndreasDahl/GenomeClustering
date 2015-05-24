/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:12
*/

#include "FastaIO.h"

FastaIO::FastaIO() :
    m_readStream(NULL),
    m_writeStream(NULL),
    m_nextLineNumber(0),
    m_readFileLength(0),
    m_bytesRead(0)
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

    m_nextLineNumber = 1;

    if(!m_readStream->is_open()) {
        closeRead();
        return -1;
    }

    m_readStream->seekg(0, m_readStream->end);
    m_readFileLength = (long)m_readStream->tellg();
    m_readStream->seekg(0, m_readStream->beg);

    m_bytesRead = 0;

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

bool FastaIO::readIsOpen() const
{
    return m_readStream != NULL;
}

bool FastaIO::writeIsOpen() const
{
    return m_writeStream != NULL;
}

void FastaIO::writeAsync()
{

}

int FastaIO::getNextLine(FastaContainer& out)
{
    std::string str;

    out.sequence = "";
    out.lineNumber = 0;

    while(std::getline(*m_readStream, str))
    {
        m_bytesRead += (long)str.size();
        if(str[0] != '>') {
            out.sequence += str;
            if(out.lineNumber == 0) {
                out.lineNumber = m_nextLineNumber;
            }
            m_nextLineNumber += 1;
        }
        else if(out.sequence != "") {
            return 0;
        }
    }

    return -1;
}

unsigned long FastaIO::getReadFileLength() const
{
    return m_readFileLength;
}

unsigned long FastaIO::getReadFileRead() const
{
    return m_bytesRead;
}
