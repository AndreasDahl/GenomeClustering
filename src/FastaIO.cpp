/*
 * Copyright 2015 Andreas Dahl, Christian Muf
 *
 * This file is part of MufDahlClust.
 *
 * MufDahlClust is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MufDahlClust is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MufDahlCLust.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * @Author: Christian Muf
 * @Author: Andreas Dahl
 * @Date:   2015-03-03 01:02:12
 */

#include "FastaIO.h"

FastaIO::FastaIO() :
    m_readStream(NULL),
    m_writeStream(NULL),
    m_nextLineNumber(0),
    m_readFileLength(0),
    m_bytesRead(0),
    m_readComment("")
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

int FastaIO::getNextLine(FastaContainer& out)
{
    if(m_readComment.size() > 0) {
        out.comment = std::string(m_readComment, 1);
    } else {
        out.comment = "";
    }

    out.sequence = "";
    out.lineNumber = 0;

    while(std::getline(*m_readStream, m_readComment))
    {
        m_bytesRead += (long)m_readComment.size();
        if(m_readComment[0] != '>') {
            out.sequence += m_readComment;
            if(out.lineNumber == 0) {
                out.lineNumber = m_nextLineNumber;
            }
            m_nextLineNumber += 1;
        }
        else if(out.sequence == "") {
            out.comment = std::string(m_readComment, 1);
        }
        else {
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
