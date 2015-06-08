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
			lineNumber(0),
            comment("")
		{}

		FastaContainer(const FastaContainer& other) = delete;
		FastaContainer& operator=(const FastaContainer& other) = delete;

		std::string sequence;
		unsigned long lineNumber;

        std::string comment;

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

        int getNextLine(FastaContainer& out);

        unsigned long getReadFileLength() const;
        unsigned long getReadFileRead() const;

    private:
        std::ifstream* m_readStream;
        std::ofstream* m_writeStream;

        unsigned int m_nextLineNumber;

        unsigned long m_readFileLength;
        unsigned long m_bytesRead;

        std::string m_readComment;
};

#endif
