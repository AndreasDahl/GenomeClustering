/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 01:02:01
*/

#ifndef FASTA_IO_H
#define FASTA_IO_H

#include <fstream>
#include <string>

class FastaContainer
{
    public:
        FastaContainer() :
            sequence(""),
            lineNumber(0),
            kMerLength(0),
            k(-1),
            m_kMer(NULL)
        {}

        FastaContainer(const FastaContainer& other) :
            sequence(other.sequence),
            lineNumber(other.lineNumber),
            kMerLength(other.kMerLength),
            k(other.k),
            m_kMer(NULL)
        {
            if(kMerLength > 0) {
                m_kMer = new int[kMerLength];
                for(unsigned int i = 0; i < kMerLength; i++)
                    m_kMer[i] = other.m_kMer[i];
            }
        }

        ~FastaContainer() {
            if(m_kMer) delete m_kMer;
        }

        FastaContainer& operator=(const FastaContainer& other) {
            if(&other != this) {
                sequence = other.sequence;
                lineNumber = other.lineNumber;
                kMerLength = other.kMerLength;
                k = other.k;
                if(m_kMer) delete m_kMer;
                m_kMer = (kMerLength > 0) ? new int[kMerLength] : NULL;
            }
            return *this;
        }

        void setKMerLength(unsigned int KMerLength) {
            kMerLength = KMerLength;
            if(m_kMer) delete m_kMer;
            m_kMer = (kMerLength > 0) ? new int[kMerLength] : NULL;
        }

        int* getKMer() {
            return m_kMer;
        }

        const int* getKMer() const {
            return m_kMer;
        }

        std::string sequence;
        unsigned long lineNumber;

        unsigned int kMerLength;
        unsigned int k; // Size of k

    private:
        int* m_kMer;
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
