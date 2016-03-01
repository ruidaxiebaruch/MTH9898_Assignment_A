//============================================================================
// Name        : Assignment_A.cpp
// Author      : Ruida Xie
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <cinttypes>
#include <utility>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <ctype.h>
#include <time.h>
#include <fstream>
#include <queue>
#include "mpi.h"

double MaxPrice = 0, MaxVolume = 0;
long long StartTime = 0;

// split a line into time, price, and volume
template<class T = std::string>
std::vector<T> splitString(const std::string& str, char delim) {
	size_t num = 1;
	const char* p;
	for (p = str.c_str(); *p != 0; p++) {
		if (*p == delim) {
			num++;
		}
	}
	std::vector<T> list(num);
	const char* last = str.c_str();
	num = 0;
	for (p = str.c_str(); *p != 0; p++) {
		if (*p == delim) {
			list[num++] = T(last, p - last);
			//list[num++] = boost::lexical_cast<T>(last, p - last);
			last = p + 1;
		}
	}
	list[num++] = T(last, p - last);
	//list[num++] = boost::lexical_cast<T>(last, p - last);
	return list;
}

// split a buffer with many lines into individual lines
template<class T = std::string>
std::vector<T> splitString(char* &str, char delim, int Count) {
	size_t num = 1;
	char* p;
    char* begin = str;
    if (*str == 0) {
        return std::vector<T>();
    }
	for (p = str; *p != 0 && num < Count; p++) {
		if (*p == delim) {
			num++;
		}
	}
    while (*p != 0 && *p != delim) {
        p++;
    }
    if (*p == delim) {
        *p = 0;
        str = p + 1;
    } else {
        str = p;
    }

	std::vector<T> list(num);
	char* last = begin;
	num = 0;
	for (p = begin; *p != 0; p++) {
		if (*p == delim) {
			list[num++] = T(last, p - last);
			last = p + 1;
		}
	}
	list[num++] = T(last, p - last);
	return list;
}

// obtain how many microseconds it has been from a set timestamp to the time in the data
long long getMicroSec(std::string const& s) {
    int i = 0;
    long totalMSinAYear = 365 * 24 * 3600 * pow(10, 6);
    long totalMSinAYearWithLeap = 366 * 24 * 3600 * pow(10, 6);
    long totalMSinADay = 24 * 3600 * pow(10, 6);
    // Set February with a leap day with index 0
    std::vector<int> daysInMonth = {29, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    long long dur = 0;
    int year, month, day, hour, minute, second, microSeconds;
    //std::string yearStr = "", monthStr = "", dayStr = "", hourStr = "", minuteStr = "", secondStr = "", microSecondsStr = "";
    try{
		while (i < 24) {
			if ((i == 8 || i == 11 || i == 14) && s[i] != ':') return 0;
			if (i == 17 && s[i] != '.') return 0;
			if (i != 8 && i != 11 && i != 14 && i != 17 && !isdigit(s[i])) return 0;
			i++;
		}
		year = atoi(s.substr(0, 4).c_str());
		if (year < 1970 || year > 2016) return 0;
		month = atoi(s.substr(4, 2).c_str());
		if (month < 1 || month > 12) return 0;
		day = atoi(s.substr(6, 2).c_str());
		if (day < 1 || day > 31) return 0;
		hour = atoi(s.substr(9, 2).c_str());
		if (hour >= 24) return 0;
		minute = atoi(s.substr(12, 2).c_str());
		if (minute >= 60) return 0;
		second = atoi(s.substr(15, 2).c_str());
		if (second >= 60) return 0;
		microSeconds = atoi(s.substr(18, 6).c_str());

		// how many microseconds in years before the year under discussion
		dur += (year - 1970) * totalMSinAYear;
		// adjust for leap days
		dur += (year - 1972) / 4 * totalMSinADay;
		// month
		for (int m = 1; m < month; m++){
			dur += daysInMonth[m] * totalMSinADay;
		}
		// day
		dur += (day - 1) * totalMSinADay;
		// hour
		dur += (hour - 1) * 3600 * pow(10, 6);
		// minute
		dur += (minute - 1) * 60 * pow(10, 6);
		// second
		dur += (second - 1) * pow(10, 6);
		// microseconds
		dur += microSeconds;

		return dur;
    } catch (std::exception& e) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << s << " at rank" << rank << " !!!!!!!!!!!!!\n";
        return 0;
    }
}

struct rec {
    std::string str;
    int status;
    int idate;
    long long ts;
    double fprice, fvol;
    //rec() : status(-1) {}
    rec(std::string& line) {
        status = -1;
        fprice = -1;
        fvol = -1;
        ts = 0;
        idate = 0;
	    auto cols = splitString(line, ',');
		if (cols.size() > 2) {
			ts = getMicroSec(cols[0]);
            if (ts == 0) {
                return;
            }
			idate = atoi(cols[0].substr(0,8).c_str());
            try {
			    fprice = atof(cols[1].c_str());
			    //fprice = boost::lexical_cast<double>(cols[1]);
			    //fvol = boost::lexical_cast<double>(cols[2]);
			    fvol = atof(cols[2].c_str());
            } catch (std::exception &e) {
                //std::cout << ln << std::endl;
                return;
            }
            if (fprice > 0 && fvol > 0 && fprice < 20000 && fvol < 10000000) {
            	// valid
                status = 0;
            } else {
            	// invalid
                status = 1;
            }
            str = std::move(line);
        }
    }
};


struct recBlock {
    std::vector<rec> vrec;
    double priceAvg, priceSquare, volumeAvg, volumeSquare;
    struct stat {
        float sum, sq;
        int n;
        stat() : sum(0), sq(0), n(0) {}
        void update(int x) {
            n++;
            sum += x;
            sq += pow(x, 2.0);
        }
        float sd() {
            if (n < 2) {
                return 20;
            }
            return sqrt(sq / float(n) - sum / float(pow(n, 2.0)));
        }
    };

    std::map<int, stat> dCount;
    int firstDate, secondDate;
    int fill(char* &pBuffer) {
        vrec.clear();
        if (*pBuffer == 0) {
            return 0;
        }

        std::vector<std::string> vstr = splitString(pBuffer, '\n', 10000);
        int Count = 0;
        priceAvg = priceSquare = volumeAvg = volumeSquare = 0;
        for (auto & s : vstr) {
            rec r(s);
            if (r.status < 0) {
                std::cout << "line parse error:" << s << std::endl;
                continue;
            }
            double w = double(Count) / double(Count + 1);
            double price = r.fprice;
            double volume = r.fvol;
            priceAvg = priceAvg * w + price * (1 - w);
			priceSquare = priceSquare * w + pow(price, 2.0) * (1 - w);
			volumeAvg = volumeAvg * w + volume * (1 - w);
			volumeSquare = volumeSquare * w + pow(volume, 2.0) * (1 - w);
            if (dCount[r.idate].n < 20) {
                dCount[r.idate].update(Count);
            }
            vrec.push_back(r);
			Count++;
        }
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MaxPrice = priceAvg + 6 * sqrt((-priceAvg * priceAvg + priceSquare) / float(Count));
	    MaxVolume = volumeAvg + 6 * sqrt((-volumeAvg * volumeAvg + volumeSquare) / float(Count));
        std::cout << "At rank:" << rank
        << " Count:" << Count << " priceAvg:" << priceAvg << " MaxPrice:" << MaxPrice << " volumeAvg:" << volumeAvg
        << " volumeMax:" << MaxVolume << std::endl;
        //std::cout << dCount << std::endl;
        int d1 = 0, d2 = 0, m1 = 0, m2 = 0;
        for (auto &p : dCount) {
            int d = p.first;
            stat& rr = p.second;
            int c = rr.n;
            float sd = rr.sd();
            //std::cout << "date:" << d << " Count:" << rr.n << " sd:" << sd << std::endl;
            if (c <= 20 && sd > 12) {
            	//std::cout << d<< "   sd = " << sd <<std::endl;
                continue;
            }
            if (c > m1) {
                m2 = m1; d2 = d1;
                m1 = c; d1 = d;
            } else if (c > m2) {
                m2 = c; d2 = d;
            }
        }
        firstDate = std::min(d1, d2);
        secondDate = std::max(d1, d2);
        return Count;
    }

    // return true if invalid, false if valid
    bool check_rec(rec &r) {
        if (r.status == 1) {
            return 1;
        } else {
            int curd = r.idate;
            if (curd != firstDate && curd != secondDate) {
                return 1;
            }
        }
        return 0;
    }


    void split(std::vector<std::string>& validResults, std::vector<std::string>& invalidResults) {
        long long lastTS = 0;
        int lastDate = 0;
        for (auto& record : vrec) {
            int valid = 1;
            if (check_rec(record)) {
                valid = 0;
            } else {
                if (lastTS == 0) {
                    lastTS = record.ts;
                    lastDate = record.idate;
                } else {

                    if (abs(lastTS - record.ts) > 5 * pow(10, 6)) {
                        if (lastDate == firstDate && record.idate == secondDate) {
                        	std::cout << "first: " << firstDate << std::endl;
                        	std::cout << "secon: " << secondDate << std::endl;
                        	std::cout << "last: " << lastDate << std::endl;
                        	std::cout << "idate: " << record.idate << std::endl;
                            lastTS = record.ts;
                            lastDate = record.idate;
                        } else {
                        	std::cout << "first, invalid: " << firstDate << std::endl;
                        	std::cout << "secon, invalid: " << secondDate << std::endl;
                        	std::cout << "last, invalid: " << lastDate << std::endl;
                        	std::cout << "idate, invalid: " << record.idate << std::endl;
                            valid = 0;
                        }
                    } else {
                        lastTS = record.ts;
                        lastDate = record.idate;
                    }
                }
            }
            if (valid) {
                validResults.push_back(record.str);
            } else {
                invalidResults.push_back(record.str);
            }
        }
    }

};


long long set_wBuffer(char* &p, std::vector<std::string>& v) {
    long long ret = 0;
    for (auto& str : v) {
        ret += str.length() + 1;
    }
    p = new char[ret];
    char* ptr = p;
    for (auto& str : v) {
        strcpy(ptr, str.c_str());
        ptr += str.length();
        *ptr = '\n';
        ptr++;
    }
    return ret;
}

int mpi_group_write(char* fileName, char* wBuffer, long long wSize, int rank, int nodes) {
    int ret;
    MPI_File fh;
    MPI_Status status;

    long long *allSize = (long long*)malloc(sizeof(long long) * nodes);
    MPI_Allgather(&wSize, 1, MPI_LONG, allSize, 1, MPI_LONG, MPI_COMM_WORLD);
    for (int i = 1; i < nodes; i++) {
        allSize[i] += allSize[i-1];
    }

    ret = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    ret = MPI_File_set_size(fh, 0);
    if (ret) {
        std::cout << "file open error, mpi_fopen ret:" << ret << std::endl;
        return 0;
    }
    ret = MPI_File_write_at_all(fh, allSize[rank] - wSize, (void*)wBuffer, wSize, MPI_BYTE, &status);
    if (ret) {
        std::cout << "file write error, mpi_fwrite ret:" << ret << std::endl;
        return 0;
    }
    MPI_File_close(&fh);
}

int main(int argc, char **argv){
	// get the time at the beginning of the program.
	clock_t startingTime, startReading, endReading, startWriting, endWriting;
	startingTime = clock();


    int rank, nodes;
	MPI_Offset fileSize;
    char fileName[] = "data.txt";
    MPI_File fh;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);

    startReading = clock();
	printf ("Beginning to read. It's been %f seconds", ((float)startReading - startingTime)/CLOCKS_PER_SEC);
	std::cout << "\n";

    int ret = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (ret) {
        std::cout << "file open error, mpi_fopen ret:" << ret << std::endl;
        return 0;
    }
    ret = MPI_File_get_size(fh, &fileSize);
    std::cout << "mpi_fileSize ret:" << ret << " file size:" << fileSize << std::endl;
    long long tmpSize = fileSize / nodes;
    long long readCount = rank==nodes-1 ? fileSize - (nodes-1)*tmpSize : tmpSize;
    std::cout << "rank is:" << rank << " nodes number is:" << nodes << " read from " << rank * tmpSize
        << " to " << rank * tmpSize + readCount << std::endl;

    char* Buffer = new char[readCount + 1];



    ret = MPI_File_read_at(fh, rank*tmpSize, Buffer, readCount, MPI_BYTE, &status);
    MPI_File_close(&fh);
    if (ret) {
        std::cout << "mpi read error code:" << ret;
        delete[] Buffer;
        return 0;
    }

    endReading = clock();
	printf ("End reading. It's been %f seconds", ((float)endReading - startingTime)/CLOCKS_PER_SEC);
	std::cout << "\n";

    Buffer[readCount] = 0;
    Buffer[readCount-1] = 0;
    recBlock block;
    char* pBuffer = Buffer;
    int Count = block.fill(pBuffer);
    std::vector<std::string> validResults;
    std::vector<std::string> invalidResults;
    int countSum = Count;
    while (Count) {
        block.split(validResults, invalidResults);
        char* orgp = pBuffer;
        Count = block.fill(pBuffer);
        countSum += Count;
    }
    delete[] Buffer;
    //std::cout << "rank:" << rank << " sum:" << countSum << " invalid:" << invalidResults.size();
    char* wBuffer;
    long long wSize = set_wBuffer(wBuffer, invalidResults);
    char fileName2[] = "noise_data.txt";

    startWriting = clock();
	printf ("Beginning to write. It's been %f seconds", ((float)startWriting - startingTime)/CLOCKS_PER_SEC);
	std::cout << "\n";

    ret = mpi_group_write(fileName2, wBuffer, wSize, rank, nodes);
    if (wBuffer) {
        delete[] wBuffer;
    }
    invalidResults.clear();
    wSize = set_wBuffer(wBuffer, validResults);
    char fileName3[] = "valid_data.txt";
    ret = mpi_group_write(fileName3, wBuffer, wSize, rank, nodes);
    if (wBuffer) {
        delete[] wBuffer;
    }

    endWriting = clock();
	printf ("End writing. It's been %f seconds", ((float)endWriting - startingTime)/CLOCKS_PER_SEC);
	std::cout << "\n";


    MPI_Finalize();
    return 0;
}
