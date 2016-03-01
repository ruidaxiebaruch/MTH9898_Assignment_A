//============================================================================
// Name        : testNormal.cpp
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
using namespace std;

long long getMicroSec(std::string const& s) {
    int i = 0;
    long totalMSinAYear = 365 * 24 * 3600 * pow(10, 6);
    //long totalMSinAYearWithLeap = 366 * 24 * 3600 * pow(10, 6);
    long totalMSinADay = 24 * 3600 * pow(10, 6);
    // Set February with a leap day with index 0
    std::vector<int> daysInMonth = {29, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    long long dur = 0;
    int year, month, day, hour, minute, second, microSeconds;
    //std::string yearStr = "", monthStr = "", dayStr = "", hourStr = "", minuteStr = "", secondStr = "", microSecondsStr = "";

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

}



int main() {



    // test of normality is conducted separately as follows
    // pick a sample of up to 100,000 data points

    std::string line;
    std::ifstream validDataFile("valid_data.txt");
    std::queue<std::string> timeQ, priceQ;
    int countTotal = 0, countSelect = 0;
	std::string time = "", price = "";
	// compute the annualized 1st, 2nd, 3rd, and 4th moments of log-returns
	double firstM, secondM, thirdM, fourthM;
	firstM = secondM = thirdM = fourthM = 0;
    if (validDataFile.is_open()){
    	getline(validDataFile, line);
    	countTotal++;  // 1
    	countSelect++; // 1
    	int j = 0;
    	while (line[j] != ','){
			time += line[j];
			j++;
		}
		timeQ.push(time);
		time = "";
		j++;
		while (line[j] != ','){
			price += line[j];
			j++;
		}
		priceQ.push(price);
		price = "";
    	while(getline(validDataFile, line)){
			countTotal++;
			// sample one data point for every 1000 data
    		if ((countTotal) % 1000 == 0){
    			countSelect++;
				j = 0;
				// obtain the time data as string
				while (line[j] != ','){
					time += line[j];
					j++;
				}
				timeQ.push(time);
				//cout << time << endl;
				time = "";
				j++;
				// obtain the price data as string
				while (line[j] != ','){
					price += line[j];
					j++;
				}
				priceQ.push(price);
				//cout << price << endl;
				price = "";



	    		firstM = firstM * (countSelect-1) / (countSelect) + (log(atof(priceQ.back().c_str()) / atof(priceQ.front().c_str()))) / (getMicroSec(timeQ.back()) - getMicroSec(timeQ.front())) * pow(10, 6) * 3600 * 24 * 365.25 / (countSelect);
	    		secondM = secondM * (countSelect-1) / (countSelect) +  pow((log(atof(priceQ.back().c_str()) / atof(priceQ.front().c_str()))) / (getMicroSec(timeQ.back()) - getMicroSec(timeQ.front())) * pow(10, 6) * 3600 * 24 * 365.25, 2.0) / (countSelect);
				thirdM = thirdM * (countSelect-1) / (countSelect) + pow((log(atof(priceQ.back().c_str()) / atof(priceQ.front().c_str()))) / (getMicroSec(timeQ.back()) - getMicroSec(timeQ.front())) * pow(10, 6) * 3600 * 24 * 365.25, 3.0) / (countSelect);
				fourthM = fourthM * (countSelect-1) / (countSelect) + pow((log(atof(priceQ.back().c_str()) / atof(priceQ.front().c_str()))) / (getMicroSec(timeQ.back()) - getMicroSec(timeQ.front())) * pow(10, 6) * 3600 * 24 * 365.25, 4.0) / (countSelect);

				// not to store too much data that are no longer useful
			    priceQ.pop();
			    timeQ.pop();
    		}


    	}
    }else{
    	std::cout << "Unable to open valid_data.txt" << std::endl;
    }
    validDataFile.close();

    // Jarque-Bera Test
    std::cout << "The first moment is " << firstM << std::endl;
    std::cout << "The second moment is " << secondM << std::endl;
    std::cout << "The third moment is " << thirdM << std::endl;
    std::cout << "The fourth moment is " << fourthM << std::endl;

    double skewness = (thirdM - 3 * firstM * secondM + 3 * pow(firstM, 2.0) * firstM - pow(firstM, 3.0))/(pow((secondM - pow(firstM, 2.0)), 1.5));
	double kurtosis = (fourthM - 8 * firstM * thirdM + 6 * secondM * secondM + pow(firstM, 4.0)) / (pow((secondM - pow(firstM, 2.0)), 2.0));


    std::cout << "The skewness is " << skewness << std::endl;
    std::cout << "The kurtosis is " << kurtosis << std::endl;

    double JB = (countSelect / 6.0) * (skewness*skewness + 0.25 * pow((kurtosis - 3), 2.0));


    std::cout << "The test statistic is " << JB << std::endl;

    if (JB > 0.002) {
    	std::cout << "The returns are not normally distributed." << std::endl;
    }else{
    	std::cout << "The returns are normally distributed." << std::endl;
    }

	return 0;
}
