#include <iostream>
#include <fstream>
#include <string>
#include "picosha2.h"

using namespace std;


int main() {
	ofstream fout("database.txt");
	int sf = 1000000000;
	for(int i=0; i<sf; i++) {//64*sf characters
		stringstream strs;
		strs << i;
		string hashed_str;
		picosha2::hash256_hex_string(strs.str(), hashed_str);
		fout << hashed_str;
	}
	return 0;
}
