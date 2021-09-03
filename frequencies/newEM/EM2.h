#ifndef EM
#define EM

#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <iterator>
#include <cmath>

void EM2pop(std::map<std::string, std::string>* params, std::vector<std::string>& poplabels, std::ofstream& logfile);
double C(const int& n, int k);

#endif
