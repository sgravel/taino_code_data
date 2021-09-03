#ifndef INITIALIZATION
#define INITIALIZATION

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

//Processing functions
void processInput(int argc, char *argv[], std::map<std::string, std::string>*& params, std::vector<std::string>& poplabels, std::ofstream& logfile);
void checkSanity(std::map<std::string, std::string>* params, std::vector<std::string>& poplabels, std::ofstream& logfile);

#endif
