#include "initialization2.h"
#include "EM2.h"

using namespace std;

//This version with EM algorithm
int main(int argc, char *argv[]){
  //0. Set up the log file
  ofstream logfile;
  cout<<"Jake EM begins"<<endl;
  cout<<"Loading data"<<endl;

  //1. Collect the required information
  map<string, string>* params=new map<string, string>();
  vector<string> poplabels;
  processInput(argc, argv, params, poplabels, logfile);

  //2. Invoke EM algorithm 1
  if(poplabels.size()==2){
    EM2pop(params, poplabels, logfile);
  }

  //3. Cleanup
  delete params;
}
