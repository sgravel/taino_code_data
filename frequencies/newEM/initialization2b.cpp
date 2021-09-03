#include "initialization2.h"

using namespace std;

//Process and store commandline arguments
void processInput(int argc, char *argv[], map<string, string>*& params, vector<string>& poplabels, ofstream& logfile){
  //Step 1: Make sure an input file has been specified
  if(argc<2){
    cout<<"Error: No parameters\n";
    exit(1);
  }

  //Step 2: Get the parameters and set default values
  (*params)["-o"]="sfs";
  (*params)["-I1"]="20"; //Max number of iterations for EM 1
  (*params)["-dim1"]="100";
  (*params)["-dim2"]="100";
  (*params)["-dim3"]="20";

  //Parse commandline  
  for(int i=1; i<argc-1; i+=2){
    string s(argv[i]);
    (*params)[argv[i]]=argv[i+1];

    //In case ancestral labels are specified
    if(s=="-anc"){
      while(i<argc-1){
        if(argv[i+1][0]=='-'){
          break;
	}
        poplabels.push_back(argv[i+1]);
	i++;
      }
      i--;
    }

  }

  //Step 2a: Setup log file
  logfile.open(((*params)["-o"]+".log.txt").c_str());
  logfile<<"Jake EM begins\n1. Loading data..."<<endl;

  //Step 3: Sanity check
  checkSanity(params, poplabels, logfile);
}

void checkSanity(map<string, string>* params, vector<string>& poplabels, ofstream& logfile){
  //Check for specification of ancestral population labels
  if(!(*params).count("-anc")){
    cout<<"Error: No ancestral population labels specified. Use -anc flag to specify"<<endl;
    logfile<<"Error: No ancestral population labels specified. Use -anc flag to specify"<<endl;
    exit(1);
  }

  //Check for misspecified population labels
  if(poplabels.size()==1){
    cout<<"Error: Only 1 ancestral population specified. Need at least two"<<endl;
    logfile<<"Error: Only 1 ancestral population specified. Need at least two"<<endl;
    exit(1);
  }
  else if(poplabels.size()>3){
    cout<<"Error: More than 3 ancestral populations specified; this algorithm supports up to 3 populations"<<endl;
    logfile<<"Error: More than 3 ancestral populations specified; this algorithm supports up to 3 populations"<<endl;
    exit(1);
  }

  //Check for specification of an input file
  if(!(*params).count("-vcf")){
    cout<<"Error: No vcf file specified"<<endl;
    logfile<<"Error: No vcf file specified"<<endl;
    exit(1);
  }

  //Check existence of the input file
  ifstream input ((*params)["-vcf"].c_str());
  if (!input.is_open() ) {
    cout<<"Error opening input file "<<(*params)["-vcf"]<<endl;
    logfile<<"Error opening input file "<<(*params)["-vcf"]<<endl;
    exit(1);
  }
  input.close();

  //Check that a field ID has been specified for ancestry
  if(!(*params).count("-A")){
    cout<<"No field ID specified for ancestry; defaulting to 'LA'"<<endl;
    logfile<<"No field ID specified for ancestry; defaulting to 'LA'"<<endl;
    (*params)["-A"]="LA";
  }

  //Check for naive approach
  if((*params).count("-c")){
    (*params)["-I2"]="0";
    cout<<"Assuming a sample size cutoff of "<<(*params)["-c"]<<"; no second EM step"<<endl;
    logfile<<"Assuming a sample size cutoff of "<<(*params)["-c"]<<"; no second EM step"<<endl;
  }

  //Check for grid size
  if(atoi((*params)["-dim1"].c_str())<=0){
    cout<<"Warning: Invalid value for -dim1; setting to 20 by default"<<endl;
    logfile<<"Warning: Invalid value for -dim1; setting to 20 by default"<<endl;
    (*params)["-dim1"]="20";
  }

  if(atoi((*params)["-dim2"].c_str())<=0){
    cout<<"Warning: Invalid value for -dim2; setting to 20 by default"<<endl;
    logfile<<"Warning: Invalid value for -dim2; setting to 20 by default"<<endl;
    (*params)["-dim2"]="20";
  }

  if(atoi((*params)["-dim3"].c_str())<=0){
    cout<<"Warning: Invalid value for -dim3; setting to 20 by default"<<endl;
    logfile<<"Warning: Invalid value for -dim3; setting to 20 by default"<<endl;
    (*params)["-dim3"]="20";
  }
}
