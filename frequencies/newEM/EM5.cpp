#include "EM2.h"

using namespace std;

//Version 5: Output ML frequencies and confidence intervals per site
//EM round 1: Inferring observed counts of allele 1 in each ancestral population
void EM2pop(map<string, string>* params, vector<string>& poplabels, ofstream& logfile){
  cout<<"2 population EM begins"<<endl;
  logfile<<"2 population EM begins"<<endl;

  //1. Pre-processing of parameters, set up temporary constructs
  map<string, int> ancpops;
  vector<int>* genotypes=new vector<int>();
  vector<int>* anc=new vector<int>[2];
  const int iter=atoi((*params)["-I1"].c_str());

  //2. Define a map of file genotypes to actual genotypes
  //Code: 0=ref/ref, 1=ref/alt, 2=alt/alt, 3=msg
  map<string, int> genotyper;
  genotyper["0/0"]=genotyper["0|0"]=0;
  genotyper["0/1"]=genotyper["0|1"]=genotyper["1/0"]=genotyper["1|0"]=1;
  genotyper["1/1"]=genotyper["1|1"]=2;

  //3. Map ancestral populations to their index in the poplabels vector
  for(int i=0; i<poplabels.size(); i++){
    ancpops[poplabels[i]]=i;
  }
  
	
  //4. Initiate the joint SFS
  const int dim1=atoi((*params)["-dim1"].c_str()); 
  const int dim2=atoi((*params)["-dim2"].c_str()); 
  double** sfs=new double*[dim1+1];
  double** sfs_temp=new double*[dim1+1];
  double** sfs_next=new double*[dim1+1];

  for(int i=0; i<=dim1; i++){
    sfs[i]=new double[dim2+1];
    sfs_temp[i]=new double[dim2+1];
    sfs_next[i]=new double[dim2+1];
    for(int j=0; j<=dim2; j++){
      //      sfs[i][j]=1/((double)((i+4)*(j+4)));
      sfs[i][j]=1; //Uniform prior
      sfs_next[i][j]=0;
    }
  }

  //5. Iteratively scan through the entire file
  ofstream frequencies(((*params)["-o"]+".freqs.txt").c_str());
  for(int it=0; it<iter; it++){
    cout<<"EM Iteration: "<<it<<endl;
   
    //Open file
    ifstream input((*params)["-vcf"].c_str());
    int M=0;

    //Skip the lines containing meta-info (these start with #)
    string s;
    while(!input.eof()){
      getline(input,s);
      if(s.find("##")==string::npos){
	break;
      }
    }
 
    //Now go over the file
    while(!input.eof()){
      getline(input, s);

      //Breakdown by space
      stringstream iss(s);
      vector<string> tokens;
      copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter< vector<string> >(tokens));

      if(tokens.size()==0){
	break;
      }       

      //If genotype and local ancestry information are not both present, skip this marker
      if(tokens[8].find("GT")==string::npos || tokens[8].find((*params)["-A"])==string::npos){
	getline(input, s);
	continue;
      }
      
      //Clear all conctructs from previous marker
      (*genotypes).clear();
      anc[0].clear();
      anc[1].clear();
      
      //Now breakdown the FORMAT part by colon (:)
      //1. Find index of the GT and local ancestry fields
      int indGT=(int)(((double)(tokens[8].find("GT")))/3);
      int indLA=(int)(((double)(tokens[8].find((*params)["-A"])))/3);

      //2. Iterate over the tokens, store genotype and local ancestry info from the file once
      bool mono=true;
		int numfailed=0;
		for(int j=9; j<tokens.size(); j++){
	int ctr=0;
	string item;
	string geno;
	string la;
	stringstream iss2(tokens[j]);
	bool failedgen=false;
	//Genotype and local ancestry fields
	while(ctr<=indGT || ctr<=indLA) {
	  getline(iss2, item, ':');
	  if(ctr==indGT){
	    geno=item;
	  }
	  else if(ctr==indLA){
	    la=item;
	  }
	  ctr++;
	}
	
	//Parse in the genotype
	if(genotyper.count(geno)){
	  (*genotypes).push_back(genotyper[geno]);
	  if(genotyper[geno]!=0){
	    mono=false;
	  }
	}
	else{
	  
	  failedgen=true;
		numfailed+=1;
		continue;//proceed to next token	
	}
	
	//Parse in the local ancestry field
	stringstream iss3(la);
	for(int i=0; i<2; i++){
	  getline(iss3, item, ',');
	  if(!ancpops.count(item)){
	    //cout<<"Error: Unrecognized local ancestry field "<<item<<endl;
	    //logfile<<"Error: Unrecognized local ancestry field "<<item<<endl;
	    failedgen=true;
		  if(i==0){
			  numfailed+=1;
			  break; //exit the loop over ancestries. The first ancestry field should trip this. If first field is ok but second isn't, stop 
		  }
		  else {
			  cout<<"Error: Unrecognized local ancestry field "<<item<<endl;
			  logfile<<"Error: Unrecognized second local ancestry field "<<item<<endl;
			  exit(0);
		  }

		  }	
	  else{
		  anc[i].push_back(ancpops[item]);}
		
	}
			
			
			
		}

      //Check for monomorphic
      if(mono){
	continue; //in principle we could keep the monomorphic sites
      }
	  //Check if too many genotypes or ancestries are failed
		if(numfailed>10){
			cout<<"warning: more than 10 failures "<<input<<endl;
			logfile<<"warning: Unrecognized genotype "<<input<<endl;
			continue;
		}

      //Now integrate over all pairs of frequencies
      double denom=0;
      for(int i1=0; i1<=dim1; i1++){
	    for(int i2=0; i2<=dim2; i2++){
	  sfs_temp[i1][i2]=1;
	  double f[]={(double)(i1)/((double)dim1), (double)(i2)/((double)dim2)};//the estimated frequency. Changed so that it could be between 0 and 1. 
	  
	//calculate g=P(D_i|f,R)		
	for(int j=0; j<(*genotypes).size(); j++){
	    double g=0;
	    switch((*genotypes)[j]){
	    case 0:
	      g=(1-f[anc[0][j]])*(1-f[anc[1][j]]); 
	      break;
	    case 1:
	      g=f[anc[0][j]]*(1-f[anc[1][j]])+f[anc[1][j]]*(1-f[anc[0][j]]);
	      break;
	    case 2:
	      g=f[anc[0][j]]*f[anc[1][j]];
	      break;
	    }

	    //Update the MLE
	    sfs_temp[i1][i2]*=g; //the overall likelihood is the product of the genotype likelihoods: This is P(D|f,R)
	  }
	  sfs_temp[i1][i2]*=sfs[i1][i2]; //this is the prior SFS, to get the numerator in the expression for P(f|DR) and the integrand in the denominator
	  denom+=sfs_temp[i1][i2]; //the denominator. sfs_temp/denom is P(f|DR)
	}
      }      

      //Next update SFS posterior at this site
      if(denom>0){
	for(int i1=0; i1<=dim1; i1++){
	  for(int i2=0; i2<=dim2; i2++){
	    sfs_next[i1][i2]+=sfs_temp[i1][i2]/denom;
	  }
	}   
	M++;
      }

      //Iterate over ancestral populations
      if(it==iter-1 && denom>0){
	//Compute ML frequencies
	bool lower[]={false,false};
	bool upper[]={false,false};
	double lower_bound[]={0,0};
	double upper_bound[]={0,0};
	double mean[]={0,0};

	for(int a=0; a<2; a++){
	  double cumsum=0;
	  if(a==0){
	    for(int i1=0; i1<=dim1; i1++){
	      double sfs_marginal=0;
	      for(int i2=0; i2<=dim2; i2++){
		sfs_marginal+=sfs_temp[i1][i2]/denom;
	      }
	      cumsum+=sfs_marginal;
	      
	      //Check lower bound
	      if((cumsum>0.025) && (lower[a]==false)){
		lower[a]=true;
		lower_bound[a]=((double)i1)/((double) dim1);
	      }
	      
	      //Check upper bound
	      if((cumsum>0.975) && (upper[a]==false)){
		upper[a]=true;
		upper_bound[a]=((double)i1)/((double) dim1);
	      }
	      
	      //Augment the mean	    
	      mean[a]+=((double)i1)/((double) dim1)*sfs_marginal;
	    }
      
	    //Degenerate cases
	    if(upper[a]==false){
	      upper_bound[a]=1;
	    }
	  }
	  else{
	    for(int i2=0; i2<=dim2; i2++){
              double sfs_marginal=0;
              for(int i1=0; i1<=dim1; i1++){
                sfs_marginal+=sfs_temp[i1][i2]/denom;
              }
	      cumsum+=sfs_marginal;	      

              //Check lower bound
              if((cumsum>0.025) && (lower[a]==false)){
                lower[a]=true;
                lower_bound[a]=((double)i2)/((double) dim2);
              }
	      
              //Check upper bound
	      if((cumsum>0.975) && (upper[a]==false)){
                upper[a]=true;
                upper_bound[a]=((double)i2)/((double) dim2);
              }
	      
              //Augment the mean 
              mean[a]+=((double)i2)/((double) dim2)*sfs_marginal;
            }

            //Degenerate cases
            if(upper[a]==false){
              upper_bound[a]=1;
            }
	  }
	}

	//Print
	frequencies<<tokens[0]<<" "<<tokens[1]<<" "<<lower_bound[0]<<" "<<mean[0]<<" "<<upper_bound[0]<<" "<<lower_bound[1]<<" "<<mean[1]<<" "<<upper_bound[1]<<endl;
      }
    }

    //Update full SFS, reset temporary variables
    for(int i1=0; i1<dim1; i1++){
      for(int i2=0; i2<dim2; i2++){
	sfs[i1][i2]=sfs_next[i1][i2]/((double) M); //note that this is not the "real" SFS, but the sfs among non-monomorphic sites. 
	sfs_next[i1][i2]=0;
      }
    }
  }

  //6. SFS output
  //Population 0
  ofstream outfile(((*params)["-o"]+"."+poplabels[0]+".txt").c_str());
  for(int i1=0; i1<dim1; i1++){
    double sfs_marginal=0;
    for(int i2=0; i2<dim2; i2++){
      sfs_marginal+=sfs[i1][i2];
    }
    outfile<<((double)i1)/((double)dim1-1)<<" "<<sfs_marginal<<endl;
  }

  //Population 1
  ofstream outfile2(((*params)["-o"]+"."+poplabels[1]+".txt").c_str());
  for(int i2=0; i2<dim2; i2++){
    double sfs_marginal=0;
    for(int i1=0; i1<dim1; i1++){
      sfs_marginal+=sfs[i1][i2];
    }
    outfile2<<((double)i2)/((double)dim2-1)<<" "<<sfs_marginal<<endl;
  }

  //Cleanup
  for(int i1=0; i1<dim1; i1++){
    delete[] sfs[i1];
    delete[] sfs_temp[i1];
    delete[] sfs_next[i1];
  }

  delete[] sfs;
  delete[] sfs_next;
  delete[] sfs_temp;
  delete genotypes;
  delete[] anc;
}
