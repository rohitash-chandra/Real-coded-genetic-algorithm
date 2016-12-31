/*
*	Created by: Dr. Rohitash Chandra (2006)(c.rohitash@gmail.com)
*	Software Foundation Fiji (Artficial Intelligence and Cybernetics Research Group--www.softwarefoundationfiji.org/aicrg) 

 complies in g++ under any linux environment. Can also compile in Windows. 

includes 1. Wrights Heuristic Crossoever, 2. Non Uniform Mutation 3. Roulette Wheel Selection 4 benchmark functions (Sphere and Rosen)
*/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <sys/time.h>  
#include <stdio.h>
#include <time.h>

time_t TicTime;
time_t TocTime;

using namespace::std;
//type delcalaration
typedef vector<double> Nodes;
typedef vector<double> Layer;
typedef vector<int> Sizes;
typedef vector<vector<double> > Weight;
typedef vector<vector<double> > Data;

const int MaxFE = 1000000;
  
class GeneticAlgorithm{ 
       protected:
            Data Chromozone;
			Data FinalChromozone;
			Data FirstChromozone;
			Data NewGeneration;
			Layer Fitness;
			Layer FitnessRatio;   
			Data Distance;
			
			int Population;
			int Stringsize;
			int chrome;       	 
			int Time;
			int function; 
			
			    int TotalEval;

                double Error;

                double Cycles;
                bool Success;
       //---------------------------
               public:
			GeneticAlgorithm(int population,int stringSize, int func)
			{ 
				chrome = 0;
				Population = population;//population size
				Stringsize = stringSize;//chromozone size
				function=func;//function to evaluate
			}    
      
			double Random() ; 
			
			int GetEval(){
              		  return TotalEval;
                      }
            		  double GetCycle(){
                         return Cycles;
                               }
              		double GetError(){
                        return Error;
                              }

             		 bool GetSucess(){
                                  return Success;
                                        } 
		 
       
			double FitnessFunc(Layer x, int ProbNum);
			
			int  Algorithmn(double Crossover,double Mutation, ofstream &output1, ofstream &output2);   
			
			double Evaluate()  ;
   
			int Select() ;
   
			void CrossoverAndMutate(int leftpair,int rightpair,double Crossover,double Mutation,int position);
      
			void InitilisePopulation(int population, int stringsize );
    
			double MaxFitness(); 
    
			void Print()  ;
   
			int MaxLocation();
  
			double RandomWeights();

			double RandomAddition();
    		 
};
   

double GeneticAlgorithm::Random()
{     
    double string;
    return string = rand()%2;    
     
}

double GeneticAlgorithm::RandomWeights()
{     
    int chance;
    double randomWeight;
    double NegativeWeight;
    chance =rand()%2;//randomise between negative and positive numbers
      
    if(chance ==0){
		return (((rand()%10000)*0.0001) * 5);
    }
     
    if(chance ==1){
		return ((rand()%10000)*0.0001) * -5;
	}
     
}

double GeneticAlgorithm::RandomAddition()
{     
    int chance;
    double randomWeight;
    double NegativeWeight;
    chance =rand()%2;//randomise between negative and positive numbers
      
    if(chance ==0){ 
		return drand48()/10000;
    }
     
    if(chance ==1){
     
		return -drand48()/10000;
    }
     
}
 

void GeneticAlgorithm::InitilisePopulation(int population, int stringsize )
{
    int count = 0;

	for(int r=0; r <  population  ; r++)
		Chromozone.push_back(vector<double> ());//create matrix
   
	for(int row = 0; row < population  ; row++) { 
		for(int col = 0; col < stringsize ; col++) 
			Chromozone[row].push_back(RandomWeights());//initialize with randome weights
	}

	for(int r=0; r <  population  ; r++)
		FinalChromozone.push_back(vector<double> ());//create matrix
   
	for(int row = 0; row < population  ; row++) { 
		for(int col = 0; col < stringsize ; col++) 
			FinalChromozone[row].push_back(0);//initialise with 0s for each row
	}
  
    for(int r=0; r <  population  ; r++)
		FirstChromozone.push_back(vector<double> ());//create matrix
   
	for(int row = 0; row < population  ; row++) { 
		for(int col = 0; col < stringsize ; col++) 
			FirstChromozone[row].push_back(RandomWeights());//intialize with random weights
	}
 
    for(int r=0; r <  population  ; r++)
		NewGeneration.push_back(vector<double> ());//create matrix to hold data
   
	for(int row = 0; row < population  ; row++) { 
		for(int col = 0; col < stringsize ; col++) 
			NewGeneration[row].push_back(0);//intialize newgeneration vector with 0s for all rows
	}
    
	for(int r=0; r <  population  ; r++)
		Fitness.push_back(0);//initialize fitness vector with 0s
         
    for(int r=0; r <  population  ; r++)
		FitnessRatio.push_back(0);  //intialize fitness_ratio vector with 0s 
      
}

 

double GeneticAlgorithm::FitnessFunc(Layer x, int ProbNum)
{

	int i,j,k;
    double z; 
	double fit = 0.0;  
	double   sumSCH; 

	if(ProbNum==1){
	// Ellipsoidal function
		for(j=0;j< x.size();j++)
			fit+=((j+1)*(x[j]*x[j]));
	}
	else if(ProbNum==2){
		// Schwefel's function
		for(j=0; j< x.size(); j++)
		{
			sumSCH=0;
			for(i=0; i<j; i++)
			sumSCH += x[i];
			fit += sumSCH * sumSCH;
		}
	}
	 
	return  1/fit;
}
 
void GeneticAlgorithm::Print()
{
	
	for(int row = 0; row < Population  ; row++) { 
		for(int col = 0; col < Stringsize ; col++) 
			cout<<Chromozone[row][col]<<" ";//output all the values of the population
			cout<<"       "<<1.0/ FitnessFunc(Chromozone[row],function );
			cout<<endl;
	}
   
	cout<<endl; 
	cout<<endl;
	cout<<endl;
	cout<<"-----------------------------------------"<<endl;
	cout<<" Fitness: "<<endl<<endl;
    for(int r=0; r <  Population  ; r++)
		cout<<Fitness[r]<<" ";//output the fitness for the entire population
		cout<<endl;   
	
	cout<<" Fitness Ratio: "<<endl<<endl;
    for(int r=0; r <  Population  ; r++)
		cout<<FitnessRatio[r]<<" ";  //output the fitness ration for the entire ratio
		cout<<endl;         
               
               
}

double GeneticAlgorithm::Evaluate() 
{
    double sum =0;
           
    for(int label = 0; label < Population  ; label++) {           
        Fitness[label] =    FitnessFunc(Chromozone[label], function );
  
        sum+=Fitness[label];
	}
          
    for(int label = 0; label < Population  ; label++) { 
        FitnessRatio[label] = Fitness[label]/sum*100;
	}
                  
}      
      

 double GeneticAlgorithm:: MaxFitness(){
    double max = 0;
        
    for(int label = 0; label <Population  ; label++) { 
        if( (Fitness[label]) >   max ){ 
			max = Fitness[label]; 
			
		}
  
    }
        
    return max;
}


int GeneticAlgorithm:: MaxLocation(){
        
    for(int label = 0; label <Population  ; label++) { 
        if( MaxFitness() ==Fitness[label]  ) 
			return label; 
    }
        
        
    return 0;
}    

int GeneticAlgorithm:: Select()
{
    Nodes Wheel(Population+1);
    double random;
    random = rand()%100;
    double sum = 0;   
    Wheel[0] = 0;
    for(int label = 1; label < Wheel.size()  ; label++) { 
        Wheel[label]  =    (FitnessRatio[label-1]+Wheel[label-1]);
	
	}
         
	for(int label = 0; label <(Wheel.size()-1 ) ; label++)  {
        if( (random >=Wheel[label])&&(random < Wheel[label+1])) {  
           
			return label;
        }
	}
}
  

void GeneticAlgorithm:: CrossoverAndMutate(int leftpair,int rightpair,double Crossover,double Mutation, int position)
{
                 
	Layer LeftChromozone(Stringsize);                           
	Layer RightChromozone(Stringsize);
	Nodes ChildChromozone(Stringsize);  
	Nodes Temporary(Stringsize);
	
	double leftfit = 0;
	double rightfit = 0;
	double childfit = 0;
	double mutatechildfit = 0;
   
	for(int gene = 0; gene < Stringsize ; gene++){
		LeftChromozone[gene] = Chromozone[leftpair][gene];
		RightChromozone[gene] = Chromozone[rightpair][gene];
	}

	leftfit = Fitness[leftpair];
	rightfit = Fitness[rightpair];


	double random = rand() % 100 ;
	int chooseparent  ;
	int choosegene  ;
	bool cross = false;      
	double newfit = 0;
	double  Lfit = 0;
	double Rfit = 0;
	
	if( random < (100 *Crossover)){
	
		cross = true;
		double randomalpha = drand48();
		Lfit = FitnessFunc(LeftChromozone, function) ;
		Rfit = FitnessFunc(RightChromozone, function) ;
		
		//Wrights Heuristic Crossover
		if ( Lfit > Rfit){
			for(int gene = 0; gene < Stringsize ; gene++)
				ChildChromozone[gene] = randomalpha *((LeftChromozone[gene]-RightChromozone[gene]))+ LeftChromozone[gene];
		}else{
			for(int gene = 0; gene < Stringsize ; gene++)
				ChildChromozone[gene] = randomalpha *((RightChromozone[gene]-LeftChromozone[gene]))+ RightChromozone[gene];   
		}

	}
	else{ 
		ChildChromozone= LeftChromozone; 

	}
	int ch = rand()%100;

	
	double randommutate;
	int seed = rand()%Stringsize;
	double MutationNonUniform = 0.1;
	int chance = rand()%1000;
	double Max = 0;
    int newseed = rand()%Stringsize;
	
	
	if( chance < (1000 *MutationNonUniform)){ 
		
		for(int gene = 0 ; gene < Stringsize ; gene++)
			ChildChromozone[ gene ] +=  RandomAddition();
      
    } 
      
    for(int gene = 0; gene < Stringsize ; gene++){
		NewGeneration[position][gene] = ChildChromozone[gene]; 
    }  
}
 
     
int GeneticAlgorithm:: Algorithmn(double Crossover,double Mutation, ofstream &output1, ofstream &output2 )   
{
	
	
	clock_t start = clock();
        int Gen = 1;
	double maxfitness = 0;                    
	InitilisePopulation( Population, Stringsize);                   
 					 
	Evaluate(); 
	maxfitness = MaxFitness();
  
   
	int leftpair = 0;
	int rightpair = 0; 
	bool stop = false;
	double train = 0;     
	 	
	
	
        Success = false;
        
      
        int FE = 0;

	while( FE < MaxFE){      
			 
		Evaluate();                     
		maxfitness =  MaxFitness(); 

		chrome = 0;
		int count =0;
		
		while((count) < Population ) {  
			do{
				leftpair = Select() ;
				rightpair = Select() ;
			} while (leftpair == rightpair);
	   
			CrossoverAndMutate(leftpair,rightpair,Crossover, Mutation,count);
   
			count ++; 
		  
		}
		count=0;
		
		 
			
		Chromozone = NewGeneration; 
 
		 Error = 1/maxfitness;
		 if (Error < 1E-10){
			 Success = true;
			 break; 
		 }

                 FE = Gen*Population; //Function Eval
		 
		  if (Gen % 1000 ==0){
		  output1<<Error<<" "<<FE<<endl;} // output convergence trend to file
		 
                 Gen++;

	}
	
             output1 <<endl<<endl; 

   clock_t finish = clock();

          Cycles = ((double)(finish - start))/CLOCKS_PER_SEC;
          cout<<Cycles<<" ----"<<endl;

          TotalEval = FE;
   
	    

		for(int x=0; x<Stringsize; x++){		
			output2<<Chromozone[MaxLocation()][x]<< " ";
		}
		output2<<endl;
		
	output2<<"   ---    "<<1/maxfitness <<" " <<  FE<< "    "<<  Cycles <<endl;
  

    return Gen;                       
}
  
  
  



int main(void)
{ 
 
    double Crossover = 0.9; 
    int Dimension = 10; //  
      
    srand (time(NULL));

	ofstream out1;
 	out1.open("out1.txt");
 	ofstream out2;
 	out2.open("out2.txt");
	ofstream out3;
  	out3.open("out3.txt");  

	double Mutation = 0.05;
	
        int Population =  50; 
	

	 for( int function =1 ; function <= 2;  function++){

		         int MeanEval=0;
	 	       	 double MeanError=0;
	 	         double MeanTime=0;
 
                         int CountSuccess = 1;
		         Sizes EvalAverage;
	 	       	 Layer ErrorAverage;
	 	       	 Layer TimeAverage; 
	 	       	 
		         int maxRun = 10;// num of exp


		for( int run =1 ; run <= maxRun; run ++){
          		
			GeneticAlgorithm ga(Population,Dimension, function);  
     		        int FE = ga.Algorithmn(Crossover, Mutation, out1 ,out2);
      
                  if (ga.GetSucess()){ 
                       CountSuccess++; 
                       EvalAverage.push_back(ga.GetEval());
	               MeanEval+=ga.GetEval();

	               ErrorAverage.push_back(ga.GetError());
	               MeanError+= ga.GetError();
	           	  
	               TimeAverage.push_back(ga.GetCycle());
	               MeanTime+= ga.GetCycle();
                       } 
          
	         cout<<" problem-Num Run-Num FE Fitness: "<< function<<" "<<run<<" "<<ga.GetEval()<<" "<<ga.GetError()<<endl;
     
	             }
 
	     MeanEval=MeanEval/CountSuccess;
	     MeanError=MeanError/CountSuccess;
	     MeanTime=MeanTime/CountSuccess; 
              //write stats to file
	 out3<< "   "<<MeanEval<<"     "<<MeanError<<"    "<<MeanTime<<"       "<<Population<<"   "<<ErrorAverage.size()<<endl;
	 EvalAverage.empty(); 
	 ErrorAverage.empty(); 
	 TimeAverage.empty();  
	
    }
 
	out1.close();
	out2.close();
	out3.close(); 
 
 return 0;

} 


