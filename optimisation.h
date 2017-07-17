//Shared information for linked optimisation
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

using namespace std;

void GetFactVectorSL (vector<str> sltrajs, vector<double>& fact_store);
void GetFactVectorML (vector< vector<mtr> > mltrajs, vector<double>& fact_store);
void FindLogFact(vector<double>& fact_store, int N);
void OptimiseCsl (run_params p, double& Csl_opt, vector<str> sltrajs, vector<double> fact_store, gsl_rng *rgen);

