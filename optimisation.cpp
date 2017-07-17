#include "shared.h"
#include "optimisation.h"
#include "likelihood.h"
#include <iostream>
#include <string>

void GetFactVectorSL (vector<str> sltrajs, vector<double>& fact_store) {
	int maxN=0;
	for (int i=0;i<sltrajs.size();i++) {
		
		for (int j=0;j<sltrajs[i].qA.size();j++) {
			if (maxN<sltrajs[i].nN[j]) {
				maxN=sltrajs[i].nN[j];
			}
		}
	}
	maxN=maxN*2;
	cout << maxN << "\n";
	FindLogFact(fact_store,maxN);
}

void GetFactVectorML (vector< vector<mtr> > mltrajs, vector<double>& fact_store) {
	int maxN=0;
	for (int i=0;i<mltrajs.size();i++) {
		for (int j=0;j<mltrajs[i].size();j++) {
			int tot=0;
			for (int k=0;k<mltrajs[i][j].n.size();k++) {
				tot=tot+mltrajs[i][j].n[k];
			}
			if (maxN<tot) {
				maxN=tot;
			}
		}
	}
	maxN=maxN+10;
	FindLogFact(fact_store,maxN);
}

void FindLogFact(vector<double>& fact_store,int N){
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=N;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
	}
}


void OptimiseCsl (run_params p, double& Csl_opt, vector<str> sltrajs, vector<double> fact_store, gsl_rng *rgen) {
	cout << "Find optimal value of C_sl\n";
	Csl_opt=100;
	double C=100;
	double logL=-1e10;
	double logL_store=-1e10;
	double changec=10;
	int tryc=0;
	int acceptc=0;
	int max_it=100000;
	int first=1;
	vector<double> progress;
	
	for (int it=0;it<max_it;it++) {
		//cout << "Iteration " << it << "\n";
		if (it==max_it-1) {
			cout << "Warning: Possible convergence failure in C\n";
		}
		//cout << "Iteration " << it << "\n";
		if (first==0) {
			if (logL-logL_store>0) {
				Csl_opt=C;
				if (p.verb==1) {
					cout << "C " << C << " logL " << logL << "\n";
				}
				logL_store=logL;
				progress.push_back(logL);
				acceptc++;
			} else {
				C=Csl_opt;
			}
		}
		first=0;
		if (progress.size()>11&&progress[progress.size()-1]-progress[progress.size()-10]<0.001) {
			break;
		}
		
		if (it%100==0&&it>0) {
			double a_rate=(acceptc+0.)/(tryc+0.);
			changec=changec*(0.95+a_rate);
			acceptc=0;
			tryc=0;
		}
		
		C=C+(gsl_rng_uniform(rgen)*changec)-(changec/2);
		tryc++;
		if (it==max_it-1) {
			C=Csl_opt;
		}
		
		//Calculate likelihood
		logL=0;
		for (int i=0;i<sltrajs.size();i++) {
			if (sltrajs[i].inc==1) {
				vector<double> inf;
				vector<int> inc (4,0);
				if (sltrajs[i].mA>0) {
					inf.push_back(sltrajs[i].mA);
					inc[0]=1;
				}
				if (sltrajs[i].mC>0) {
					inf.push_back(sltrajs[i].mC);
					inc[1]=1;
				}
				if (sltrajs[i].mG>0) {
					inf.push_back(sltrajs[i].mG);
					inc[2]=1;
				}
				if (sltrajs[i].mT>0) {
					inf.push_back(sltrajs[i].mT);
					inc[3]=1;
				}
				for (int j=0;j<sltrajs[i].qA.size();j++) {
					vector<int> obs;
					int N=sltrajs[i].nN[j];
					if (N>=p.dep_cut) {
						if (inc[0]==1) {
							obs.push_back(sltrajs[i].nA[j]);
						}
						if (inc[1]==1) {
							obs.push_back(sltrajs[i].nC[j]);
						}
						if (inc[2]==1) {
							obs.push_back(sltrajs[i].nG[j]);
						}
						if (inc[3]==1) {
							obs.push_back(sltrajs[i].nT[j]);
						}
						logL=logL+DirichletMultiCalc(N,C,obs,inf,fact_store);
					}
				}
			}
		}
	}
	
	cout << "Optimal C is " << Csl_opt << "\n";

}

