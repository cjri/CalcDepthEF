//Program to process .sam file data

#include <iostream>
#include <vector>
#include <string>
#include <sstream>


using namespace std;

#include "shared.h"
#include "alignment.h"
#include "call_snps.h"
#include "ddups.h"
#include "io.h"
#include "matchpairs.h"
#include "optimisation.h"
#include "utilities_sam.h"


int main(int argc, const char **argv) {

	run_params p;
	GetOptions (p,argc,argv);
	
	int seed=p.seed;
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);

	ifstream col_file;
	col_file.open("column.in");
	int n;
	col_file >> n;
	p.plines=n;
	
	
	cout << "Processing .sam files\n";
	
	//Import reference sequence
	ifstream ref_file;
	ref_file.open(p.ref.c_str());
	if (p.ref=="") {
		cout << "Error: No reference sequence specified.  Use e.g. --ref reference.fa.  Note that the sequence must be contained on a single line of the FASTA file\n";
		return 0;
	}
	rseq refseq;
	GetRefSeq (ref_file,refseq);
	string ref=refseq.seq;
		
	//Import list of .sam file names
	vector<string> sam_files;
	ImportSamFileNames(p,sam_files);
	if (sam_files.size()==0) {
		cout << "Error: No input .sam files specified.  Specify these files in Input_files.in\n";
		return 0;
	}
		
	//Construct quality reference vector
	vector<char> qual;
	makequal(p,qual);

	//Process files one by one...
	
	ofstream tim_file;
	tim_file.open("Times.in");
	for (int i=0;i<sam_files.size();i++) {
		tim_file << i << "\n";
		//Read .sam file
		alldat a;
		ImportSamFile (p,i,refseq,sam_files,a);
			
		//Detect and remove duplicate reads
		int found_pairs=0;
		if (p.ddup>=0) {
			DelDupSequences (p,found_pairs,qual,a);
		}

		AlignSequencesSam (p,a.s_length,qual,refseq,a.data);

		//Identify pairs
		if (p.pairs==1) {
			cout << "Match paired-end reads\n";
			if (found_pairs==0) {
				MatchPairs(p,1,a.data);
			} else {
				CheckOldPairs(a.data);
			}
		}
		
		//Join pairs if they exist
		vector<joined> jseqs;
		JoinPairs (i,p,jseqs,a.data);

	}
	tim_file.close();
	cout << "\nFind Single Locus polymorphisms\n";
	

	//Load sam file data
//	cout << sam_files.size() << "\n";
	
	
	//Get joined and aligned sequence reads
	vector<nuc> ref_counts;
	for (int i=0;i<sam_files.size();i++) {
		vector<joined> t_read;
		InputJnData (i,t_read);
	
		//Construct individual allele frequencies
		nuc r_count;
		CountNucs (p,refseq,t_read,r_count);
		ref_counts.push_back(r_count);
		
		//Output variant file
		OutputVarFile (i,refseq,r_count);
	}
	
	//Identify frequencies above a given cutoff
	vector<poly> polys;
	CallPolymorphisms (p,refseq,ref_counts,polys);
	
	//Construct single-locus tranjectories.  Contains locus, times, four nucleotide counts
	vector<str> sltrajs;
	ConstructSLTrajs(p,polys,ref_counts,sltrajs);
	
	//Output sltrajs data
	OutputSLTData(p,p.out_file.c_str(),sltrajs);
	
	
	
	
	cout << "\nCalculate noise parameter\n";
	
	sltrajs.clear();
	p.in_file="Single_locus_trajectories.out";
	
	ImportSLTData(p,sltrajs);  //N.B. Times are encoded in the trajectories
	
	SLTFreqs(p,sltrajs); //Convert observations to frequencies
	
	//Remove trajectories that move by more than a cutoff amount per day
	if (sltrajs[0].times.size()>1) {
		FilterSLTrajs(p,sltrajs);
	}
	
	//Remove non-polymorphic time-points from trajectories
	FilterSLTrajs2(p,sltrajs);
	
	//Calculate mean frequencies of trajectories - approximate fit under assumption of neutrality
	SLTMeanFreqs(p,sltrajs);
	
	//Calculate vector of log factorials
	vector<double> fact_store;
	GetFactVectorSL(sltrajs,fact_store);
	
	double Csl_opt=0;
	OptimiseCsl(p,Csl_opt,sltrajs,fact_store,rgen);
	
	ofstream csl_file;
	csl_file.open("Csl.out");
	csl_file << Csl_opt << "\n";
	csl_file.close();
	

	
	cout << "\nCalculate effective depths\n";
	
	//Read in Variant data
	vector< vector<int> > all_tots;
	GetVariantTotals(sam_files,all_tots);
	
	//Convert depth to effective depth
	ifstream csli_file;
	double c;
	csli_file.open("Csl.out");
	csli_file >> c;
	vector< vector<int> > all_etots;
	for (int i=0;i<all_tots.size();i++) {
		vector<int> et;
		for (int j=0;j<all_tots[i].size();j++) {
			double cd=(all_tots[i][j]*(1+c))/(all_tots[i][j]+c);
			int cdi=floor(cd);
			et.push_back(cdi);
		}
		cout << i << " " << et.size() << "\n";
		all_etots.push_back(et);
	}
		
		
	for (int i=0;i<sam_files.size();i++) {
		ofstream dep_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name = "Depths"+temp+".out";
		dep_file.open(name.c_str());
		for (int j=0;j<all_tots[i].size();j++) {
			dep_file << j << " " << all_tots[i][j] << " " <<  all_etots[i][j] << "\n";
		}
	}
	
	return 0;

}
							
 




