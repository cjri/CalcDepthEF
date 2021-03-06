#include "shared.h"
#include "alignment.h"
#include "utilities_sam.h"
#include <iostream>
#include <sstream>
#include <string>

void AlignSequencesSam (run_params& p, int s_length, vector<char> qual, rseq refseq, vector<rd>& data) {
	int q=0;
	//cout << "Length " << s_length << "\n";
	string q0;
	int xalq=0;
	MinBaseQual(qual,q0); //Find code for minimum base quality
	for (int i=0;i<s_length;i++) {
		if (p.verb==1) {
			cout << "Sequences " << i << "\n";
			cout << data[i].seq << "\n";
			cout << data[i].qual << "\n";
		}
		data[i].inc=0;
		if (data[i].del==0) { //Not previously deleted
			//cout << data[i].ref << " " << refseq.name << "\n";
			if (data[i].ref==refseq.name) {
				ReadCigar(i,data); //Read CIGAR string
				q=findqual(p,i,p.min_qual,p.max_qual,qual,data);  //Assess read by median nucleotide quality
				//cout << "Median qual " << q << "\n";
				if (q>=p.min_qual) {  //Minimum sequence quality test
					//cout << "Align qual " << data[i].alq << " " << p.ali_qual << "\n";
					if (data[i].alq>=p.ali_qual) { //Minumum alignment quality
						if (p.ali_inc==1||(p.ali_inc==0&&data[i].alq!=255)) { //Minumum alignment quality
							//cout << data[i].alpos << " " << refseq.seq.size()-p.min_rlen << "\n";
							if (data[i].alpos<refseq.seq.size()-p.min_rlen) { //Sequence must be aligned to have an overlap of at least the minimum number of alleles reported by a read; this weeds out misaligned sequences at the end of the reference sequence.
								//cout << "Included\n";
								data[i].inc=1;
								RemoveInitialSoftClipping(i,data);  //Remove initial soft clipping
								FixDeletions(i,q0,data);  //Include blank data for deletions
								FixInsertions(i,data);  //Remove insertions
								RemoveSoftClipping(i,data);  //Remove remaining soft clipping
								ProcessReadQual (i,p,qual,data); //Process data by individual nucleotide quality
							}
						} else if (data[i].alq==255) {
							xalq++;
						}
					}
				}
			}
		}
	}
	if (xalq>200) {
		cout << "Warning: " << xalq << " sequences excluded due to unknown alignment quality\n";
	}
}


void MinBaseQual (vector<char> qual, string& q0) {
	stringstream ss;
	char c = qual[0];
	ss << c;
	ss >> q0;
}

//Read and process CIGAR string
void ReadCigar (int i, vector<rd>& data) {
	int pos=0;
	for (int j=0;j<data[i].cigar.size();j++) {
		if (data[i].cigar.compare(j,1,"M")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "M " << i_num << "\n";
			pos=j+1;
			string app (i_num,'M');
			data[i].cig_string.append(app);
		}
		if (data[i].cigar.compare(j,1,"S")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "S " << i_num << "\n";
			pos=j+1;
			string app (i_num,'S');
			data[i].cig_string.append(app);
		}
		
		if (data[i].cigar.compare(j,1,"H")==0) {
			pos=j+1;
		}

		if (data[i].cigar.compare(j,1,"I")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "I " << i_num << "\n";
			pos=j+1;
			string app (i_num,'I');
			data[i].cig_string.append(app);
		}
		if (data[i].cigar.compare(j,1,"D")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "D " << i_num << "\n";
			pos=j+1;
			string app (i_num,'D');
			data[i].cig_string.append(app);
		}
		if (data[i].cigar.compare(j,1,"P")==0) {  //Padding
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "P " << i_num << "\n";
			pos=j+1;
			string app (i_num,'P');
			data[i].cig_string.append(app);
		}
	}
	FilterCigar(i,data);
}

//Remove padding from cigar string
void FilterCigar (int i, vector<rd>& data) {
	int j=0;
	while (j<data[i].cig_string.length()) {
		if (data[i].cig_string.compare(j,1,"P")==0) {
			data[i].cig_string.erase(j,1);
		} else {
			j++;
		}
	}
}

//Check median quality of sequence.  Reduce sequence to achieve median quality if required.  Reduction is carried out from both ends; get the longest qualifying sequence
int findqual (run_params p, int i, int min_qual, int max_qual, vector<char> qual, vector<rd>& data) {
	string q=data[i].qual;
	string s=data[i].seq;
	string c=data[i].cig_string;
	//	cout << s << "\n";
	//	cout << q << "\n";
	vector<int> qvec;
	for (int i=0;i<s.size();i++) {
		int done=0;
		for (int j=0;j<=max_qual;j++) {
			if (q[i]==qual[j]) {
				qvec.push_back(j);
				done=1;
				break;
			}
		}
		if (done==0) {
			cout << "Error: No match to base quality score: " << q[i] << "\n";
			qvec.push_back(0);
		}
	}
	
	int median = GetMedian(0,0,qvec);
	//cout << "Median " << median << " ";
	int qo=median;
	
	if (qvec.size()>=p.min_rlen) {
		if (median<min_qual) {  //Edit sequence to get high quality part
			int s1=999;
			int s2=999;
			//Reduce sequence by removing nucleotides from the end of the read
			for (int a=1;a<qvec.size()-p.min_rlen;a++) {
				median=GetMedian(a,0,qvec);
				if (median==min_qual) {
					s1=a;
					break;
				}
			}
			for (int a=1;a<qvec.size()-p.min_rlen;a++) {
				median=GetMedian(0,a,qvec);
				if (median==min_qual) {
					s2=a;
					break;
				}
			}
			int qs=q.size();
			
			if (s1<999&&s1<=s2) {
				q=q.substr(s1,qs-s1+1);
				s=s.substr(s1,qs-s1+1);
				if (c.length()>0) {
					c=c.substr(s1,qs-s1+1);
				}
				data[i].qual=q;
				data[i].seq=s;
				data[i].cig_string=c;
				qo=min_qual;
			} else if (s2<999&&s2<s1) {
				q=q.substr(0,qs-s2+1);
				s=s.substr(0,qs-s2+1);
				if (c.length()>0) {
					c=c.substr(0,qs-s2+1);
				}
				data[i].qual=q;
				data[i].seq=s;
				data[i].cig_string=c;
				qo=min_qual;
			} else {
				qo=-1;
			}
		}
	} else {
		qo=-1;
	}
	//	cout << "Quality " << qo << "\n";
	
	if (qo>-1) { //Check number of minimum quality reads
		int nm=0;
		for (int i=0;i<s.size();i++) {
			for (int j=min_qual;j<=max_qual;j++) {
				if (q[i]==qual[j]) {
					nm++;
					break;
				}
			}
		}
		if (nm<p.min_rlen) {
			qo=-1;
		}
	}
	return qo;
}

int GetMedian (int a, int b, vector<int> qvec) {
	vector<int> qvec2=qvec;
	sort(qvec2.begin()+a,qvec2.end()-b);
	double median;
	double size=qvec2.size();
	if (qvec2.size()==0) {
		median=(qvec2[size/2-1]+qvec2[size/2])/2;
	} else {
		median = qvec2[size/2];
	}
	return median;
}


void RemoveInitialSoftClipping (int i, vector<rd>& data) {
	if (data[i].cig_string.compare(0,1,"S")==0) {
		string s;
		string q;
		string c;
		int j=0;
		while (data[i].cig_string.compare(j,1,"S")==0) {
			j++;
		}
		data[i].seq=data[i].seq.substr(j,data[i].seq.length()-j);
		data[i].qual=data[i].qual.substr(j,data[i].qual.length()-j);
		data[i].cig_string=data[i].cig_string.substr(j,data[i].cig_string.length()-j);
	}
}

void FixDeletions (int i, string q0, vector<rd>& data) {
	for (int j=0;j<data[i].cig_string.length();j++) {
		if (data[i].cig_string.compare(j,1,"D")==0) {
			data[i].seq.insert(j,"-");
			data[i].qual.insert(j,q0);
		}
	}
}

void FixInsertions (int i, vector<rd>& data) {
	for (int j=0;j<data[i].cig_string.length();j++) {
		if (data[i].cig_string.compare(j,1,"I")==0) {
			while (data[i].cig_string.compare(j,1,"I")==0&&j<data[i].cig_string.length()) {
				data[i].seq.erase(j,1);
				data[i].qual.erase(j,1);
				data[i].cig_string.erase(j,1);
			}
		}
	}
}

void RemoveSoftClipping (int i, vector<rd>& data) {
	for (int j=0;j<data[i].cig_string.length();j++) {
		if (data[i].cig_string.compare(j,1,"S")==0) {
			while (j<data[i].cig_string.length()&&data[i].cig_string.compare(j,1,"S")==0) {
				data[i].seq.erase(j,1);
				data[i].qual.erase(j,1);
				data[i].cig_string.erase(j,1);
			}
		}
	}
}

//Process short read, removing alleles that do not have sufficient nucleotide quality
void ProcessReadQual (int i, run_params p, vector<char> qual, vector<rd>& data) {
	int incl=0;
	for (int k=0;k<data[i].seq.size();k++) {
		int keep=0;
		for (int j=p.min_qual;j<=p.max_qual;j++) {
			if (data[i].qual[k]==qual[j]) {
				keep=1;
				break;
			}
		}
		if (keep==0) {
			data[i].seq.replace(k,1,"-");
		} else {
			incl++;
		}
	}
	if (incl<p.min_rlen) {
		data[i].inc=0;
	}
}
