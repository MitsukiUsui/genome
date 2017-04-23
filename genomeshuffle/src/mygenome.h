#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <random>
#include <cassert>

using namespace std;

vector<string> split(string const &str, char sep)
{
	vector<string> v;
	stringstream ss(str);
	string buffer;
	while (getline(ss, buffer, sep)) {
		v.push_back(buffer);
	}
	return v;
}

template <class Content> bool is_in(Content c, vector<Content> &v){
	for (int i = 0; i < v.size(); i++){
		if (c == v[i]){
			return true;
		}
	}
	return false;
}

string reverse_complement(string &s){
	string before = "AaCcGgTt";
	string after = "TtGgCcAa";
	assert(before.length() == after.length());

	string r;
	for (int i = s.length() - 1; i >= 0; i--){
		char c = s[i];
		bool found = false;
		for (int j = 0; j<before.length(); j++){
			if (before[j] == c){
				found = true;
				r += after[j];
				break;
			}
		}
		if (!found){
			r += 'N';
		}
	}

	return r;
}

class GeneticCode{
	string aa_str = "*ACDEFGHIKLMNPQRSTVWY";//length 1+20
	int codonToAa[64];
	vector<int> aaToCodon[21];

	int codonCount[64];
	double codonFreq[64];

	int codon_encode(string codon){
		if (codon.length() != 3){
			cerr << "The length of codon is not 3." << endl;
			exit(1);
		}

		int sum = 0;
		for (int i = 0; i < 3; i++){
			int base = 1;
			for (int j = i; j < 2; j++){
				base *= 4;
			}

			switch (codon[i]){
			case 'T':
				sum += 0 * base;
				break;
			case 'C':
				sum += 1 * base;
				break;
			case 'A':
				sum += 2 * base;
				break;
			case 'G':
				sum += 3 * base;
				break;
			default:
				cerr << "Unusual nucleotide " << codon[i] << " is included." << endl;
				exit(1);
			}
		}
		return sum;
	}
	string codon_decode(int id){
		if (id < 0 | id > 63){
			cerr << "The id " << id << " is not in the range of 0 ~ 64." << endl;
			exit(1);
		}

		int r_int[3];
		string r;

		r_int[0] = id / 16;
		r_int[1] = (id % 16) / 4;
		r_int[2] = (id % 16) % 4;

		for (int i = 0; i < 3; i++){
			switch (r_int[i]){
			case 0:
				r += 'T';
				break;
			case 1:
				r += 'C';
				break;
			case 2:
				r += 'A';
				break;
			case 3:
				r += 'G';
				break;
			}
		}
		return r;
	}

public:
	string code_str;//length 64

	GeneticCode(string str){
		code_str = str;
		if (code_str.length() != 64){
			cerr << "The length of genetic code is not 64." << endl;
			exit(1);
		}

		//create codonToAa and aaToCodon;
		for (int codonId = 0; codonId < 64; codonId++){
			char c = code_str[codonId];
			int aaId = -1;
			for (int i = 0; i < 21; i++){
				if (c == aa_str[i]){
					aaId = i;
					break;
				}
			}

			if (aaId == -1){
				cerr << "The genetic code includes unusual Amino Acid \""<<c<< " \"."<< endl;
				exit(1);
			}

			codonToAa[codonId] = aaId;
			aaToCodon[aaId].push_back(codonId);
		}
	}

	~GeneticCode(){
	}

	void clear_count(){
		for (int i = 0; i < 64; i++){
			codonCount[i] = 0;
		}
		return;
	}

	void update_count(string str){
		if (str.length() % 3 != 0){//str must be multiple of 3
			cerr << "The length of str is not multiple of3." << endl;
			exit(1);
		}

		int loopCount = str.length() / 3;
		for (int i = 0; i < loopCount; i++){
			string codon = str.substr(3 * i, 3);
			int codonId = codon_encode(codon);
			codonCount[codonId] += 1;
		}
	}

	//update codonFreq
	void convert_to_freq(){
		for (int aaId = 0; aaId < 21; aaId++){
			int size = aaToCodon[aaId].size();

			int sum = 0;
			for (int i = 0; i < size; i++){
				int codonId = aaToCodon[aaId][i];
				sum += codonCount[codonId];
			}

			if (sum>0){
				for (int i = 0; i < size; i++){
					int codonId = aaToCodon[aaId][i];
					codonFreq[codonId] = (double)codonCount[codonId] / sum;
				}
			}
			else if(sum==0){
				for (int i = 0; i < size; i++){
					int codonId = aaToCodon[aaId][i];
					codonFreq[codonId] = 0;
				}
			}
		}
	}

	//true...freq, false...count;
	void show_freq(bool mode){
		for (int aaId = 0; aaId < 21; aaId++){
			cout << aa_str[aaId] << "   ";
			int size = aaToCodon[aaId].size();
			for (int i = 0; i < size; i++){
				int codonId = aaToCodon[aaId][i];
				if (mode){
					cout << codon_decode(codonId) << ":" << codonFreq[codonId] << ", ";
				}
				else{
					cout << codon_decode(codonId) << ":" << codonCount[codonId] << ", ";
				}
			}
			cout << endl;
		}
		cout << endl;
	}

	bool is_stop_codon(string codon){
		int aaId = codonToAa[codon_encode(codon)];
		if (aaId == 0){
			return true;
		}
		else{
			return false;
		}
	}


	//return synonymous codon according to codonFreq
	string synonymous_sub(string codon, mt19937 &mt){
		int aaId = codonToAa[codon_encode(codon)];

		uniform_real_distribution<double> dist(0.0, 1.0);
		double d = dist(mt);

		double cumsum = 0;
		int codonId = -1;
		for (int i = 0; i < aaToCodon[aaId].size(); i++){
			codonId = aaToCodon[aaId][i];
			cumsum += codonFreq[codonId];
			if (cumsum > d){
				break;
			}
		}
		return codon_decode(codonId);
	}

	string translate(string &str){
		int length = str.length();
		if (length % 3 != 0){
			cerr << "string to be translated is not multiple of 3." << endl;
			exit(1);
		}
		int aaLength = length / 3;

		string r;
		for (int i = 0; i < aaLength; i++){
			string codon = str.substr(3 * i, 3);
			int codonId = codon_encode(codon);
			int aaId = codonToAa[codonId];
			r += aa_str[aaId];
		}
		return r;
	}
	
};
