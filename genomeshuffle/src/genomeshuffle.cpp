#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <numeric>
#include <map>
#include <sstream>
#include <cstdlib>
#include <cmath>
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

//vectorÇÃíÜêgÇ…ä÷ÇÌÇÁÇ∏åüçıÇ≈Ç´ÇÈÇÊÇ§Ç…ÇµÇΩÇ¢
bool is_in(string &s, vector<string> &v){
	for (int i = 0; i < v.size(); i++){
		if (s == v[i]){
			return true;
		}
	}
	return false;
}

string reverse_complement(string &s){
	string r;
	for (int i = s.length() - 1; i >= 0; i--){
		char c;
		switch (s[i]){
		case 'A':
			c = 'T';
			break;
		case 'C':
			c = 'G';
			break;
		case 'G':
			c = 'C';
			break;
		case 'T':
			c = 'A';
			break;
		default:
			cerr << "Not a normal base." << endl;
			exit(1);
		}
		r += c;
	}
	return r;
}


class GeneticCode{
	pair<int, int> geneticCode[64];
	vector<pair<int, int> > codonCount[21];//first: codonId, second: count
	vector<pair<int, double> > codonFreq[21];//first: codonId, second: frequency

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
			cerr << "The id " << id << " is not in the range of 0~64" << endl;
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
	string code_str;
	string aa = "*ACDEFGHIKLMNPQRSTVWY";//length 1+20

	GeneticCode(string str){
		code_str = str;
		if (code_str.length() != 64){
			cerr << "The length of genetic code is not 64." << endl;
			exit(1);
		}

		//create geneticCode
		for (int codonId = 0; codonId < 64; codonId++){//initialization
			geneticCode[codonId] = make_pair(-1, -1);
		}
		for (int aaId = 0; aaId < 21; aaId++){
			char c = aa[aaId];
			int count = 0;
			for (int codonId = 0; codonId < 64; codonId++){
				if (code_str[codonId] == c){
					geneticCode[codonId] = make_pair(aaId, count++);
					codonCount[aaId].push_back(make_pair(codonId, 0));
					codonFreq[aaId].push_back(make_pair(codonId, 0));
				}
			}
		}
		for (int codonId = 0; codonId < 64; codonId++){//check error
			if (geneticCode[codonId].first == -1 | geneticCode[codonId].second == -1){
				cerr << "Unusual amino acid " << code_str[codonId] << " is included." << endl;
				exit(2);
			}
		}
	}

	~GeneticCode(){

	}

	//str must be multiple of 3
	void update_count(string str){
		if (str.length() % 3 != 0){
			cerr << "The length of str is not multiple of3." << endl;
			exit(1);
		}

		int loopCount = str.length() / 3;
		for (int i = 0; i < loopCount; i++){
			string codon = str.substr(3 * i, 3);
			int codonId = codon_encode(codon);
			pair<int, int> p = geneticCode[codonId];
			codonCount[p.first][p.second].second += 1;
		}
	}

	//update codonFreq
	void convert_to_freq(){
		for (int aaId = 0; aaId < 21; aaId++){
			int size = codonCount[aaId].size();

			int sum = 0;
			for (int i = 0; i < size; i++){
				sum += codonCount[aaId][i].second;
			}

			for (int i = 0; i < size; i++){
				codonFreq[aaId][i].second = (double)codonCount[aaId][i].second / sum;
			}
		}
	}

	//return synonymous codon according to codonFreq
	string synonymous_sub(string codon, mt19937 &mt){
		pair<int, int> p = geneticCode[codon_encode(codon)];
		int aaId = p.first;

		uniform_real_distribution<double> dist(0.0, 1.0);
		double d = dist(mt);

		double cumsum = 0;
		int codonId;
		for (int i = 0; i < codonFreq[aaId].size(); i++){
			cumsum += codonFreq[aaId][i].second;
			if (cumsum > d){
				codonId = codonFreq[aaId][i].first;
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
			r += aa[geneticCode[codon_encode(codon)].first];
		}

		return r;
	}

	void show_freq(){
		for (int aaId = 0; aaId < 21; aaId++){
			cout << aa[aaId] << "   ";
			int size = codonFreq[aaId].size();
			for (int i = 0; i < size; i++){
				cout << codon_decode(codonFreq[aaId][i].first) << ":" << codonFreq[aaId][i].second << ", ";
			}
			cout << endl;
		}
		cout << endl;
	}

	void show_count(){
		for (int aaId = 0; aaId < 21; aaId++){
			cout << aa[aaId] << "   ";
			for (int i = 0; i < codonCount[aaId].size(); i++){
				cout << codon_decode(codonCount[aaId][i].first) << ":" << codonCount[aaId][i].second << ", ";
			}
			cout << endl;
		}
		cout << endl;
	}

	void clear_count(){
		for (int aaId = 0; aaId < 21; aaId++){
			for (int i = 0; i < codonCount[aaId].size(); i++){
				codonCount[aaId][i].second = 0;
			}
		}
	}
};

class Region{
public:
	int start;
	int end;
	bool isForword;
	int type;

	static bool compare(Region r1, Region r2) { return (r1.start < r2.start); }

	Region(string str){
		vector<string> strSplited = split(str, ' ');
		start = stoi(strSplited[0]);
		end = stoi(strSplited[1]);
		isForword = (strSplited[2] == "+" ? true : false);
	}

	Region(int s, int e, bool b){
		start = s;
		end = e;
		isForword = b;
	}

	~Region(){
		return;
	}
};

class Regions{
	// 0... normal
	// 1... shorter than 6 bp
	// 2...not multiple of 3
	// 3...stop codon is inserted
	// 4...stop codon not exists
	// 5...start codon not exists
	int judge_type(string str, vector<string> &startCodon, vector<string> &stopCodon){
		int bpLength = str.length();
		int aaLength = str.length() / 3;

		if (bpLength < 6){//if shorter than 6 bp
			return 1;
		}
		if (bpLength % 3 != 0){//if not multiple of 3
			return 2;
		}


		for (int i = 1; i < aaLength - 1; i++){
			string codon = str.substr(3 * i, 3);
			if (is_in(codon, stopCodon)){//stop codon is inserted
				return 3;
			}
		}

		string codon = str.substr(3 * (aaLength - 1), 3);
		if (!is_in(codon, stopCodon)){
			return 4;
		}

		codon = str.substr(0, 3);
		if (!is_in(codon, startCodon)){
			return 5;
		}

		return 0;
	}

public:
	vector<Region> regions;


	//sentinal at end
	Regions(string &seq, string inputFilepath){
		ifstream ifs(inputFilepath);
		if (!ifs){
			cerr << "Cannot open inputfile (" << inputFilepath << ")." << endl;
			exit(1);
		}

		string buf;
		while (getline(ifs, buf)){
			Region r(buf);
			regions.push_back(r);
		}
		Region sentinal(seq.length() - 1, seq.length() - 1, true);
		regions.push_back(sentinal);

	}

	~Regions(){

	}

	bool check_inclusion(){
		bool r = true;

		sort(regions.begin(), regions.end(), Region::compare);
		//check if it has inclusion relation or not
		int tmp = -1;
		for (int i = 1; i < regions.size() - 1; i++){
			if (regions[i].end < tmp){
				cerr << "Inclusion relation occured." << endl;
				cerr << "region: " << regions[i - 1].start << "-" << regions[i - 1].end << ", type: " << regions[i - 1].type << " & " << "region: " << regions[i].start << "-" << regions[i].end << ", type: " << regions[i].type << endl;
				r = false;
			}
			tmp = regions[i].end;
		}
		return r;
	}

	void set_types(string &seq, vector<string> &startCodon, vector<string> &stopCodon){
		int size = regions.size();
		for (int i = 0; i < size; i++){
			Region r = regions[i];
			string tmp = seq.substr(r.start, r.end - r.start);
			string subSeq = (r.isForword ? tmp : reverse_complement(tmp));
			regions[i].type = judge_type(subSeq, startCodon, stopCodon);
		}
	}

	void show(){
		int size = regions.size();
		int count[6] = { 0, 0, 0, 0, 0, 0 };
		string exp[6] = { "Normal", "shorter than 6 bp", "not multiple of 3", "stop codon is inserted", "stop codon not exists", "start codon not exists" };

		for (int i = 0; i < size; i++){
			count[regions[i].type]++;
		}

		cout << "Total: " << size << " units." << endl;
		for (int i = 0; i < 6; i++){
			cout << "   Type " << i << " (" << exp[i] << "): " << count[i] << " units." << endl;
		}
		cout << endl;
		return;
	}
};




//shuffle synonymously str by the codon.
//str must be multiple of 3.
string shuffle_synonymous(string &str, bool isForword, mt19937 &mt, GeneticCode &gc){
	string s = (isForword ? str : reverse_complement(str));

	int aaLength = s.length() / 3;
	string tmp;
	for (int i = 0; i < aaLength; i++){
		string codon = s.substr(i * 3, 3);
		string codonNew = gc.synonymous_sub(codon, mt);
		tmp += codonNew;
	}
	assert(gc.translate(s) == gc.translate(tmp));

	string r = (isForword ? tmp : reverse_complement(tmp));
	return r;
}

//shuffle str by the codon
//str must be multiple of 3.
string shuffle_codon(string &str, mt19937 &mt){
	//create vector in order and randomize
	int aaLength = str.length() / 3;
	vector<int> v(aaLength);
	iota(v.begin(), v.end(), 0);
	vector<int>::iterator itr = v.begin();
	shuffle(v.begin(), v.end(), mt);

	//create shuffled gene 
	string r;
	for (int i = 0; i < aaLength; i++){
		r += str.substr(v[i] * 3, 3);
	}
	return r;
}

//shuffle str by the base
//str must be longer than 0 bp.
string shuffle_base(string &str, mt19937 &mt){
	//create vector in order and randomize
	int bpLength = str.length();
	vector<int> v(bpLength);
	iota(v.begin(), v.end(), 0);
	vector<int>::iterator itr = v.begin();
	shuffle(v.begin(), v.end(), mt);

	//create shuffled gene 
	string r;
	for (int i = 0; i < bpLength; i++){
		r += str[v[i]];
	}
	return r;
}


//shuffle region and return shuffled region length. when shuffle not occured return 0.
//mode=1,2, or 3
int shuffle_region(string &seq, int start, int end, bool isForword, int mode, mt19937 &mt, GeneticCode &gc){
	int length = end - start;

	if (mode == 1){
		if (length <= 0){
			return 0;
		}
	}
	else{//mode==2 or 3
		if (length % 3 != 0){
			cerr << "shuffle regin is not multiple of three." << endl;
			cerr << "start: " << start << ", end: " << end << endl;
			exit(1);
		}
		if (length <= 3){
			return 0;
		}
	}

	string subStr = seq.substr(start, length);
	string shuffled;
	switch (mode){
	case 1:
		shuffled = shuffle_base(subStr, mt);
		break;
	case 2:
		shuffled = shuffle_codon(subStr, mt);
		break;
	case 3:
		shuffled = shuffle_synonymous(subStr, isForword, mt, gc);
	}

	//overwrite
	for (int j = 0; j < length; j++){
		seq[start + j] = shuffled[j];
	}
	return length;
}


//shuffle genome according to interval end report
//called ones
void shuffle_genome(string &seq, Regions &cds, GeneticCode &gc, pair<int,int> mode){
	random_device rd;
	mt19937 mt(rd());
	int count_shuffle[4] = { 0, 0, 0, 0 };

	int endMax = 0;
	for (int i = 0; i < cds.regions.size() - 1; i++){//last element of cds.regions is sentinel
		int start = cds.regions[i].start;
		int end = cds.regions[i].end;
		int startNext = cds.regions[i + 1].start;
		assert(start <= startNext);//regions should be sorted

		bool isForword = cds.regions[i].isForword;
		int type = cds.regions[i].type;

		//shuffle intergenic region
		if (mode.first == 1){
			count_shuffle[1] += shuffle_region(seq, endMax, start, isForword, 1, mt, gc);
		}
		

		//shuffle genomic region
		if (type == 0){
			int cand1 = start + 3 * ceil((double)(endMax - start) / 3);
			int cand2 = start + 3;
			int shuffleStart = max(cand1, cand2);
			cand1 = end - 3 * ceil((double)(end - startNext) / 3);
			cand2 = end - 3;
			int shuffleEnd = min(cand1, cand2);
			count_shuffle[mode.second] += shuffle_region(seq, shuffleStart, shuffleEnd, isForword, mode.second, mt, gc);
		}

		//update maxEnd
		if (end > endMax){
			endMax = end;
		}
	}


	count_shuffle[0] = seq.length();
	for (int i = 1; i < 4; i++){
		count_shuffle[0] -= count_shuffle[i];
	}

	for (int i = 0; i < 4; i++){
		cout << "mode " << i << " shuffled : " << count_shuffle[i] << " bp." << endl;
	}
}


void count_orf_all(string &seq, vector<string> &stopCodon, vector<int> &count){
    int seqLen = seq.length();

    vector<string> compStopCodon;
    for (int i = 0; i<stopCodon.size(); i++){
        compStopCodon.push_back(reverse_complement(stopCodon[i]));
    }

    for (int lane = 1; lane <= 6; lane++){
        bool isForword = true;
        if (lane > 3){
            isForword = false;
        }

        int gap = (lane - 1) % 3;
        int numOfLoop = (seqLen - gap) / 3;

        int stopPos = -1;

        //find fist stopcodon
        int processNum = 0;
        for (; processNum < numOfLoop; processNum++){
            int pos = (isForword ? (seqLen - 3) - (3 * processNum + gap) : 3 * processNum + gap);//if 1~3, read backwords
            string codon = seq.substr(pos, 3);
            if (is_in(codon, (isForword ? stopCodon : compStopCodon))){
                stopPos = processNum;
                break;
            }
        }

        if (stopPos != -1){//if found first stop codon
            for (processNum++; processNum < numOfLoop; processNum++){
                int pos = (isForword ? (seqLen - 3) - (3 * processNum + gap) : 3 * processNum + gap);//if 1~3, read backwords
                string codon = seq.substr(pos, 3);
                if (is_in(codon, (isForword ? stopCodon : compStopCodon))){//reached next stop codon
                    int aaLength = processNum - stopPos;
                    if(aaLength<10000){
                        count[aaLength]++;
                    }else{
                        count[9999]++;
                    }
                    stopPos=processNum;
                }
            }
            int aaLength = processNum - stopPos;//update last pair if exists
            if(aaLength<10000){
                count[aaLength]++;
            }else{
                count[9999]++;
            }
        }
    }
}


void count_orf_stop(string &seq, vector<string> &startCodon, vector<string> &stopCodon, vector<int> &count){
	int seqLen = seq.length();

	vector<string> compStartCodon;
	for (int i = 0; i < startCodon.size(); i++){
		compStartCodon.push_back(reverse_complement(startCodon[i]));
	}
	vector<string> compStopCodon;
	for (int i = 0; i<stopCodon.size(); i++){
		compStopCodon.push_back(reverse_complement(stopCodon[i]));
	}


    int succession=0;
	for (int lane = 1; lane <= 6; lane++){
		bool isForword = true;
		if (lane > 3){
			isForword = false;
		}

		int gap = (lane - 1) % 3;
		int numOfLoop = (seqLen - gap) / 3;

		int stopPos = -1;
		int startPos = -1;


		//find fist stopcodon
		int processNum = 0;
		for (; processNum < numOfLoop; processNum++){
			int pos = (isForword ? (seqLen - 3) - (3 * processNum + gap) : 3 * processNum + gap);//if 1~3, read backwords
			string codon = seq.substr(pos, 3);
			if (is_in(codon, (isForword ? stopCodon : compStopCodon))){
				stopPos = processNum;
				break;
			}
		}

		if (stopPos != -1){//if found first stop codon
			startPos = stopPos;
			
			for (processNum++; processNum < numOfLoop; processNum++){
				int pos = (isForword ? (seqLen - 3) - (3 * processNum + gap) : 3 * processNum + gap);//if 1~3, read backwords
				string codon = seq.substr(pos, 3);
				if (is_in(codon, (isForword ? stopCodon : compStopCodon))){//reached next stop codon
					int aaLength = startPos - stopPos;
					if (aaLength > 0){//if startPos was updated, count
						count[aaLength]++;
					}else{
                        succession++;
                    }
					stopPos = startPos = processNum;//update stopPos and startPos
				}
				else if (is_in(codon, (isForword ? startCodon : compStartCodon))){//update startPos
					startPos = processNum;
				}
			}
			int aaLength = startPos - stopPos;//update last pair if exists
			if (aaLength > 0){
				count[aaLength]++;
			}
		}
	}
    cout<<"total succession : "<<succession<<endl;
}


void count_orf_start(string &seq, vector<string> &startCodon, vector<string> &stopCodon, vector<int> &count){
	int seqLen = seq.length();

	vector<string> compStartCodon;
	for (int i = 0; i < startCodon.size(); i++){
		compStartCodon.push_back(reverse_complement(startCodon[i]));
	}
	vector<string> compStopCodon;
	for (int i = 0; i<stopCodon.size(); i++){
		compStopCodon.push_back(reverse_complement(stopCodon[i]));
	}

	vector<int> v;//remember start codon position
	for (int lane = 1; lane <= 6; lane++){
		bool isForword = true;
		if (lane > 3){
			isForword = false;
		}

		int gap = (lane - 1) % 3;
		int numOfLoop = (seqLen - gap) / 3;

		for (int i = 0; i < numOfLoop; i++){
			int pos = (isForword ? 3 * i + gap : (seqLen - 3) - (3 * i + gap));//if 3~6, read backwords
			string codon = seq.substr(pos, 3);

			if (is_in(codon, (isForword ? startCodon : compStartCodon))){
				v.push_back(i);
			}
			else if (is_in(codon, (isForword ? stopCodon : compStopCodon))){
				for (int j = 0; j < v.size(); j++){
					int aaLength = i - v[j];//amino acid length(not including stop codon)
					count[aaLength]++;
				}
				v.clear();//found stopCodon corresponding
			}
		}
		v.clear();//assume linear genome so that some start codon has no stop codon corresponding
	}

}

string report_properties(string &seq, vector<string> &startCodon, vector<string> &stopCodon, int countMode){

	vector<int> count(10000, 0);

    if(countMode==0){
	    count_orf_start(seq, startCodon, stopCodon, count);
    }else if(countMode==1){
	    count_orf_stop(seq, startCodon, stopCodon, count);
    }else if (countMode==2){
	    count_orf_all(seq, stopCodon, count);
    }
        

	for (int i = 0; i < 10; i++){
		int sum = 0;
		for (int j = 0; j < 99; j++){
			sum += count[100 * i + j];
		}
		cout << "sum at " << 100 * (i+1) << " = " << sum << endl;
	}

    int cumsum=0;
	stringstream ss;
	for (int i = 0; i < 10000; i++){
		ss << count[i] << ",";
        cumsum+=count[i];
	}
    cout<<"TOTAL: "<<cumsum<<" orfs."<<endl;
    cout<<endl;

	string s = ss.str();
	s[s.length() - 1] = '\n';
	return s;
}


void test_count_orf(){
	vector<string> seqs(3);
	seqs[0] = "AATGAAAATGTAG";
	seqs[1] = reverse_complement(seqs[0]);
	seqs[2] = seqs[0] + "ATG";
	vector<string> startCodon = { "ATG" };
	vector<string> stopCodon = { "TAA", "TAG", "TGA" };
	vector<int> count(5);

	//***count_orf_start test***//

	//basic case
	string seq = "ATGAAAATGTAG";
	vector<int> answer{ 0, 1, 0, 1, 0 };
	fill(count.begin(), count.end(), 0);
	count_orf_start(seq, startCodon, stopCodon, count);
	assert(count == answer);

	//check lane 2
	seq = "AATGAAAATGTAG";
	fill(count.begin(), count.end(), 0);
	count_orf_start(seq, startCodon, stopCodon, count);
	assert(count == answer);

	//check reverse(lane 4)
	seq = reverse_complement(seq);
	fill(count.begin(), count.end(), 0);
	count_orf_start(seq, startCodon, stopCodon, count);
	assert(count == answer);

	//check start codon without corresponding stop
	seq = "AATGAAAATGTAGATG";
	fill(count.begin(), count.end(), 0);
	count_orf_start(seq, startCodon, stopCodon, count);
	assert(count == answer);
	


	//***count_orf_stop test***//
	//string without stop codon
	seq = "ATGAAAATG";
	fill(answer.begin(), answer.end(), 0);
	fill(count.begin(), count.end(), 0);
	count_orf_stop(seq, startCodon, stopCodon, count);
	assert(count == answer);

	//basic case
	seq = "ATGAAAATGTAG";
	answer[3] = 1;
	fill(count.begin(), count.end(), 0);
	count_orf_stop(seq, startCodon, stopCodon, count);
	assert(count == answer);

	//check lane 2
	seq = "AATGAAAATGTAG";
	fill(count.begin(), count.end(), 0);
	count_orf_stop(seq, startCodon, stopCodon, count);
	assert(count == answer);

	//check reverse(lane 4)
	seq = reverse_complement(seq);
	fill(count.begin(), count.end(), 0);
	count_orf_stop(seq, startCodon, stopCodon, count);
	assert(count == answer);

	//check stop codon without corresponding start
	seq = "TAAATGAAAATGTAG";
	fill(count.begin(), count.end(), 0);
	count_orf_stop(seq, startCodon, stopCodon, count);
	assert(count == answer);

	return;
}



int main(int argc, char *argv[]){
	string seqFilepath = argv[1];
	string cdsFilepath = argv[2];
	string outputFilepath = argv[3];

	string mode1 = argv[4];
	string mode2 = argv[5];
	pair<int, int> shuffleMode = make_pair(stoi(mode1), stoi(mode2));
    string mode=argv[6];
    int countMode=stoi(mode);


//	vector<string> startCodon = { "ATG" }; 
	vector<string> startCodon = { "ATG" ,"TTG", "CTG", "GTG"};
	vector<string> stopCodon = { "TAA", "TAG", "TGA" };
	string code_str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	GeneticCode gc(code_str);


	//read seq
	ifstream ifs(seqFilepath);
	if (!ifs){
		cerr << "Cannot open inputfile (" << seqFilepath << ")." << endl;
		return 1;
	}
	string buf;
	getline(ifs, buf);//skip header
	string seq;
	getline(ifs, seq);
	ifs.close();


    //output property
    cout<<"Processeing "<<seqFilepath<<"..."<<endl;
	cout << "seq length : " << seq.length() << " bp." << endl;
    cout<<"start codons : ";
    for(int i=0;i<startCodon.size();i++){
        cout<<startCodon[i]<<", ";
    }
    cout<<endl;
    cout<<"stop codons : ";
    for(int i=0;i<stopCodon.size();i++){
        cout<<stopCodon[i]<<", ";
    }
    cout<<endl;
    cout<<"shuffle mode : "<<shuffleMode.first<<", "<<shuffleMode.second<<endl;
    cout<<"count mode : "<<countMode<<endl;

	cout << endl;



	//create Regions;
	Regions cds(seq, cdsFilepath);
	cds.set_types(seq, startCodon, stopCodon);
	cds.show();
	

	//update codon usage count
	for (int i = 0; i < cds.regions.size(); i++){
		Region r = cds.regions[i];
		if (r.type == 0){
			string tmp = seq.substr(r.start + 3, r.end - r.start - 6);//exclude start codon and stop codon
			string subSeq = (r.isForword ? tmp : reverse_complement(tmp));
			gc.update_count(subSeq);
		}
	}
	gc.convert_to_freq();
	gc.show_count();
     
    gc.clear_count();
    shuffle_genome(seq, cds, gc, shuffleMode);
	//update codon usage count
	for (int i = 0; i < cds.regions.size(); i++){
		Region r = cds.regions[i];
		if (r.type == 0){
			string tmp = seq.substr(r.start + 3, r.end - r.start - 6);//exclude start codon and stop codon
			string subSeq = (r.isForword ? tmp : reverse_complement(tmp));
			gc.update_count(subSeq);
		}
	}
	gc.convert_to_freq();
	gc.show_count();
	
    
	

	//shuffle_genome(seq, cds, gc, make_pair(1, 2));
	//shuffle_genome(seq, cds, gc, make_pair(0, 2));
	

//	report_properties(seq, startCodon, stopCodon);
//	shuffle_genome(seq, cds, gc, mode);
//	report_properties(seq, startCodon, stopCodon);
	
	return 0;
}
