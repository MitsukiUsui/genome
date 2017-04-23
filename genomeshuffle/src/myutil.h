#include <vector>
#include <string>
#include <sstream>

using namespace std;

//split string into vector of string
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

//check if the designated element is in the vector
template <class Content> bool is_in(Content c, vector<Content> &v){
	for (int i = 0; i < v.size(); i++){
		if (c == v[i]){
			return true;
		}
	}
	return false;
}

