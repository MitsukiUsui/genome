#ifndef MYUTIL_H_
#define MYUTIL_H_

#include <vector>
#include <string>
#include <sstream>

//------------------------------------------------------------
//split string into vector of string according to separator
//------------------------------------------------------------
std::vector<std::string> split(std::string const &str, char sep)
{
    std::vector<std::string> v;
    std::stringstream ss(str);
    std::string buffer;
    while (std::getline(ss, buffer, sep)) {
        v.push_back(buffer);
    }
    return v;
}


//------------------------------------------------------------
//check if the designated element is in the vector
//------------------------------------------------------------
template <typename TContent> 
bool is_in(TContent c, std::vector<TContent> &v){
	for (int i = 0; i < v.size(); i++){
		if (c == v[i]){
			return true;
		}
	}
	return false;
}

#endif //MYUTIL_H_
