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


//--------------------------------------------------------------------------------
//Easy log class
//--------------------------------------------------------------------------------
class Logger{
    std::string filepath;
    std::ofstream ofs;
    std::stringstream ss;

public:
    Logger(std::string const & filepath){
        this->filepath=filepath;
        this->ofs.open(filepath);
        if (!ofs) {
            std::cerr << "ERROR: Could not open " << filepath << std::endl;
        }
    }

    ~Logger() {
        ofs.close();
    }

    template<typename T>
    Logger& operator<< (const T& data)
    {
        std::string type = typeid(data).name();
        if (type == typeid(std::string).name() || type[0]=='A'){ //if string, quote with ""
            ss << "\"" <<data << "\",";
        } else {
            ss << data << ",";
        }
        return *this;
    }

    void flush() {
        std::string message = ss.str();
        ss.str(std::string()); //clear
        if (message.length() > 0) {
            message.pop_back();
        }
        ofs<<message<<std::endl;
    }

    void __show() {
        std::cout << ss.str()<<std::endl;
    }
};

#endif //MYUTIL_H_
