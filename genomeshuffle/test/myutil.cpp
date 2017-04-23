//test src/myutil.h

#include <iostream>
#include <string>
#include <vector>

#include "../src/myutil.h"
using namespace std;

int main(){
	//test split function
	string str="str,to,be,splited";
	vector<string> strSplited=split(str,',');
	for(string s:strSplited){
		cout<<s<<endl;
	}
	cout<<endl;

	//test is_in function
	vector<string> check{"to","test"};
	for (string s:check){
		if (is_in(s,strSplited)){
			cout<<s<<" is included."<<endl;
		}else{
			cout<<s<<" is not included."<<endl;
		}
	}
	return 0;
}
