// Updated at Sept 30, 2022
// Updated by Mingqian Zhang
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// Code by: Mingqian Zhang, Xiaoquan Su
// Flex Meta-Storms version 1.0

#pragma once
#include<iostream>
#include<string>
#include<vector>
#include<sstream>
using namespace std;

string Check_FMS_Env() {
	
	if (getenv("FlexMetaStorms") == NULL) {

		cerr << "Error: Please set the environment variable \"FlexMetaStorms\" to the directory" << endl;
		exit(0);

	}
	string path;
	if (getenv("FlexMetaStorms") != NULL)
		path = getenv("FlexMetaStorms");
	return path;

	//debug
	//return "/opt/tools/parallel-meta/";
}

string Check_FMS_OTU(string otu) {

	string a_otu = otu;
	if (a_otu.size() > 4) {
		string prefix_4 = a_otu.substr(0, 4);
		if ((prefix_4 == "otu_") || (prefix_4 == "OTU_"))
			a_otu = a_otu.substr(4, a_otu.size() - 4);
	}
	return a_otu;
}

vector<string> split(const string& s, const string& seperator) {
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while (i != s.size()) {
		int flag = 0;
		while (i != s.size() && flag == 0) {
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x]) {
					++i;
					flag = 0;
					break;
				}
		}


		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0) {
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x]) {
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j) {
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}
template <class type>
type stringtonum(const string& str)
{
	istringstream iss(str);
	type num;
	iss >> num;
	return num;
}
