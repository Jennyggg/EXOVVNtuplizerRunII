#include <iostream>
#include <sstream>
#include <string>
using namespace std;
void splitstring(string origin, vector<string>& split ){
  istringstream iss(origin);
  string word;
  while (iss >> word){
    cout<<word<<"\n";
    split.push_back(word);
  }
  cout<<endl;
}
