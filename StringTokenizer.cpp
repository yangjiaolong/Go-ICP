#include <iostream>
#include "StringTokenizer.hpp"

StringTokenizer::StringTokenizer() {
};

StringTokenizer::StringTokenizer(std::string str, char delim) {
  int delim_loc;
  
  while ((delim_loc = str.find_first_of(delim,0)) != std::string::npos) {

    if (str.substr(0, delim_loc).length() > 0) {
      tokens.push_back(str.substr(0,delim_loc));
    }
    str = str.substr(delim_loc+1, str.length());
  }

  if (str.length() > 0) {
    tokens.push_back(str);
  }
};

StringTokenizer::StringTokenizer(std::string str, char * delims) {
  int delim_loc;
  
  while ((delim_loc = str.find_first_of(delims,0)) != std::string::npos) {
    if (str.substr(0, delim_loc).length() > 0) {
      tokens.push_back(str.substr(0,delim_loc));
    }
    str = str.substr(delim_loc+1, str.length());
  }
  if (str.length() > 0) tokens.push_back(str);
};

StringTokenizer::~StringTokenizer() {
};

std::string StringTokenizer::nextToken() {
  if (!hasMoreTokens()) return "";
  std::string return_str(tokens.front());
  tokens.pop_front();
  return return_str;
};

bool StringTokenizer::hasMoreTokens() {
  return !tokens.empty();
};

int StringTokenizer::numberOfTokens() {
  return tokens.size();
};
