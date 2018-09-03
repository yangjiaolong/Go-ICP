
#ifndef STRINGTOKENIZER_H
#define STRINGTOKENIZER_H

#include <list>
#include <string>

class StringTokenizer {

private:
  std::list<std::string> tokens;

public:
  StringTokenizer();
  StringTokenizer(std::string str, char delim = ' ');
  StringTokenizer(std::string str, char * delims);
  ~StringTokenizer();

  std::string nextToken();
  bool hasMoreTokens();
  int numberOfTokens();
};

#endif // for STRINGTOKENIZER_H
