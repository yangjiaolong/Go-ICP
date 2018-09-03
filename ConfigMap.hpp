
#ifndef CONFIGMAP_H
#define CONFIGMAP_H

#include <map>
#include <list>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "StringTokenizer.hpp"

class ConfigMap {

private:
  std::map<std::string,std::string> mappings;
  std::list<float*> allocated_memory_collector;

public:
  ConfigMap();
  ConfigMap(const char * config_file);
  ~ConfigMap();

  // adds a config line to the map
  //  (auto handles tokenising around '=' character)
  void addLine(std::string line_);

  // adds a key->value mapping
  void addPair(std::string key_, std::string value_);

  // gets the string value at a given key
  char * get(char * key_);
  char * get(const char * key_);

  // gets the integer value at a given key
  //  (uses the stdlib function 'atoi' to convert)
  int getI(char * key_);
  int getI(const char * key_);

  // gets the double value at a given key
  //  (uses the stdlib function 'atof' to convert)
  double getF(char * key_);
  double getF(const char * key_);

  // gets a vector of floats at the given key
  float * getVector(char * key_);
  float * getVector(const char * key_);

  // prints the contents of the map list to stdout
  void print();

};

#endif // for CONFIGMAP_H
