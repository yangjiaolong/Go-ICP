#include "ConfigMap.hpp"

ConfigMap::ConfigMap() {  
};

ConfigMap::ConfigMap(const char * config_file) {
  
  std::ifstream fin(config_file);
  
  if (!fin.is_open()) {
    std::cout << "Unable to open config file '" << config_file << "'" << std::endl;
	exit(-2);
  }
  else {
  
    std::string buffer;
    
    while(std::getline(fin, buffer)) {
	  // Ignore comments
      if(buffer.c_str()[0] == '#') {
      }
      else {
        this->addLine(buffer);
      }
    }
    
    fin.close();
  }
};

ConfigMap::~ConfigMap() {

  for (std::list<float*>::iterator iter = allocated_memory_collector.begin();
       iter != allocated_memory_collector.end(); iter++) {
    delete [] (*iter);
  }
};

void ConfigMap::addLine(std::string line_) {

  // Swallow last character (carriage return: ASCII 13)
  if(line_.size() > 0)
  {
    if((int)line_.c_str()[line_.size() - 1] == 13)
    {
      line_.resize (line_.size () - 1);
    }
  }

  StringTokenizer st(line_, const_cast<char*>(" =;"));
  
  if (st.numberOfTokens() != 2) {
    return;
  }

  std::string key = st.nextToken();
  std::string val = st.nextToken();

  addPair(key, val);
};

void ConfigMap::addPair(std::string key_, std::string value_) {
  
  mappings[key_] = value_;
  
};

char * ConfigMap::get(char * key_) {

  std::string key(key_);

#ifdef VERBOSE
  std::cout << "DEBUG::ConfigMap.get()::key is `" << key_ << "'" << std::endl;

  std::string val = mappings[key];
  std::cout << "DEBUG::Requesting (" << key_ << ")->(" << val << ")" << std::endl;
  
  return const_cast<char*>(val.c_str());
#else
  return const_cast<char*>(mappings[key].c_str());
#endif  

};

char * ConfigMap::get(const char * key_) {
  return get(const_cast<char*>(key_));
};

int ConfigMap::getI(char * key_) {
  
  char * str_val = get(key_);
  
  if (str_val == NULL) {
    return 0;
  }
  else {
    return atoi(str_val);
  } 
};

int ConfigMap::getI(const char * key_) {
  return getI(const_cast<char*>(key_));
};

double ConfigMap::getF(char *key_) {
  
  char * str_val = get(key_);
  
  if (str_val == NULL) {
    return 0;
  } else {
    return atof(str_val);
  }
};

double ConfigMap::getF(const char * key_) {
  return getF(const_cast<char*>(key_));
};

float * ConfigMap::getVector(char * key_) {
  
  std::string key(key_);
  std::string val = mappings[key];
  
  if (val.empty()) {
    return NULL;
  }
  
  StringTokenizer st(val, const_cast<char*>("(,)"));
  
  float * return_vector = new float[st.numberOfTokens()];
  allocated_memory_collector.push_back(return_vector);

  for (int i = 0; st.numberOfTokens() > 0; i++) {
    return_vector[i] = (float)atof(st.nextToken().c_str());
  }
  
  return return_vector;
};

float * ConfigMap::getVector(const char * key_) {
  return getVector(const_cast<char*>(key_));
};

void ConfigMap::print() {
  
  for (std::map<std::string,std::string>::const_iterator iter = mappings.begin(); iter != mappings.end(); iter++)
  {
    std::cout << "(" << iter->first << ")->(" << iter->second << ")" << std::endl;
  }
};
