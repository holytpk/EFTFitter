// -*- C++ -*-
// author: afiq anuar
// short: non-plotting related ROOT utility methods

#ifndef ROOTUTIL_H
#define ROOTUTIL_H

#include "TSystem.h"
#include "TString.h"

// poor man's version of file finder by extension
std::vector<std::string> file_by_ext(const std::string &dir, const std::string &ext)
{
  // which really relies on ROOT's ability to run shell commands aha
  TString allfile = gSystem->GetFromPipe(("find " + dir + " -type f -name '*" + ext + "'").c_str());
  TString file;
  Ssiz_t index = 0;

  std::vector<std::string> v_file;
  while (allfile.Tokenize(file, index, "\n"))
    v_file.push_back(file.Data());

  return v_file;
}

#endif
