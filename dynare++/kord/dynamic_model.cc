// Copyright 2005, Ondra Kamenik

#include "dynamic_model.hh"

#include <iostream>
#include <algorithm>

void
NameList::print() const
{
  for (int i = 0; i < getNum(); i++)
    std::cout << getName(i) << '\n';
}

void
NameList::writeMat(mat_t *fd, const std::string &vname) const
{
  int maxlen = 0;
  for (int i = 0; i < getNum(); i++)
    maxlen = std::max(maxlen, static_cast<int>(getName(i).size()));

  if (maxlen == 0)
    return;

  auto m = std::make_unique<char[]>(getNum()*maxlen);

  for (int i = 0; i < getNum(); i++)
    for (int j = 0; j < maxlen; j++)
      if (j < static_cast<int>(getName(i).size()))
        m[j*getNum()+i] = getName(i)[j];
      else
        m[j*getNum()+i] = ' ';

  size_t dims[2];
  dims[0] = getNum();
  dims[1] = maxlen;

  matvar_t *v = Mat_VarCreate(vname.c_str(), MAT_C_CHAR, MAT_T_UINT8, 2, dims, m.get(), 0);

  Mat_VarWrite(fd, v, MAT_COMPRESSION_NONE);

  Mat_VarFree(v);
}

void
NameList::writeMatIndices(mat_t *fd, const std::string &prefix) const
{
  TwoDMatrix aux(1, 1);
  for (int i = 0; i < getNum(); i++)
    {
      aux.get(0, 0) = i+1;
      aux.writeMat(fd, prefix + "_i_" + getName(i));
    }
}
