// Copyright 2005, Ondra Kamenik

#include "dynamic_model.hh"

void
NameList::print() const
{
  for (int i = 0; i < getNum(); i++)
    printf("%s\n", getName(i));
}

void
NameList::writeMat(mat_t *fd, const char *vname) const
{
  int maxlen = 0;
  for (int i = 0; i < getNum(); i++)
    if (maxlen < (int) strlen(getName(i)))
      maxlen = (int) strlen(getName(i));

  if (maxlen == 0)
    return;

  auto *m = new char[getNum()*maxlen];

  for (int i = 0; i < getNum(); i++)
    for (int j = 0; j < maxlen; j++)
      if (j < (int) strlen(getName(i)))
        m[j*getNum()+i] = getName(i)[j];
      else
        m[j*getNum()+i] = ' ';

#if MATIO_MAJOR_VERSION > 1 || (MATIO_MAJOR_VERSION == 1 && MATIO_MINOR_VERSION >= 5)
  size_t dims[2];
  const matio_compression compression = MAT_COMPRESSION_NONE;
#else
  int dims[2];
  const int compression = COMPRESSION_NONE;
#endif
  dims[0] = getNum();
  dims[1] = maxlen;

  matvar_t *v = Mat_VarCreate(vname, MAT_C_CHAR, MAT_T_UINT8, 2, dims, m, 0);

  Mat_VarWrite(fd, v, compression);

  Mat_VarFree(v);
  delete[] m;
}

void
NameList::writeMatIndices(mat_t *fd, const char *prefix) const
{
  char tmp[100];
  TwoDMatrix aux(1, 1);
  for (int i = 0; i < getNum(); i++)
    {
      sprintf(tmp, "%s_i_%s", prefix, getName(i));
      aux.get(0, 0) = i+1;
      aux.writeMat(fd, tmp);
    }
}
