/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvException.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLV_EXCEPTION_H
#define SYLV_EXCEPTION_H

#include "SylvMemory.hh"

class SylvException : public MallocAllocator
{
protected:
  char file[50];
  int line;
public:
  SylvException(const char *f, int l);
  virtual ~SylvException() = default;
  virtual int printMessage(char *str, int maxlen) const;
  void printMessage() const;
};

class SylvExceptionMessage : public SylvException
{
  char message[500];
public:
  SylvExceptionMessage(const char *f, int l, const char *mes);
  int printMessage(char *str, int maxlen) const override;
};

// define macros:
#define SYLV_MES_EXCEPTION(mes) (SylvExceptionMessage(__FILE__, __LINE__, mes))

#endif /* SYLV_EXCEPTION_H */
