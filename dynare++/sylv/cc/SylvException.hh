/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvException.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLV_EXCEPTION_H
#define SYLV_EXCEPTION_H

#include <string>

#include "SylvMemory.hh"

class SylvException : public MallocAllocator
{
protected:
  std::string file;
  int line;
public:
  SylvException(std::string f, int l);
  virtual ~SylvException() = default;
  void printMessage() const;
  virtual std::string getMessage() const;
};

class SylvExceptionMessage : public SylvException
{
  std::string message;
public:
  SylvExceptionMessage(std::string f, int l, std::string mes);
  std::string getMessage() const override;
};

// define macros:
#define SYLV_MES_EXCEPTION(mes) (SylvExceptionMessage(__FILE__, __LINE__, mes))

#endif /* SYLV_EXCEPTION_H */
