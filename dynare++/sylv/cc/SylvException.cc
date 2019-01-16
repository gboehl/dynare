/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvException.cpp,v 1.2 2004/10/01 10:30:40 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvException.hh"

#include <iostream>
#include <utility>

SylvException::SylvException(std::string f, int l)
  : file{std::move(f)}, line{l}
{
}

void
SylvException::printMessage() const
{
  std::cout << getMessage();
}

std::string
SylvException::getMessage() const
{
  return "From " + file + ':' + std::to_string(line) + '\n';
}

SylvExceptionMessage::SylvExceptionMessage(std::string f, int i,
                                           std::string mes)
  : SylvException{std::move(f), i}, message{std::move(mes)}
{
}

std::string
SylvExceptionMessage::getMessage() const
{
  return "At " + file + ':' + std::to_string(line) + ':' + message + '\n';
}
