/*
 * Copyright © 2004-2011 Ondra Kamenik
 * Copyright © 2019 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "SylvException.hh"

#include <iostream>
#include <utility>

SylvException::SylvException(std::string f, int l) : file {std::move(f)}, line {l}
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

SylvExceptionMessage::SylvExceptionMessage(std::string f, int i, std::string mes) :
    SylvException {std::move(f), i}, message {std::move(mes)}
{
}

std::string
SylvExceptionMessage::getMessage() const
{
  return "At " + file + ':' + std::to_string(line) + ':' + message + '\n';
}
