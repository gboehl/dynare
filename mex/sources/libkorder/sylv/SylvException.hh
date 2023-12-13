/*
 * Copyright © 2004-2011 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

#ifndef SYLV_EXCEPTION_HH
#define SYLV_EXCEPTION_HH

#include <string>

class SylvException
{
protected:
  std::string file;
  int line;

public:
  SylvException(std::string f, int l);
  virtual ~SylvException() = default;
  void printMessage() const;
  [[nodiscard]] virtual std::string getMessage() const;
};

class SylvExceptionMessage : public SylvException
{
  std::string message;

public:
  SylvExceptionMessage(std::string f, int l, std::string mes);
  [[nodiscard]] std::string getMessage() const override;
};

// define macros:
#define SYLV_MES_EXCEPTION(mes) (SylvExceptionMessage(__FILE__, __LINE__, mes))

#endif
