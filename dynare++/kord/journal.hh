/*
 * Copyright © 2004 Ondra Kamenik
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

// Resource usage journal

#ifndef JOURNAL_H
#define JOURNAL_H

#include "int_sequence.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>

/* Implement static methods for accessing some system resources. An instance of
   this class is a photograph of these resources at the time of instantiation. */
struct SystemResources
{
  // The starting time of the executable
  static const std::chrono::time_point<std::chrono::high_resolution_clock> start;

#ifndef _WIN32
  static long pageSize();
#endif
  static long availableMemory();

  double load_avg;
  long mem_avail;
  double utime;
  double stime;
  double elapsed;
  long idrss;
  long majflt;

  SystemResources();
  void diff(const SystemResources &pre);
};

class Journal : public std::ofstream
{
  int ord;
  int depth;
public:
  explicit Journal(const std::string &fname)
    : std::ofstream(fname), ord(0), depth(0)
  {
    printHeader();
  }
  ~Journal() override
  {
    flush();
  }
  void printHeader();
  void
  incrementOrd()
  {
    ord++;
  }
  int
  getOrd() const
  {
    return ord;
  }
  void
  incrementDepth()
  {
    depth++;
  }
  void
  decrementDepth()
  {
    depth--;
  }
  int
  getDepth() const
  {
    return depth;
  }
};

class JournalRecord;
JournalRecord &endrec(JournalRecord &);

class JournalRecord
{
protected:
  char recChar;
  int ord;
public:
  Journal &journal;
  std::string prefix;
  std::string mes;
  SystemResources flash;
  using _Tfunc = JournalRecord &(*)(JournalRecord &);

  explicit JournalRecord(Journal &jr, char rc = 'M')
    : recChar(rc), ord(jr.getOrd()), journal(jr)
  {
    writePrefix(flash);
  }
  virtual ~JournalRecord() = default;
  JournalRecord &operator<<(const IntSequence &s);
  JournalRecord &
  operator<<(_Tfunc f)
  {
    (*f)(*this);
    return *this;
  }
  JournalRecord &
  operator<<(char c)
  {
    mes += c;
    return *this;
  }
  JournalRecord &
  operator<<(const std::string &s)
  {
    mes += s;
    return *this;
  }
  JournalRecord &
  operator<<(int i)
  {
    mes += std::to_string(i);
    return *this;
  }
  JournalRecord &
  operator<<(double d)
  {
    mes += std::to_string(d);
    return *this;
  }
protected:
  void writePrefix(const SystemResources &f);
  /* Writes a floating point number as a field of exactly ‘width’ characters
     large. Note that the width will not be respected if the integer part is
     too large. */
  static void writeFloatTabular(std::ostream &s, double d, int width);
};

/*
  Constructs a ‘pair’ of symmetric records with a RAII-like logic:
   — when fed with ‘endrec’, print the opening record
   — subsequent records will have depth increased by 1
   — when deleted, prints the symmetric closing record, and decrease depth
*/
class JournalRecordPair : public JournalRecord
{
  std::string prefix_end;
public:
  explicit JournalRecordPair(Journal &jr)
    : JournalRecord(jr, 'S')
  {
    journal.incrementDepth();
  }
  ~JournalRecordPair() override;
private:
  void writePrefixForEnd(const SystemResources &f);
};

#endif
