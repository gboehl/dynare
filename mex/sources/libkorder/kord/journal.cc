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

#include "journal.hh"
#include "kord_exception.hh"

#include <iomanip>
#include <cmath>
#include <ctime>
#include <limits>
#include <thread>

#ifndef _WIN32
# include <sys/time.h>     // For getrusage()
# include <sys/resource.h> // For getrusage()
# include <sys/utsname.h>  // For uname()
# include <cstdlib>        // For getloadavg()
# include <unistd.h>       // For sysconf()
# ifdef __APPLE__
#  include <sys/types.h>
#  include <sys/sysctl.h>
# endif
#else
# ifndef NOMINMAX
#  define NOMINMAX         // Do not define "min" and "max" macros
# endif
# include <windows.h>      // For GlobalMemoryStatus()
#endif

const std::chrono::time_point<std::chrono::high_resolution_clock> SystemResources::start = std::chrono::high_resolution_clock::now();

#ifndef _WIN32
long
SystemResources::pageSize()
{
  return sysconf(_SC_PAGESIZE);
}
#endif

long
SystemResources::availableMemory()
{
#if !defined(_WIN32) && !defined(__APPLE__)
  return sysconf(_SC_AVPHYS_PAGES)*pageSize();
#elif defined(__APPLE__)
  unsigned long usermem = 0;
  size_t len = sizeof usermem;
  static int mib[2] = { CTL_HW, HW_USERMEM };
  int retval = sysctl(mib, 2, &usermem, &len, NULL, 0);
  if (retval == 0)
    return static_cast<long>(usermem);
  return 0;
#else // _WIN32
  MEMORYSTATUS memstat;
  GlobalMemoryStatus(&memstat);
  return memstat.dwAvailPhys;
#endif
}

SystemResources::SystemResources()
{
  auto now = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = now - start;
  elapsed = duration.count();

#ifndef _WIN32
  struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);
  utime = rus.ru_utime.tv_sec+rus.ru_utime.tv_usec*1.0e-6;
  stime = rus.ru_stime.tv_sec+rus.ru_stime.tv_usec*1.0e-6;
  idrss = rus.ru_idrss;
  majflt = rus.ru_majflt * pageSize();
#else
  utime = std::numeric_limits<double>::quiet_NaN();
  stime = std::numeric_limits<double>::quiet_NaN();
  idrss = -1;
  majflt = -1;
#endif

#ifndef _WIN32
  getloadavg(&load_avg, 1);
#else
  load_avg = std::numeric_limits<double>::quiet_NaN();
#endif

  mem_avail = availableMemory();
}

void
SystemResources::diff(const SystemResources &pre)
{
  utime -= pre.utime;
  stime -= pre.stime;
  elapsed -= pre.elapsed;
  idrss -= pre.idrss;
  majflt -= pre.majflt;
}

// JournalRecord::operator<<() symmetry code
JournalRecord &
JournalRecord::operator<<(const IntSequence &s)
{
  operator<<('[');
  for (int i = 0; i < s.size(); i++)
    {
      operator<<(s[i]);
      if (i < s.size()-1)
        operator<<(',');
    }
  operator<<(']');
  return *this;
}

void
JournalRecord::writeFloatTabular(std::ostream &s, double d, int width)
{
  // Number of digits of integer part
  int intdigits = std::max(static_cast<int>(std::floor(log10(d))+1), 1);

  int prec = std::max(width - 1 - intdigits, 0);
  s << std::fixed << std::setw(width) << std::setprecision(prec) << d;
}

void
JournalRecord::writePrefix(const SystemResources &f)
{
  constexpr double mb = 1024*1024;
  std::ostringstream s;
  s << std::setfill('0');
  writeFloatTabular(s, f.elapsed, 7);
  s << "│" << recChar << std::setw(5) << ord << "│";
  writeFloatTabular(s, f.load_avg, 3);
  s << "│";
  writeFloatTabular(s, f.mem_avail/mb, 5);
  s << "│      │ ";
  for (int i = 0; i < 2*journal.getDepth(); i++)
    s << ' ';
  prefix = s.str();
}

void
JournalRecordPair::writePrefixForEnd(const SystemResources &f)
{
  constexpr double mb = 1024*1024;
  SystemResources difnow;
  difnow.diff(f);
  std::ostringstream s;
  s << std::setfill('0');
  writeFloatTabular(s, f.elapsed+difnow.elapsed, 7);
  s << "│E" << std::setw(5) << ord << "│";
  writeFloatTabular(s, difnow.load_avg, 3);
  s << "│";
  writeFloatTabular(s, difnow.mem_avail/mb, 5);
  s << "│";
  writeFloatTabular(s, difnow.majflt/mb, 6);
  s << "│ ";
  for (int i = 0; i < 2*journal.getDepth(); i++)
    s << ' ';
  prefix_end = s.str();
}

JournalRecordPair::~JournalRecordPair()
{
  journal.decrementDepth();
  writePrefixForEnd(flash);
  journal << prefix_end;
  journal << mes;
  journal << std::endl;
  journal.flush();
}

JournalRecord &
endrec(JournalRecord &rec)
{
  rec.journal << rec.prefix;
  rec.journal << rec.mes;
  rec.journal << std::endl;
  rec.journal.flush();
  rec.journal.incrementOrd();
  return rec;
}

void
Journal::printHeader()
{
  *this << "Dynare, version " << PACKAGE_VERSION << '\n'
        << '\n'
        << "System info: ";
#ifndef _WIN32
  utsname info;
  uname(&info);
  *this << info.sysname << " " << info.release << " " << info.version << " "
        << info.machine;
#else
  *this << "Windows";
#endif
  *this << ", processors online: " << std::thread::hardware_concurrency()
        << "\n\nStart time: ";
  std::time_t t = std::time(nullptr);
  // NB: in the date/time string, we avoid using locale-specific strings (#1751)
  *this << std::put_time(std::localtime(&t),
#ifndef _WIN32
                         "%F %T %z"
#else
                         // Specifiers introduced in C++11 don’t work under Windows…
                         "%Y-%m-%d %H:%M:%S"
#endif
                         )
        << "\n\n"
        << "  ┌────╼ elapsed time (seconds)                     \n"
        << "  │       ┌────╼ record unique identifier           \n"
        << "  │       │     ┌────╼ load average                 \n"
        << "  │       │     │    ┌────╼ available memory (MB)   \n"
        << "  │       │     │    │     ┌─────╼ major faults (MB)\n"
        << "  │       │     │    │     │                        \n"
        << "  ╽       ╽     ╽    ╽     ╽                        \n";
}
