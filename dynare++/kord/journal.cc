// Copyright (C) 2004-2011, Ondra Kamenik

#include "journal.hh"
#include "kord_exception.hh"

#include <iomanip>
#include <cmath>
#include <ctime>

#ifndef __MINGW32__
# include <sys/time.h>     // For getrusage()
# include <sys/resource.h> // For getrusage()
# include <sys/utsname.h>  // For uname()
# include <cstdlib>        // For getloadavg()
# include <unistd.h>       // For sysconf()
#else
# ifndef NOMINMAX
#  define NOMINMAX         // Do not define "min" and "max" macros
# endif
# include <windows.h>      // For GlobalMemoryStatus()
#endif

const std::chrono::time_point<std::chrono::high_resolution_clock> SystemResources::start = std::chrono::high_resolution_clock::now();

/* The pagesize is set to 1024 bytes on Windows. Real pagesize can differ but
   it is not important. We can do this since Windows kernel32
   GlobalMemoryStatus() call returns a number of bytes. */
long
SystemResources::pageSize()
{
#ifndef __MINGW32__
  return sysconf(_SC_PAGESIZE);
#else
  return 1024;
#endif
}

long
SystemResources::physicalPages()
{
#ifndef __MINGW32__
  return sysconf(_SC_PHYS_PAGES);
#else
  MEMORYSTATUS memstat;
  GlobalMemoryStatus(&memstat);
  return memstat.dwTotalPhys/1024;
#endif
}

long
SystemResources::availablePhysicalPages()
{
#ifndef __MINGW32__
  return sysconf(_SC_AVPHYS_PAGES);
#else
  MEMORYSTATUS memstat;
  GlobalMemoryStatus(&memstat);
  return memstat.dwAvailPhys/1024;
#endif
}

long
SystemResources::onlineProcessors()
{
#ifndef __MINGW32__
  return sysconf(_SC_NPROCESSORS_ONLN);
#else
  return -1;
#endif
}

long
SystemResources::availableMemory()
{
  return pageSize()*availablePhysicalPages();
}

SystemResources::SystemResources()
{
  auto now = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = now - start;
  elapsed = duration.count();

#ifndef __MINGW32__
  struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);
  utime = rus.ru_utime.tv_sec+rus.ru_utime.tv_usec*1.0e-6;
  stime = rus.ru_stime.tv_sec+rus.ru_stime.tv_usec*1.0e-6;
  idrss = rus.ru_idrss;
  majflt = rus.ru_majflt;
#else
  utime = -1.0;
  stime = -1.0;
  idrss = -1;
  majflt = -1;
#endif

#ifndef __MINGW32__
  getloadavg(&load_avg, 1);
#else
  load_avg = -1.0;
#endif

  pg_avail = availablePhysicalPages();
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

// |JournalRecord::operator<<| symmetry code
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
  s  << ':' << recChar << std::setw(5) << ord << ':';
  writeFloatTabular(s, f.load_avg, 3);
  s << ':';
  writeFloatTabular(s, f.pg_avail*SystemResources::pageSize()/mb, 5);
  s << ":      : ";
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
  s << ":E" << std::setw(5) << ord << ':';
  writeFloatTabular(s, difnow.load_avg, 3);
  s << ':';
  writeFloatTabular(s, difnow.pg_avail*SystemResources::pageSize()/mb, 5);
  s << ':';
  writeFloatTabular(s, difnow.majflt*SystemResources::pageSize()/mb, 6);
  s << ": ";
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
  *this << "This is Dynare++, Copyright (C) 2004-2011, Ondra Kamenik\n"
        << "Dynare++ comes with ABSOLUTELY NO WARRANTY and is distributed under\n"
        << "GPL: modules integ, tl, kord, sylv, src, extern and documentation\n"
        << "LGPL: modules parser, utils\n"
        << " for GPL  see http://www.gnu.org/licenses/gpl.html\n"
        << " for LGPL see http://www.gnu.org/licenses/lgpl.html\n"
        << "\n\n"
        << "System info: ";
#ifndef __MINGW32__
  utsname info;
  uname(&info);
  *this << info.sysname << " " << info.release << " " << info.version << " "
        << info.machine << ", processors online: " << SystemResources::onlineProcessors();
#else
  *this << "(not implemented for MinGW)";
#endif
  *this << "\n\nStart time: ";
  std::time_t t = std::time(nullptr);
  *this << std::put_time(std::localtime(&t), "%c %Z")
        << "\n\n"
        << "  ------ elapsed time (seconds)                     \n"
        << "  |       ------ record unique identifier           \n"
        << "  |       |     ------ load average                 \n"
        << "  |       |     |    ------ available memory (MB)   \n"
        << "  |       |     |    |     ------  major faults (MB)\n"
        << "  |       |     |    |     |                        \n"
        << "  V       V     V    V     V                        \n"
        << "\n";
}
