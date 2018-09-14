#!/bin/bash

declare -i total=0
declare -i total_xfail=0
declare -i failed=0
declare -i xpassed=0
declare -a failed_tests=()
declare -a xpassed_tests=()

# Parse TRS Files
tosort=""
for file in $1 ; do
  # Find number of tests run in trs file
  ((total += $(grep ^:number-tests: "$file" | cut -d: -f3)))

  # Find number of tests failed in trs file
  numfailed=$(grep ^:number-failed-tests: "$file" | cut -d: -f3)
  if ((numfailed != 0)) ; then
      ((failed += numfailed))
      failed_tests+=($(grep ^:list-of-failed-tests: "$file" | cut -d: -f3))
  fi

  time=$(grep ^:elapsed-time: "$file" | cut -d: -f3)
  tosort="$tosort| $file - $time:"
done
((passed = total - failed))

# Parse XFAIL TRS Files
for file in $2 ; do
  # Find number of tests run in xfail trs file
  xfail=$(grep ^:number-tests: "$file" | cut -d: -f3)
  ((total_xfail += xfail))

  # Find number of tests failed in trs file
  numfailed=$(grep ^:number-failed-tests: "$file" | cut -d: -f3)
  if ((numfailed != xfail)) ; then
    ((xpassed += xfail - numfailed))
    xpassed_tests+=($(grep ^:list-of-passed-tests: "$file" | cut -d: -f3))
  fi

  time=$(grep ^:elapsed-time: "$file" | cut -d: -f3)
  tosort="$tosort| $file - $time:"
done
((xfailed = total_xfail - xpassed))
((total += total_xfail))

timing=$(echo "$tosort" | tr ":" "\n" | sort -rn -k4 | sed -e 's/$/:/' | head -n10)

# Determine if we are parsing Matlab or Octave trs files
if (($(grep -c '.m.trs' <<< "$1") == 0)); then
  prg=OCTAVE
  outfile=run_test_octave_output.txt
else
  prg=MATLAB
  outfile=run_test_matlab_output.txt
fi

# Print Output (to stdout and to a file)
{
    echo '================================'
    echo "DYNARE MAKE CHECK $prg RESULTS"
    echo '================================'
    echo "| TOTAL: $total"
    echo "|  PASS: $passed"
    echo "|  FAIL: $failed"
    echo "| XFAIL: $xfailed"
    echo "| XPASS: $xpassed"
    if ((failed > 0)) ; then
        echo '|'
        echo '| LIST OF FAILED TESTS:'
        for file in "${failed_tests[@]}" ; do
            echo "|     * $file"
        done
    fi
    if ((xpassed > 0)) ; then
        echo '|'
        echo '| LIST OF XPASSED TESTS:'
        for file in "${xpassed_tests[@]}" ; do
            echo "|     * $file"
        done
    fi
    echo '|'
    echo '| LIST OF 10 SLOWEST TESTS:'
    if [[ $prg == MATLAB ]]; then
        timing=$(sed 's/\.m\.trs/\.mod/g' <<< "$timing")
    else
        timing=$(sed 's/\.o\.trs/\.mod/g' <<< "$timing")
    fi
    echo "$timing" | tr ':' '\n' | sed -e 's/^[ \t]*//;/^$/d;s/^|[ ]/|     * /'
    echo
} | tee $outfile

# Exit with error code if some tests failed
((failed + xpassed == 0))
