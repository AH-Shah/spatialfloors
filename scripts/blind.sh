#!/bin/zsh
# P4 blind wrapper (PLAN §17): run a job with all output sealed in runs/p4/sealed/<id>.log (read only after p4-final);
# print only the exit code, and on failure the Python traceback (from the last "Traceback" line), so crashes can be fixed
# without seeing scores. Usage: scripts/blind.sh <id> <command ...>
id=$1; shift; mkdir -p runs/p4/sealed; L=runs/p4/sealed/$id.log
"$@" > $L 2>&1; rc=$?
echo "[$id] exit $rc"
if [ $rc -ne 0 ]; then awk '/^Traceback/{buf=""; on=1} on{buf=buf $0 "\n"} END{printf "%s", buf}' $L | tail -40; fi
exit $rc
