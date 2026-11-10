#!/bin/zsh
# One pool job, started by sched.sh in its own session (setsid) so it survives a pool restart. Usage: sched_job.sh <id> <cap_gb> <cmd...>
cd ${0:A:h}/..; source scripts/mem_guard.sh; setopt NO_BG_NICE
S=runs/sched; id=$1; cap=$2; shift 2
guarded_run $cap $S/$id.log "$@"; rc=$?; [ $rc = 0 ] && touch $S/$id.done || { echo $rc > $S/$id.fail; }; rm -f $S/$id.run $S/$id.paused
