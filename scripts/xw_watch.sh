#!/bin/zsh
# Whole-slide arm stop rule (results/WHOLESLIDE_PLAN.md §4), with the investigators' system limits (2026-11-09):
# RAM 23 GB (every job's mem_guard cap) and swap > 2 GB for more than 15 min.
# - 2026-11-09 22:15 user order: a failing job (any slide, incl. 120) no longer stops the arm; slide 120 ran its
#   heaviest jobs, so a failure is treated as a bug: the job is set aside (runs/xw/SET_ASIDE.txt; <id>.done touched so
#   dependents and the gate proceed; xw_analyze skips missing outputs) and debugged separately.
# - Swap > 2 GB continuously for > 15 min, at any time -> the arm is declared failed and stopped (the running jobs are
#   named in runs/xw/ARM_FAILED).
# Failures of single jobs on the other slides do not stop the arm (reported as failed on that slide, §4).
# Exits when runs/xw/X_wholeslide.json exists (analysis done).
cd ${0:A:h}/..
kt() { local c; for c in $(pgrep -P $1); do kt $c; done; kill -TERM $1 2>/dev/null; kill -CONT $1 2>/dev/null; }
stop_arm() {
  { echo "[$(date +%T)] ARM FAILED: $1"; echo "running: ${(j: :)${(@)$(print -l runs/sched/xw_*.run(N))##*/}}"; } | tee runs/xw/ARM_FAILED
  launchctl bootout gui/$UID/com.prun.sfxw 2>/dev/null; p=$(cat runs/sched/sched.pid 2>/dev/null); [ -n "$p" ] && kill $p 2>/dev/null
  for r in runs/sched/xw_*.run(N); do kt $(awk '{print $1}' $r); mv $r ${r%.run}.stopped; done
  echo "[$(date +%T)] pool stopped"; exit 0
}
hi_since=0
while true; do
  [ -e runs/xw/X_wholeslide.json ] && { echo "[$(date +%T)] analysis written; watcher exits"; exit 0; }
  for f in runs/sched/xw_*.fail(N); do id=${${f:t}%.fail}; grep -qx "$id" runs/xw/SET_ASIDE.txt 2>/dev/null && continue
    echo $id >> runs/xw/SET_ASIDE.txt; touch runs/sched/$id.done
    echo "[$(date +%T)] set aside $id: $(tail -2 ${f%.fail}.log | tr '\n' ' ' | cut -c1-300)"; done
  s=$(sysctl -n vm.swapusage | awk '{print $6}' | tr -d M | cut -d. -f1); now=$(date +%s)
  if [ "$s" -gt 2048 ]; then [ $hi_since -eq 0 ] && { hi_since=$now; echo "[$(date +%T)] swap ${s}M > 2048M"; }; [ $((now - hi_since)) -gt 900 ] && stop_arm "swap > 2 GB for > 15 min (now ${s}M)"
  else [ $hi_since -ne 0 ] && echo "[$(date +%T)] swap back to ${s}M"; hi_since=0; fi
  sleep 30
done
