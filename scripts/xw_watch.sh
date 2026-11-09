#!/bin/zsh
# Whole-slide arm stop rule (results/WHOLESLIDE_PLAN.md §4): if any slide-120 job fails (mem_guard cap kill or crash),
# the arm is declared failed and stopped: pool booted out, running jobs killed, runs/xw/ARM_FAILED written.
# Exits quietly once runs/xw/GATE120.ok exists (slide 120 complete).
cd /Users/Hemato/project_3/antig/spatialfloors
kt() { local c; for c in $(pgrep -P $1); do kt $c; done; kill -TERM $1 2>/dev/null; kill -CONT $1 2>/dev/null; }
while true; do
  [ -e runs/xw/GATE120.ok ] && { echo "[$(date +%T)] slide 120 complete; watcher exits"; exit 0; }
  f=$(ls runs/sched/xw_*_120.fail 2>/dev/null | head -1)
  if [ -n "$f" ]; then
    echo "[$(date +%T)] ARM FAILED: $f" | tee runs/xw/ARM_FAILED; tail -5 ${f%.fail}.log >> runs/xw/ARM_FAILED
    launchctl bootout gui/$UID/com.prun.sfxw 2>/dev/null; p=$(cat runs/sched/sched.pid 2>/dev/null); [ -n "$p" ] && kill $p 2>/dev/null
    for r in runs/sched/xw_*.run(N); do kt $(head -1 $r | awk '{print $1}'); mv $r ${r%.run}.stopped; done
    echo "[$(date +%T)] pool stopped"; exit 0
  fi
  sleep 30
done
