#!/bin/bash
# Memory guard for GPU jobs. Source it, then: guarded_run <cap_gb> <log> <cmd...>
#  - before start and every 5 s: swap > 6 GB -> politely quit heavy GUI apps (an app that asks to save is left open)
#  - every 5 s: job footprint (phys_footprint, includes MPS/unified memory) > cap -> kill the job (exit 137 recorded)
QUIT_APPS=("Safari" "Google Chrome" "Microsoft Word" "Microsoft Excel" "Microsoft PowerPoint" "Telegram" "Preview" "Mail" "Notes" "Music" "Spotify" "Slack" "zoom.us" "WhatsApp" "Discord")
swap_mb() { sysctl -n vm.swapusage | awk '{print $6}' | tr -d M | cut -d. -f1; }
free_apps() {
  for app in "${QUIT_APPS[@]}"; do
    if pgrep -xq "$app" || osascript -e "application \"$app\" is running" 2>/dev/null | grep -q true; then
      osascript -e 'with timeout of 15 seconds' -e "tell application \"$app\" to quit" -e 'end timeout' >/dev/null 2>&1 \
        && echo "[$(date +%T)] mem_guard: quit $app" || echo "[$(date +%T)] mem_guard: $app did not quit (prompt?) - left open"
    fi
  done
}
SWAP_LAST=0
check_swap() {  # act only while swap is above 6 GB AND growing: macOS swap stays high after pressure ends (audit 2026-10-30)
  local s=$(swap_mb); if [ "$s" -gt 6144 ] && [ "$s" -gt $((SWAP_LAST + 256)) ] && [ "$SWAP_LAST" -gt 0 ]; then echo "[$(date +%T)] mem_guard: swap ${s}M > 6G and growing -> closing apps"; free_apps; fi; SWAP_LAST=$s; }
descendants() { local c; for c in $(pgrep -P $1); do echo $c; descendants $c; done; }  # whole tree (blind.sh adds a level: zsh -> pixi -> python)
footprint_gb() {  # sum over the process tree rooted at $1
  local total=0
  for p in $1 $(descendants $1); do
    local f=$(footprint -p $p 2>/dev/null | sed -nE 's/.*Footprint: ([0-9.]+) (KB|MB|GB).*/\1 \2/p' | head -1 | awk '{if($2=="GB")print $1*1024; else if($2=="MB")print $1; else print $1/1024}')
    total=$(echo "$total + ${f:-0}" | bc)
  done; echo "scale=2; $total/1024" | bc
}
guarded_run() {
  local cap=$1 log=$2; shift 2; check_swap
  rm -f "$log.rc"; { "$@"; echo $? > "$log.rc"; } > "$log" 2>&1 & local pid=$! peak=0  # the command records its own exit code: zsh `wait` returned 145 (128 + SIGSTOP) for a job once paused by sched.sh, though it exited 0 (2026-11-04)
  while kill -0 $pid 2>/dev/null; do
    local fp=$(footprint_gb $pid); echo "$(date +%T) $fp" >> "$log.mem"; peak=$(echo "if ($fp > $peak) $fp else $peak" | bc)
    if [ $(echo "$fp > $cap" | bc) -eq 1 ]; then
      echo "[$(date +%T)] mem_guard: footprint ${fp}G > cap ${cap}G -> kill $pid" | tee -a "$log"; local tree="$(descendants $pid)"; kill -TERM $tree $pid 2>/dev/null; sleep 3; kill -KILL $tree $pid 2>/dev/null
    fi
    check_swap >> "$log"; sleep 5
  done
  wait $pid; local rc=$?; [ -s "$log.rc" ] && rc=$(cat "$log.rc"); rm -f "$log.rc"; echo "[$(date +%T)] mem_guard: exit $rc peak_footprint ${peak}G" | tee -a "$log"; return $rc
}
