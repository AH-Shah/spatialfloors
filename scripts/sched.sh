#!/bin/zsh
# Resource-aware job pool (replaces one-at-a-time queues). Usage: nohup scripts/sched.sh <jobs.txt> > runs/sched/sched.log &
# jobs.txt lines: id|gpu or cpu|est_gb|cap_gb|done_file|after (comma ids or -)|command   (# = comment)
# Every 20 s: avail = free + inactive + speculative + purgeable pages; gpu% = IOAccelerator "Device Utilization %"
# (mean of the last 3 reads); load = 1-min load average. Launch the next ready job (>= 60 s after the last launch) when
#   memory: avail - (est not yet reached by running jobs) - est(job) >= RESERVE_GB and swap not growing (> 256 MB/round)
#   gpu job: gpu% < GPU_BUSY (underused) and running gpu jobs < GPU_MAX;  cpu job: load < NCPU - 1 and cpu jobs < CPU_MAX
# Remove: avail < PANIC_GB or swap > 6 GB and still growing (> 256 MB per round; macOS swap is sticky, 2026-10-30) -> kill the youngest job (requeued with est x 1.5); load > 2.5 x NCPU for 3 rounds -> pause (SIGSTOP) the youngest cpu job,
#   resumed once load < NCPU (2026-11-04; killing it lost hours of work). Status line every 5 min.
# Stop: kill $(cat runs/sched/sched.pid) - NOT pkill -f. Jobs run in their own session (sched_job.sh via setsid), so they survive
# a restart, even a launchctl bootout (their .run pids are re-adopted; jobs running outside the pool are detected by pgrep).
# Every job runs under mem_guard (footprint cap). State: runs/sched/<id>.{run,done,fail}. No output unless something happens.
cd ${0:A:h}/..; source scripts/mem_guard.sh
setopt NO_BG_NICE  # 2026-11-04: zsh gave each & subshell +5 (jobs ran at nice 10), which pinned every job to the E cores
JOBS=$1; S=runs/sched; mkdir -p $S; echo $$ > $S/sched.pid
RESERVE_GB=${RESERVE_GB:-3}; PANIC_GB=${PANIC_GB:-1.5}; GPU_MAX=${GPU_MAX:-4}; CPU_MAX=${CPU_MAX:-8}; GPU_BUSY=${GPU_BUSY:-85}; NCPU=$(sysctl -n hw.ncpu)
gpu_now() { ioreg -r -c IOAccelerator -d 1 2>/dev/null | grep -o '"Device Utilization %"=[0-9]*' | head -1 | cut -d= -f2; }
load1() { sysctl -n vm.loadavg | awk '{print $2}'; }
G1=0; G2=0; G3=0; hot=0; tick=0
PG=16384
avail_gb() { vm_stat | awk -v pg=$PG '/Pages (free|inactive|speculative|purgeable)/{gsub("\\.","",$NF); s+=$NF} END{printf "%.2f", s*pg/2^30}'; }
swap0=$(swap_mb); swap_prev=$swap0; last=0
typeset -A EST
log() { echo "[$(date +%T)] $*" >> $S/sched.log; }
killtree() { local c; for c in $(pgrep -P $1); do killtree $c; done; kill -TERM $1 2>/dev/null; kill -CONT $1 2>/dev/null; }
sigtree() { local c; for c in $(pgrep -P $2); do sigtree $1 $c; done; kill -$1 $2 2>/dev/null; }
log "start reserve=${RESERVE_GB} panic=${PANIC_GB} gpu_max=${GPU_MAX} gpu_busy=${GPU_BUSY}% cpu_max=${CPU_MAX} ncpu=${NCPU} swap0=${swap0}M avail=$(avail_gb)G"
while true; do
  pending=0; ngpu=0; ncpu=0; growth=0; youngest=""; yt=0; ycpu=""; yct=0
  for f in $S/*.run(N); do  # running jobs: count, pending growth, youngest
    id=${f:t:r}; read pid cls est t0 < $f
    if ! kill -0 $pid 2>/dev/null; then rm -f $f $S/$id.paused; continue; fi
    if [ $cls = gpu ]; then ngpu=$((ngpu + 1)); else ncpu=$((ncpu + 1)); [ $t0 -gt $yct ] && [ ! -e $S/$id.paused ] && { yct=$t0; ycpu=$id; ycpid=$pid; }; fi
    fp=$(tail -1 $S/$id.log.mem 2>/dev/null | awk '{print $2}'); g=$(echo "$est - ${fp:-0}" | bc); [ $(echo "$g > 0" | bc) -eq 1 ] && growth=$(echo "$growth + $g" | bc)
    [ $t0 -gt $yt ] && { yt=$t0; youngest=$id; ypid=$pid; }
  done
  av=$(avail_gb); sw=$(swap_mb); G1=$G2; G2=$G3; G3=$(gpu_now); G3=${G3:-0}; gp=$(( (G1 + G2 + G3) / 3 )); ld=$(load1)
  tick=$((tick + 1)); [ $((tick % 15)) = 0 ] && log "status avail=${av}G swap=${sw}M load=${ld} gpu=${gp}% running gpu=$ngpu cpu=$ncpu"
  if [ $(echo "$ld > 2.5 * $NCPU" | bc) -eq 1 ]; then hot=$((hot + 1)); else hot=0; fi
  # overload: pause (SIGSTOP) the youngest running cpu job, resume paused jobs once load < NCPU; no work is lost
  # (2026-11-04: killing wasted 6 half-done NCEM fits in 12 h; NCEM's TF thread pools spike load to 40+ mid-run)
  if [ $hot -ge 3 ] && [ -n "$ycpu" ]; then
    log "CPU overload load=${ld} -> pause $ycpu"; sigtree STOP $ycpid; touch $S/$ycpu.paused; hot=0
  elif [ $(echo "$ld < $NCPU" | bc) -eq 1 ]; then
    for f in $S/*.paused(N); do id=${f:t:r}; read pid x 2>/dev/null < $S/$id.run && sigtree CONT $pid; rm -f $f; log "resume $id load=${ld}"; done
  fi
  if [ -n "$youngest" ] && { [ $(echo "$av < $PANIC_GB" | bc) -eq 1 ] || { [ $sw -gt 6144 ] && [ $sw -gt $((swap_prev + 256)) ]; }; }; then
    log "PANIC avail=${av}G swap=${sw}M -> kill $youngest (requeued)"; killtree $ypid
    sleep 30; rm -f $S/$youngest.run $S/$youngest.fail; touch $S/$youngest.requeue; continue
  fi
  now=$(date +%s); launched=0; gpuwait=0
  for pass in gpu cpu; do  # gpu jobs first; 2026-11-05: 8 GB Cellina jobs starved for hours while 3-6 GB cpu jobs took every freed slot
  while IFS='|' read -r id cls est cap done after cmd; do
    [[ -z $id || $id == \#* || $cls != $pass ]] && continue
    [ -e $S/$id.done ] || [ -e $S/$id.fail ] || [ -e $S/$id.run ] && continue
    if [ -e "$done" ]; then touch $S/$id.done; continue; fi
    pgrep -qf "${cmd##*python }" && continue  # already running outside the pool
    ((pending++))
    if [ "$after" != "-" ]; then ok=1; for d in ${(s:,:)after}; do [ -e $S/$d.done ] || ok=0; done; [ $ok = 0 ] && continue; fi
    [ $launched = 1 ] || [ $((now - last)) -lt 60 ] && continue
    [ -e $S/$id.requeue ] && est=$(echo "$est * 1.5" | bc)
    if [ $cls = gpu ]; then [ $ngpu -ge $GPU_MAX ] || [ $gp -ge $GPU_BUSY ] && continue; fi
    [ $cls = cpu ] && [ $gpuwait = 1 ] && continue  # a ready gpu job waits for memory: let running cpu jobs drain
    if [ $cls = cpu ]; then [ $ncpu -ge $CPU_MAX ] || [ $(echo "$ld >= $NCPU - 1" | bc) -eq 1 ] && continue; fi
    [ $sw -gt $((swap_prev + 256)) ] && continue
    [ $(echo "$av - $growth - $est >= $RESERVE_GB" | bc) -eq 1 ] || { [ $cls = gpu ] && [ $ncpu -gt 0 ] && gpuwait=1; continue; }
    python3 -c 'import os,sys; os.setsid(); os.execvp(sys.argv[1], sys.argv[1:])' zsh scripts/sched_job.sh $id $cap ${(z)cmd} &
    echo "$! $cls $est $now" > $S/$id.run; last=$now; launched=1
    log "launch $id ($cls est ${est}G) avail=${av}G growth=${growth}G gpu=$ngpu/${gp}% cpu=$ncpu/load ${ld} swap=${sw}M"
  done < $JOBS
  done
  swap_prev=$sw
  if [ $pending = 0 ] && [ -z "$(ls $S/*.run 2>/dev/null)" ]; then log "all done"; break; fi
  sleep 20
done
