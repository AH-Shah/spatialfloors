#!/bin/zsh
# Resumable downloads for an unstable network. Manifest runs/dl/manifest.tsv: <url> <local path> [md5]
#   zsh scripts/dl.sh start    background worker: each file resumes (curl -C -) until its size matches the server, then md5 check
#   zsh scripts/dl.sh status   per file: done / partial (% of server size) / missing; worker running or not
#   zsh scripts/dl.sh pause    stop the worker and its curl (partial files are kept)
#   zsh scripts/dl.sh resume   = start (continues where the partial files stopped)
# Log: runs/dl/dl.log. Retries forever with a 10 s pause; a failed md5 renames the file to .bad and starts it again.
cd ${0:a:h}/..; D=runs/dl; M=$D/manifest.tsv; mkdir -p $D
remote_size() { curl -sIL --http1.1 --max-time 30 "$1" | awk 'tolower($1)=="content-length:"{s=$2} END{gsub("\r","",s); print s}'; }
worker() {
  echo $$ > $D/worker.pid
  while IFS=$'\t' read -r url out md5; do
    [[ -z $url || $url == \#* ]] && continue; mkdir -p ${out:h}
    while true; do
      want=$(remote_size "$url"); have=$(stat -f %z "$out" 2>/dev/null || echo 0)
      if [ -n "$want" ] && [ "$have" = "$want" ]; then
        if [ -n "$md5" ] && [ "$(md5 -q "$out")" != "$md5" ]; then echo "[$(date +%T)] md5 mismatch $out -> .bad, restart" >> $D/dl.log; mv "$out" "$out.bad"; continue; fi
        echo "[$(date +%T)] done $out" >> $D/dl.log; break
      fi
      curl -sS --http1.1 -L -C - --retry 5 --retry-delay 5 --speed-limit 1000 --speed-time 60 -o "$out" "$url" 2>> $D/dl.log || sleep 10
    done
  done < $M
  echo "[$(date +%T)] ALL_DONE" >> $D/dl.log; rm -f $D/worker.pid
}
case $1 in
  start|resume) [ -f $D/worker.pid ] && kill -0 $(cat $D/worker.pid) 2>/dev/null && { echo "already running"; exit 0; }
                nohup zsh ${0:a} _worker > /dev/null 2>&1 &; echo "started";;
  _worker) worker;;
  pause) [ -f $D/worker.pid ] && { p=$(cat $D/worker.pid); pkill -P $p; kill $p; rm -f $D/worker.pid; echo "paused"; } || echo "not running";;
  status) while IFS=$'\t' read -r url out md5; do [[ -z $url || $url == \#* ]] && continue; h=$(stat -f %z "$out" 2>/dev/null || echo 0); w=$(remote_size "$url")
            printf "%-60s %s\n" ${out:t} "$( [ "$h" = "$w" ] && echo done || echo "partial $((h * 100 / (w > 0 ? w : 1)))% ($h / $w)")"; done < $M
          [ -f $D/worker.pid ] && kill -0 $(cat $D/worker.pid) 2>/dev/null && echo "worker running (pid $(cat $D/worker.pid))" || echo "worker not running"; tail -2 $D/dl.log 2>/dev/null;;
  *) sed -n 2,7p $0;;
esac
