"""Replace the P4 block of the pool's job file with runs/sched/p4_jobs.txt (atomic rename; the job file before the swap is
kept as <jobs>.prev). The P4 block runs from its '# P4 confirmation' header to the next '# ' block header after the P4 ids.
Usage: python scripts/swap_p4_block.py [runs/sched/v_jobs.txt]"""
import os, sys
from pathlib import Path

J = Path(sys.argv[1] if len(sys.argv) > 1 else "runs/sched/v_jobs.txt"); L = J.read_text().splitlines(keepends=True)
a = next(i for i, x in enumerate(L) if x.startswith("# P4 confirmation"))
b = next((i for i in range(a + 1, len(L)) if L[i].startswith("# ") and not L[i].startswith(("# Tier 2",))), len(L))
new = L[:a] + Path("runs/sched/p4_jobs.txt").read_text().splitlines(keepends=True) + L[b:]
J.with_suffix(".txt.prev").write_text("".join(L)); tmp = J.with_suffix(".tmp"); tmp.write_text("".join(new)); os.replace(tmp, J)
print(f"{J}: P4 block lines {a + 1}-{b} ({b - a}) -> {len(new) - len(L) + (b - a)} lines; total {len(new)}")
