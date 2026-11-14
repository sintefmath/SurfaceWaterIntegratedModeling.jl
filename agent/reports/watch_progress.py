#!/usr/bin/env python3
"""Watch fill_sequence's verbose progress lines in a log and kill the run when it stops making
real progress.  Usage: watch_progress.py LOG PID [OUT]

Progress lines look like
    Fill sequence iteration: 120, t = 3.51, changed: [(1234, true)]
and appear every 10 iterations.  The run is killed (SIGKILL) as soon as one of these holds:
- crawl:  over the last 2000 iterations simulated time advanced by less than 0.01;
- flip:   a trap appears both filling and unfilling, at least 3 times each, in the last 2000
          iterations (sampled every 10);
- repeat: one trap is what changed in every sampled line of the last 2000 iterations (an
          alternating pair sampled every 10th event always shows the same kind);
- stall:  no new progress line for 15 minutes (one solve that long);
- runaway: more than 100,000 iterations.
The reason and the last progress lines are written to OUT (default LOG + ".watch").
"""
import collections, os, re, signal, sys, time

LINE = re.compile(r"Fill sequence iteration: (\d+), t = ([-\d.eE+]+), changed: \[(.*)\]")
PAIR = re.compile(r"\((\d+), (true|false|1|0)\)")
WINDOW, CRAWL_DT, FLIPS, STALL_S, MAX_IT = 2000, 0.01, 3, 15 * 60, 100_000


def main(log, pid, out):
    seen = collections.deque()           # (iteration, t, [(trap, filled)]) within the window
    pos, last_new = 0, time.time()

    def stop(reason):
        os.kill(pid, signal.SIGKILL)
        with open(out, "w") as f:
            f.write(reason + "\n")
            for it, t, ch in list(seen)[-20:]:
                f.write(f"{it} t={t} {ch}\n")
        sys.exit(0)

    while True:
        try:
            os.kill(pid, 0)
        except OSError:
            with open(out, "w") as f:
                f.write("run ended by itself\n")
            return
        with open(log, errors="replace") as f:
            f.seek(pos)
            new = f.read()
            pos = f.tell()
        for m in LINE.finditer(new):
            changed = [(int(i), v in ("true", "1")) for i, v in PAIR.findall(m.group(3))]
            seen.append((int(m.group(1)), float(m.group(2)), changed))
            last_new = time.time()
        while seen and seen[-1][0] - seen[0][0] > WINDOW:
            seen.popleft()
        if seen:
            it, t = seen[-1][0], seen[-1][1]
            if it > MAX_IT:
                stop(f"runaway: {it} iterations")
            if seen[-1][0] - seen[0][0] >= WINDOW - 10 and t - seen[0][1] < CRAWL_DT:
                stop(f"crawl: t {seen[0][1]} -> {t} over iterations {seen[0][0]}-{it}")
            counts = collections.Counter((tr, v) for _, _, ch in seen for tr, v in ch)
            flippers = sorted({tr for tr, v in counts
                               if counts[(tr, True)] >= FLIPS and counts[(tr, False)] >= FLIPS})
            if flippers:
                stop(f"flip: traps {flippers} fill and unfill repeatedly (t = {t}, iteration {it})")
            traps = {tr for _, _, ch in seen for tr, _ in ch}
            if len(seen) >= WINDOW // 10 - 1 and len(traps) == 1:
                stop(f"repeat: only trap {traps.pop()} changes, iterations {seen[0][0]}-{it}, "
                     f"t {seen[0][1]} -> {t}")
        if time.time() - last_new > STALL_S:
            stop("stall: no progress line for 15 minutes")
        time.sleep(20)


if __name__ == "__main__":
    log, pid = sys.argv[1], int(sys.argv[2])
    main(log, pid, sys.argv[3] if len(sys.argv) > 3 else log + ".watch")
