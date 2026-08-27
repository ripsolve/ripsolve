#!/usr/bin/env python3
"""Render the tractable-set results as a page.

Generated from the measurements rather than transcribed, so the page cannot drift from
`bench/out/tractable_results.json`.

Usage:  bench/tractable_page.py <out.html>
"""

import html
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = json.loads((ROOT / "bench" / "out" / "tractable_results.json").read_text())
DEST = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "bench" / "out" / "tractable.html"

SOLVERS = ["ripsolve", "HiGHS", "SCIP", "CBC", "commercial"]
LABEL = {"ripsolve": "ripsolve", "HiGHS": "HiGHS", "SCIP": "SCIP",
         "CBC": "CBC", "commercial": "commercial"}
LIMIT = DATA["limit"]


def width(seconds, solved):
    """Bar length. Log-scaled: solve times span 0.05s to the full limit."""
    if not solved:
        return 100.0
    capped = min(max(seconds, 0.05), LIMIT)
    return 4.0 + 96.0 * math.log(capped / 0.05) / math.log(LIMIT / 0.05)


def main():
    rows = DATA["rows"]
    tally = DATA["tally"]
    total = len(rows)

    summary = "\n".join(
        f'''        <div class="tally-row{' is-subject' if s == 'ripsolve' else ''}">
          <span class="tally-name">{html.escape(LABEL[s])}</span>
          <span class="tally-track">
            <span class="tally-fill" style="width:{100.0 * tally[s] / total:.1f}%"></span>
          </span>
          <span class="tally-count">{tally[s]}<span class="of">/{total}</span></span>
        </div>'''
        for s in sorted(SOLVERS, key=lambda k: -tally[k])
    )

    body = []
    for name, result in rows:
        cells = []
        for s in SOLVERS:
            entry = result[s]
            solved = entry["status"] == "optimal"
            seconds = entry["seconds"]
            klass = "cell" + (" solved" if solved else " missed")
            klass += " subject" if s == "ripsolve" else ""
            figure = f"{seconds:.1f}s" if solved else "—"
            cells.append(
                f'''          <td class="{klass}">
            <span class="bar" style="width:{width(seconds, solved):.1f}%"></span>
            <span class="figure">{figure}</span>
          </td>'''
            )
        ours = result["ripsolve"]["status"] == "optimal"
        body.append(
            f'''        <tr class="{'we-solved' if ours else 'we-missed'}">
          <th scope="row">{html.escape(name)}</th>
{chr(10).join(cells)}
        </tr>'''
        )

    fastest = [n for n, r in rows
               if r["ripsolve"]["status"] == "optimal"
               and r["ripsolve"]["seconds"] < min(
                   r[s]["seconds"] for s in ("HiGHS", "SCIP", "CBC")
                   if r[s]["status"] == "optimal")]

    page = PAGE.format(
        total=total,
        ours=tally["ripsolve"],
        limit=int(LIMIT),
        threads=DATA["threads"],
        summary=summary,
        rows="\n".join(body),
        headers="\n".join(
            f'          <th scope="col"{" class=\"subject\"" if s == "ripsolve" else ""}>'
            f'{html.escape(LABEL[s])}</th>' for s in SOLVERS),
        fastest=", ".join(html.escape(n) for n in fastest) or "none",
        fastest_count=len(fastest),
        missed=total - tally["ripsolve"],
    )
    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text(page)
    print(f"wrote {DEST}")


PAGE = r"""<title>The Tractable Twenty</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;700&family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500&display=swap">
<style>
  :root {{
    --paper: #F4F6F7;
    --card: #FFFFFF;
    --ink: #131C24;
    --muted: #5B6A75;
    --rule: #D7DEE2;
    --rip: #B4531A;
    --rip-soft: #E8D2C2;
    --field: #2F6F8F;
    --field-soft: #CBDBE4;
    --ref: #6B7B85;
    --miss: #E2E7EA;
    --shadow: 0 1px 2px rgba(19, 28, 36, .06), 0 8px 24px rgba(19, 28, 36, .05);
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --paper: #0E151B;
      --card: #16202A;
      --ink: #E5EBEF;
      --muted: #93A3AE;
      --rule: #26333E;
      --rip: #E08A4C;
      --rip-soft: #4A2E1B;
      --field: #6FB0CE;
      --field-soft: #1E3947;
      --ref: #8496A2;
      --miss: #1E2A34;
      --shadow: 0 1px 2px rgba(0, 0, 0, .4), 0 8px 24px rgba(0, 0, 0, .3);
    }}
  }}
  :root[data-theme="dark"] {{
    --paper: #0E151B;
    --card: #16202A;
    --ink: #E5EBEF;
    --muted: #93A3AE;
    --rule: #26333E;
    --rip: #E08A4C;
    --rip-soft: #4A2E1B;
    --field: #6FB0CE;
    --field-soft: #1E3947;
    --ref: #8496A2;
    --miss: #1E2A34;
    --shadow: 0 1px 2px rgba(0, 0, 0, .4), 0 8px 24px rgba(0, 0, 0, .3);
  }}

  body {{
    margin: 0;
    background: var(--paper);
    color: var(--ink);
    font: 400 16px/1.6 "IBM Plex Sans", system-ui, sans-serif;
    -webkit-font-smoothing: antialiased;
  }}
  .wrap {{
    max-width: 1080px;
    margin: 0 auto;
    padding: clamp(2rem, 5vw, 4.5rem) clamp(1rem, 4vw, 2.5rem) 5rem;
    display: flex;
    flex-direction: column;
    gap: clamp(2rem, 4vw, 3.25rem);
  }}

  .eyebrow {{
    font: 500 .74rem/1 "IBM Plex Mono", ui-monospace, monospace;
    letter-spacing: .14em;
    text-transform: uppercase;
    color: var(--muted);
    margin: 0 0 1rem;
  }}
  h1 {{
    font: 700 clamp(2.1rem, 5.5vw, 3.4rem)/1.05 Archivo, system-ui, sans-serif;
    letter-spacing: -.022em;
    text-wrap: balance;
    margin: 0 0 1rem;
  }}
  .standfirst {{
    max-width: 62ch;
    color: var(--muted);
    margin: 0;
    font-size: 1.05rem;
  }}
  .standfirst strong {{ color: var(--ink); font-weight: 500; }}

  .panel {{
    background: var(--card);
    border: 1px solid var(--rule);
    border-radius: 10px;
    box-shadow: var(--shadow);
    padding: clamp(1.25rem, 3vw, 2rem);
  }}
  h2 {{
    font: 700 1.15rem/1.2 Archivo, system-ui, sans-serif;
    letter-spacing: -.01em;
    margin: 0 0 .4rem;
  }}
  .note {{ color: var(--muted); font-size: .92rem; margin: 0 0 1.5rem; max-width: 64ch; }}

  .tally {{ display: flex; flex-direction: column; gap: .7rem; }}
  .tally-row {{
    display: grid;
    grid-template-columns: 6.5rem 1fr 4rem;
    align-items: center;
    gap: .9rem;
  }}
  .tally-name {{
    font: 500 .92rem/1 "IBM Plex Mono", ui-monospace, monospace;
    color: var(--muted);
  }}
  .is-subject .tally-name {{ color: var(--rip); }}
  .tally-track {{
    background: var(--miss);
    border-radius: 3px;
    height: 1.4rem;
    overflow: hidden;
  }}
  .tally-fill {{
    display: block;
    height: 100%;
    background: var(--field);
    border-radius: 3px;
  }}
  .is-subject .tally-fill {{ background: var(--rip); }}
  .tally-count {{
    font: 500 1rem/1 "IBM Plex Mono", ui-monospace, monospace;
    font-variant-numeric: tabular-nums;
    text-align: right;
  }}
  .of {{ color: var(--muted); font-weight: 400; }}

  .scroller {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; width: 100%; min-width: 720px; }}
  thead th {{
    font: 500 .74rem/1 "IBM Plex Mono", ui-monospace, monospace;
    letter-spacing: .1em;
    text-transform: uppercase;
    color: var(--muted);
    text-align: left;
    padding: 0 .5rem .7rem;
    border-bottom: 1px solid var(--rule);
  }}
  thead th.subject {{ color: var(--rip); }}
  tbody th {{
    font: 400 .88rem/1.3 "IBM Plex Mono", ui-monospace, monospace;
    text-align: left;
    font-weight: 400;
    padding: .5rem .9rem .5rem 0;
    white-space: nowrap;
  }}
  tbody tr {{ border-bottom: 1px solid var(--rule); }}
  tbody tr:last-child {{ border-bottom: 0; }}
  .we-solved th {{ color: var(--rip); }}
  .cell {{
    position: relative;
    padding: .5rem .5rem;
    width: 17%;
    vertical-align: middle;
  }}
  .bar {{
    display: block;
    height: .5rem;
    border-radius: 2px;
    background: var(--field-soft);
    margin-bottom: .25rem;
  }}
  .cell.solved .bar {{ background: var(--field); }}
  .cell.subject.solved .bar {{ background: var(--rip); }}
  .cell.missed .bar {{
    background: repeating-linear-gradient(
      -45deg, var(--miss), var(--miss) 3px, transparent 3px, transparent 6px);
    border: 1px solid var(--miss);
    box-sizing: border-box;
  }}
  .figure {{
    font: 400 .78rem/1 "IBM Plex Mono", ui-monospace, monospace;
    font-variant-numeric: tabular-nums;
    color: var(--muted);
  }}
  .cell.missed .figure {{ opacity: .55; }}

  .legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 1.25rem;
    margin-top: 1.25rem;
    font: 400 .82rem/1 "IBM Plex Sans", sans-serif;
    color: var(--muted);
  }}
  .key {{ display: inline-flex; align-items: center; gap: .45rem; }}
  .swatch {{ width: 1.6rem; height: .5rem; border-radius: 2px; }}
  .swatch.solved {{ background: var(--field); }}
  .swatch.subject {{ background: var(--rip); }}
  .swatch.missed {{
    background: repeating-linear-gradient(
      -45deg, var(--miss), var(--miss) 3px, transparent 3px, transparent 6px);
    border: 1px solid var(--rule);
    box-sizing: border-box;
  }}

  footer {{
    color: var(--muted);
    font-size: .88rem;
    border-top: 1px solid var(--rule);
    padding-top: 1.5rem;
    max-width: 68ch;
  }}
  code {{
    font: 400 .88em/1 "IBM Plex Mono", ui-monospace, monospace;
    background: var(--miss);
    padding: .12em .38em;
    border-radius: 3px;
  }}
</style>

<div class="wrap">
  <header>
    <p class="eyebrow">MIPLIB 2017 · {limit}s limit · {threads} threads</p>
    <h1>Twenty instances the field can close</h1>
    <p class="standfirst">Screened from the MIPLIB <em>easy</em> list by one rule: at least two of
    HiGHS, SCIP and CBC prove optimality within the budget. On a set chosen that way, every
    failure belongs to the solver rather than the instance, and every one is known to be
    reachable. <strong>ripsolve closes {ours} of {total}.</strong></p>
  </header>

  <section class="panel">
    <h2>Instances proved optimal</h2>
    <p class="note">Out of {total}. CBC sets the low bar among the three that defined the
    set, since qualifying needed only two of them to agree.</p>
    <div class="tally">
{summary}
    </div>
  </section>

  <section class="panel">
    <h2>Time to prove optimality, by instance</h2>
    <p class="note">Bars are log-scaled from 0.05s to the {limit}s limit, so a bar twice as
    long is far more than twice the time. A hatched bar means the solver ran out the clock.</p>
    <div class="scroller">
      <table>
        <thead>
          <tr>
            <th scope="col">instance</th>
{headers}
          </tr>
        </thead>
        <tbody>
{rows}
        </tbody>
      </table>
    </div>
    <div class="legend">
      <span class="key"><span class="swatch subject"></span> ripsolve, proved</span>
      <span class="key"><span class="swatch solved"></span> other solver, proved</span>
      <span class="key"><span class="swatch missed"></span> hit the time limit</span>
    </div>
  </section>

  <footer>
    <p>Where ripsolve does finish it is competitive: it is the fastest of the four
    open-source solvers on {fastest_count} of the {ours} it closes ({fastest}). The gap is not
    that it is uniformly slow, but that {missed} of {total} instances hit something it handles
    badly. Reproduce with <code>bench/tractable.py</code> to select the set and
    <code>bench/tractable_chart.py</code> to measure it; every solver gets the same limit and
    thread count, and only the solve is timed.</p>
  </footer>
</div>
"""

main()
