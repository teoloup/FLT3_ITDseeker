#!/usr/bin/env python3
"""
Render GUIDE.md as a standalone styled GUIDE.html.

The HTML is generated rather than maintained by hand so it cannot drift from the
markdown. No external dependency: this handles the subset of markdown GUIDE.md
actually uses -- headings, tables, fenced code, lists, inline code, bold, links,
horizontal rules -- rather than pulling in a full parser.

Usage:  python build_guide_html.py [in.md] [out.html]
"""

import html
import re
import sys


def inline(text):
    """Inline spans.

    Code spans are lifted out to placeholders rather than handled by splitting,
    so that bold and links can still span them -- `**A caveat on `POS`.**` and
    ``[`benchmarks/`](benchmarks/)`` both occur in the guide and both break if
    the string is split on backticks first.
    """
    codes = []

    def stash(m):
        codes.append(f"<code>{html.escape(m.group(1))}</code>")
        return f"\x00{len(codes) - 1}\x00"

    p = re.sub(r"`([^`]+)`", stash, text)
    p = html.escape(p)
    p = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', p)
    p = re.sub(r"\*\*([^*]+?)\*\*", r"<strong>\1</strong>", p)
    p = re.sub(r"(?<![\w*])\*([^*]+?)\*(?![\w*])", r"<em>\1</em>", p)
    return re.sub(r"\x00(\d+)\x00", lambda m: codes[int(m.group(1))], p)


def slug(text):
    s = re.sub(r"[^\w\s-]", "", text.lower()).strip()
    return re.sub(r"[\s_]+", "-", s)


def convert(md):
    lines = md.split("\n")
    out, toc = [], []
    i = 0
    in_code = False
    list_stack = []

    def close_lists():
        while list_stack:
            out.append(f"</{list_stack.pop()}>")

    while i < len(lines):
        line = lines[i]

        # fenced code
        m = re.match(r"^```(\w*)\s*$", line)
        if m:
            if in_code:
                out.append("</code></pre>")
                in_code = False
            else:
                close_lists()
                out.append(f'<pre class="lang-{m.group(1) or "text"}"><code>')
                in_code = True
            i += 1
            continue
        if in_code:
            out.append(html.escape(line))
            i += 1
            continue

        # table: a header row followed by a separator row
        if (line.strip().startswith("|") and i + 1 < len(lines)
                and re.match(r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1])):
            close_lists()
            def cells(row):
                return [c.strip() for c in row.strip().strip("|").split("|")]
            head = cells(line)
            aligns = ["right" if c.strip().endswith(":") and not c.strip().startswith(":")
                      else "left" for c in cells(lines[i + 1])]
            out.append("<table><thead><tr>")
            for j, h in enumerate(head):
                out.append(f'<th style="text-align:{aligns[j]}">{inline(h)}</th>')
            out.append("</tr></thead><tbody>")
            i += 2
            while i < len(lines) and lines[i].strip().startswith("|"):
                out.append("<tr>")
                for j, c in enumerate(cells(lines[i])):
                    a = aligns[j] if j < len(aligns) else "left"
                    out.append(f'<td style="text-align:{a}">{inline(c)}</td>')
                out.append("</tr>")
                i += 1
            out.append("</tbody></table>")
            continue

        # headings
        m = re.match(r"^(#{1,4})\s+(.*)$", line)
        if m:
            close_lists()
            level, text = len(m.group(1)), m.group(2)
            anchor = slug(text)
            if level == 2:
                toc.append((text, anchor))
            out.append(f'<h{level} id="{anchor}">{inline(text)}</h{level}>')
            i += 1
            continue

        if re.match(r"^---+\s*$", line):
            close_lists()
            out.append("<hr>")
            i += 1
            continue

        # lists
        m = re.match(r"^(\s*)([-*]|\d+\.)\s+(.*)$", line)
        if m:
            indent = len(m.group(1))
            tag = "ul" if m.group(2) in "-*" else "ol"
            depth = indent // 2
            while len(list_stack) > depth + 1:
                out.append(f"</{list_stack.pop()}>")
            if len(list_stack) == depth:
                out.append(f"<{tag}>")
                list_stack.append(tag)
            out.append(f"<li>{inline(m.group(3))}</li>")
            i += 1
            continue

        if not line.strip():
            close_lists()
            i += 1
            continue

        close_lists()
        para = [line]
        i += 1
        while (i < len(lines) and lines[i].strip()
               and not re.match(r"^(#{1,4}\s|```|\s*([-*]|\d+\.)\s|\||---+\s*$)", lines[i])):
            para.append(lines[i])
            i += 1
        out.append(f"<p>{inline(' '.join(para))}</p>")

    if in_code:
        out.append("</code></pre>")
    close_lists()
    return "\n".join(out), toc


CSS = """
:root { --bg:#f6f7f9; --panel:#fff; --ink:#1c2430; --muted:#6b7684;
  --line:#e3e7ec; --accent:#2f6f4f; --accent-soft:#e7f2ec; --code:#f4f6f8; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
  font:16px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }
.layout { display:grid; grid-template-columns:250px minmax(0,1fr); gap:40px;
  max-width:1160px; margin:0 auto; padding:32px 20px 80px; }
nav { position:sticky; top:32px; align-self:start; font-size:14px; }
nav h2 { font-size:11px; text-transform:uppercase; letter-spacing:.09em;
  color:var(--muted); margin:0 0 10px; }
nav a { display:block; padding:5px 0; color:var(--ink); text-decoration:none;
  border-left:2px solid var(--line); padding-left:12px; }
nav a:hover { border-left-color:var(--accent); color:var(--accent); }
main { background:var(--panel); border:1px solid var(--line); border-radius:12px;
  padding:34px 42px 48px; min-width:0; }
h1 { font-size:30px; margin:0 0 6px; letter-spacing:-.01em; }
h2 { font-size:22px; margin:44px 0 14px; padding-top:18px;
  border-top:1px solid var(--line); }
h3 { font-size:17px; margin:28px 0 10px; }
h4 { font-size:15px; margin:22px 0 8px; color:var(--muted); }
p { margin:12px 0; }
a { color:var(--accent); }
code { background:var(--code); border:1px solid var(--line); border-radius:4px;
  padding:1px 5px; font-size:.88em;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
pre { background:var(--code); border:1px solid var(--line); border-radius:8px;
  padding:14px 16px; overflow-x:auto; line-height:1.5; }
pre code { background:none; border:none; padding:0; font-size:13px; }
table { border-collapse:collapse; width:100%; margin:18px 0; font-size:14.5px;
  display:block; overflow-x:auto; }
th { text-align:left; font-size:11px; text-transform:uppercase; letter-spacing:.05em;
  color:var(--muted); border-bottom:2px solid var(--line); padding:0 12px 8px 0; }
td { padding:9px 12px 9px 0; border-bottom:1px solid var(--line);
  vertical-align:top; }
tr:last-child td { border-bottom:none; }
ul,ol { padding-left:22px; margin:12px 0; }
li { margin:5px 0; }
hr { border:none; border-top:1px solid var(--line); margin:34px 0; }
strong { font-weight:650; }
@media (max-width:900px) {
  .layout { grid-template-columns:1fr; gap:20px; }
  nav { position:static; }
  main { padding:24px 20px 32px; }
}
@media print {
  body { background:#fff; } nav { display:none; }
  .layout { display:block; max-width:none; padding:0; }
  main { border:none; padding:0; } h2 { break-after:avoid; }
  pre,table { break-inside:avoid; }
}
"""


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "GUIDE.md"
    dst = sys.argv[2] if len(sys.argv) > 2 else "GUIDE.html"
    md = open(src, encoding="utf-8").read()

    # the H1 becomes the page title; the markdown's own Contents list is dropped
    # in favour of the sticky sidebar
    md = re.sub(r"^## Contents\n(?:.*\n)*?(?=^---$)", "", md, flags=re.M)
    body, toc = convert(md)

    m = re.search(r"<h1[^>]*>(.*?)</h1>", body)
    title = re.sub(r"<[^>]+>", "", m.group(1)) if m else "Guide"
    nav = "".join(f'<a href="#{a}">{html.escape(t)}</a>' for t, a in toc)

    out = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="layout">
  <nav><h2>Contents</h2>{nav}</nav>
  <main>{body}</main>
</div>
</body>
</html>
"""
    with open(dst, "w", encoding="utf-8", newline="\n") as f:
        f.write(out)
    print(f"wrote {dst} ({len(out):,} bytes, {len(toc)} sections)")


if __name__ == "__main__":
    main()
