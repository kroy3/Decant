#!/usr/bin/env python3
"""Copy the user guide from wiki/ into rtd/guide/ for the Read the Docs build.

GitHub wiki links name a page without an extension ("](03-Quick-Start)");
MkDocs needs the file ("](03-Quick-Start.md)"). Only those links are
rewritten; everything else is copied as is. Run from the repository root.
"""
import pathlib
import re
import shutil

root = pathlib.Path(__file__).resolve().parent.parent
src, dst = root / "wiki", root / "rtd" / "guide"
if dst.exists():
    shutil.rmtree(dst)
dst.mkdir(parents=True)
pages = sorted(p for p in src.glob("*.md") if not p.name.startswith("_"))
names = {p.stem for p in pages}
link = re.compile(r"\]\(([A-Za-z0-9-]+)(#[^)]*)?\)")

def fix(m):
    name, anchor = m.group(1), m.group(2) or ""
    return f"]({name}.md{anchor})" if name in names else m.group(0)

for p in pages:
    (dst / p.name).write_text(link.sub(fix, p.read_text(encoding="utf-8")), encoding="utf-8")
print(f"copied {len(pages)} pages to {dst.relative_to(root)}")
