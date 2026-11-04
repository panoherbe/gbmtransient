"""Bericht als PDF: nachweis/<Datum>/BERICHT.md (und PROTOKOLL.md als letzter Abschnitt) -> BERICHT.pdf.

Markdown -> HTML mit stil.css -> PDF über Google Chrome ohne Oberfläche; Abbildungen werden eingebettet.
Braucht die Pakete markdown und pypdf und Google Chrome (Pfad über --chrome).

    python werkzeuge/pdf.py nachweis/2026-11-04 --chrome "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
"""

import argparse
import base64
import html
import re
import subprocess
import tempfile
from pathlib import Path

HIER = Path(__file__).resolve().parent
AUTOR = "Panagiotis Herbe"


def einbetten(md: str, basis: Path) -> str:
    def ersetzen(m):
        daten = base64.b64encode((basis / m.group(2)).read_bytes()).decode()
        return f"![{m.group(1)}](data:image/png;base64,{daten})"
    return re.sub(r"!\[([^\]]*)\]\(([^)]+\.png)\)", ersetzen, md)


def nach_html(md: str) -> str:
    import markdown
    body = markdown.markdown(md, extensions=["tables", "attr_list", "sane_lists", "fenced_code"])
    body = re.sub(r"<table>\s*<thead>\s*<tr>(\s*<th[^>]*>\s*</th>)+\s*</tr>\s*</thead>", '<table class="kv">', body)
    body = re.sub(r"(\d) (%|σ|s\b|h\b|GB|keV)", "\\1&nbsp;\\2", body)
    return body


def seite(body: str, titel: str, kopf: str, fuss: str) -> str:
    css = (HIER / "stil.css").read_text(encoding="utf-8")
    rand = f"""
@page {{
  @top-left {{ content: "{kopf}"; }}
  @top-right {{ content: "{AUTOR}"; }}
  @bottom-left {{ content: "{fuss}"; }}
  @bottom-right {{ content: "Seite " counter(page) " von " counter(pages); }}
}}
p > img {{ max-width: 100%; }}
"""
    return (f'<!doctype html><html lang="de"><head><meta charset="utf-8"><title>{html.escape(titel)}</title>'
            f"<style>{css}{rand}</style></head><body>{body}</body></html>")


def drucken(html_text: str, pdf: Path, titel: str, thema: str, chrome: str) -> int:
    from pypdf import PdfReader, PdfWriter
    with tempfile.TemporaryDirectory() as tmp:
        quelle = Path(tmp) / "bericht.html"
        quelle.write_text(html_text, encoding="utf-8")
        roh = Path(tmp) / "roh.pdf"
        subprocess.run([chrome, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                        "--generate-pdf-document-outline", f"--print-to-pdf={roh}", quelle.as_uri()],
                       check=True, capture_output=True, timeout=180)
        w = PdfWriter(clone_from=PdfReader(str(roh)))
        w.add_metadata({"/Title": titel, "/Author": AUTOR, "/Subject": thema})
        with open(pdf, "wb") as f:
            w.write(f)
    return len(PdfReader(str(pdf)).pages)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("nachweis", type=Path)
    ap.add_argument("--chrome", default="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
    ap.add_argument("--fuss", default="")
    a = ap.parse_args()
    md = (a.nachweis / "BERICHT.md").read_text(encoding="utf-8")
    prot = a.nachweis / "PROTOKOLL.md"
    if prot.exists():
        p = prot.read_text(encoding="utf-8")
        p = re.sub(r"^# (.*)$", r"## 8 \1", p, count=1, flags=re.M)
        p = re.sub(r"^## (?!8 )(.*)$", r"### \1", p, flags=re.M)
        md = md.rstrip() + "\n\n" + p
    titel = re.search(r"^# (.*)$", md, re.M).group(1)
    body = nach_html(einbetten(md, a.nachweis))
    n = drucken(seite(body, titel, "gbmtransient · Technischer Bericht", a.fuss), a.nachweis / "BERICHT.pdf", titel,
                "Vorhaben von Panagiotis Herbe", a.chrome)
    print(f"{n} Seiten: {a.nachweis / 'BERICHT.pdf'}")


if __name__ == "__main__":
    main()
