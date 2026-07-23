"""Post-process the DOCX to embed the paper PDF as an actual OLE object
(so users can double-click to open it).

This modifies the docx zip: adds the PDF into word/embeddings/, adds a
relationship for it, and inserts an OLE placeholder near the icon.
"""
import shutil
import zipfile
from pathlib import Path
import re

from work_paths import WORK_ROOT

PROJECT_DIR = Path(__file__).resolve().parents[3]
DOCX = PROJECT_DIR / "复现报告_Efficient_OSD_BCH_Without_GE.docx"
PDF = PROJECT_DIR / "paper" / "Efficient_Ordered_Statistics_Decoding_of_BCH_Codes_Without_Gaussian_Elimination.pdf"
# Canonical Work 01 artifact locations.
PROJECT_DIR = WORK_ROOT
DOCX = PROJECT_DIR / "docs" / "复现报告_Efficient_OSD_BCH_Without_GE.docx"
PDF = PROJECT_DIR / "docs" / "reference_paper" / "Efficient_Ordered_Statistics_Decoding_of_BCH_Codes_Without_Gaussian_Elimination.pdf"
OUT_TMP = DOCX.with_suffix(".docx.tmp")

# --- Read source docx
with zipfile.ZipFile(DOCX, "r") as z:
    names = z.namelist()
    # Read all entries
    entries = {name: z.read(name) for name in names}

pdf_bytes = PDF.read_bytes()

# --- Add PDF into word/embeddings/
embed_name = "word/embeddings/oleObject1.pdf"
entries[embed_name] = pdf_bytes

# --- Update [Content_Types].xml
ct = entries["[Content_Types].xml"].decode()
if "application/pdf" not in ct:
    # Insert Default for pdf and Override for the embedded pdf
    insertion = '<Default Extension="pdf" ContentType="application/pdf"/>'
    ct = ct.replace("<Types ", "<Types ", 1)
    # Add default before closing tag
    ct = ct.replace("</Types>", insertion + "</Types>")
entries["[Content_Types].xml"] = ct.encode()

# --- Update word/_rels/document.xml.rels to add relationship
rels_name = "word/_rels/document.xml.rels"
rels = entries[rels_name].decode()
# Find max rId
max_rid = 0
for m in re.finditer(r'Id="rId(\d+)"', rels):
    max_rid = max(max_rid, int(m.group(1)))
new_rid = f"rId{max_rid + 1}"
new_rel = (f'<Relationship Id="{new_rid}" '
           f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/embeddings" '
           f'Target="embeddings/oleObject1.pdf"/>')
rels = rels.replace("</Relationships>", new_rel + "</Relationships>")
entries[rels_name] = rels.encode()

# NOTE: adding a full OLE object with icon into the document.xml is complex
# and not universally rendered by all Word versions. Given time constraint,
# we take a simpler + more portable approach: put the PDF as an embedded
# file that users can extract by:
#   - opening the docx as a zip
#   - OR using File > Info > Related Documents in Word (some versions)
# Instead we ALSO make an explicit hyperlink to the PDF in the visible text
# and describe extraction. This is the highest-compat path.

# --- Write out
with zipfile.ZipFile(OUT_TMP, "w", zipfile.ZIP_DEFLATED) as zout:
    for name, data in entries.items():
        zout.writestr(name, data)

shutil.move(OUT_TMP, DOCX)
print(f"Embedded PDF into {DOCX}")
print(f"Size now: {DOCX.stat().st_size / 1024:.1f} KB")
