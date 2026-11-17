# mergePdfs.py
#
#   Merges all pdfs in cwd and puts it into outputPdf filename.
#
import os
import re
from PyPDF2 import PdfMerger
import time                     # time tracker

outputPdf = "combined.pdf"

timeStartScript = time.time()

def sort_key(filename):
    m = re.match(r"plt(\d+)_", filename, re.IGNORECASE)
    return int(m.group(1)) if m else 999999

pdf_files = [
    f for f in os.listdir(".")
    if f.lower().endswith(".pdf")
    and f.lower() != outputPdf.lower()
]

pdf_files.sort(key=sort_key)

merger = PdfMerger()

for pdf in pdf_files:
    print(f"Adding {pdf}")
    merger.append(pdf)

merger.write(outputPdf)
merger.close()

print(f"\nCreated {outputPdf}")

elapsed = time.time() - timeStartScript

appLog.info(f"Total runtime: {elapsed:.1f} sec ({elapsed/60:.1f} min)")
