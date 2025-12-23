#!/usr/bin/env python3
"""
normalize_genes.py

Reads lines (genes) from stdin or a file, normalises unicode variants and
quotes/primes/hyphens, removes surrounding double-quotes, and writes a
sorted unique list to stdout (one gene per line).

Usage:
  # from file:
  python3 normalize_genes.py GENES_NOT_FOUND.txt > GENES_NOT_FOUND.normalized.txt

  # or from a pipe:
  cat GENES_NOT_FOUND.txt | python3 normalize_genes.py > GENES_NOT_FOUND.normalized.txt

  # Do:
  cut -f6 3_run_stats.tsv | tr "," "\n"   | grep -v '^GENES_NOT_FOUND$' | grep -v '^$'   | python3 normalize_genes.py > GENES_NOT_FOUND.normalized.txt
  # Instead of:
  cut -f6 3_run_stats.tsv | tr "," "\n"   | grep -v '^GENES_NOT_FOUND$' | grep -v '^$' > GENES_NOT_FOUND.txt
"""

import sys
import re

# characters to replace with ASCII hyphen
HYPHENS = [
    "\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2212",  # various dashes
    "\u00AD"  # soft hyphen (rare)
]

# curly quotes and other quote-like characters
DOUBLE_QUOTES = [
    '\u0022',  # "
    '\u201C',  # “
    '\u201D',  # ”
    '\u201E',  # „
    '\u201F',  # ‟
]
SINGLE_QUOTES = [
    "\u2018", "\u2019", "\u201B",  # ‘ ’ ‛
    "\u02BC", "\u2032", "\u2035"   # modifier apostrophes/primes (some also map to ')
]

DOUBLE_PRIME_CHARS = [
    "\u2033", "\u2036",  # ″ ″ variants
    "\u301E", "\u301D",  # other quote-like
]

def normalise_hyphens(s: str) -> str:
    for ch in HYPHENS:
        s = s.replace(ch, "-")
    return s

def normalise_primes_and_quotes(s: str) -> str:
    # First: replace double-prime unicode with ASCII double single-quotes -> ''
    for ch in DOUBLE_PRIME_CHARS:
        s = s.replace(ch, "''")

    # Replace two ASCII double-quotes "": -> ''
    s = s.replace('""', "''")

    # Replace remaining double-prime-like sequences: a pair of double-quotes or two double-prime chars already handled.
    # Now convert any remaining unicode single-primes/curly single quotes -> ASCII apostrophe '
    for ch in SINGLE_QUOTES:
        s = s.replace(ch, "'")

    # Convert any remaining curly double quotes to plain double-quote (we'll trim leading/trailing later)
    for ch in DOUBLE_QUOTES:
        s = s.replace(ch, '"')

    # Convert `′` (u2032) if present (it is in SINGLE_QUOTES above) -> '
    # Already handled.

    return s

def strip_edge_double_quotes(s: str) -> str:
    # Remove leading/trailing double quote characters (ASCII " ) only,
    # keep internal quotes which may be meaningful (but we've converted "" -> '')
    s = s.strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1].strip()
    # also handle leading or trailing stray quotes
    s = s.lstrip('"').rstrip('"').strip()
    return s

def normalize_one(line: str) -> str:
    s = line.rstrip("\n\r")
    s = s.strip()
    if not s:
        return ""

    # Normalise hyphens first
    s = normalise_hyphens(s)

    # Normalise primes/apostrophes/double-primes, curly quotes etc
    s = normalise_primes_and_quotes(s)

    # Remove leading/trailing ASCII double quote characters (preserve interior quotes already converted)
    s = strip_edge_double_quotes(s)

    # Remove non-printable control characters (keep common punctuation)
    s = re.sub(r'[\x00-\x1f\x7f]', '', s)

    # Normalise spacing
    s = re.sub(r'\s+', ' ', s).strip()

    # Some tidy-ups: replace any sequence of hyphens/spaces with single hyphen if it looks like a gene token part
    # (don't over-normalise; keep original tokens mostly intact)
    # But ensure we don't accidentally join tokens like "OXA - 114" -> "OXA-114"
    s = re.sub(r'\s*-\s*', '-', s)

    # Replace any remaining weird double apostrophes like: (3'') but sometimes appears as (3’') or with spaces
    s = s.replace("’''", "''").replace("''", "''")

    return s

def main(argv):
    infile = None
    if len(argv) > 1:
        infile = argv[1]

    lines = []
    if infile:
        with open(infile, 'r', encoding='utf-8', errors='replace') as fh:
            lines = fh.readlines()
    else:
        lines = sys.stdin.readlines()

    out_set = set()
    for ln in lines:
        for g in ln.strip().split(","):
            norm = normalize_one(ln)
            if not norm:
                continue
            out_set.add(norm)

    # Sort nicely: put ASCII-order but keep case-sensitivity
    out_list = sorted(out_set, key=lambda x: x.lower())

    # Write to stdout
    for v in out_list:
        print(v)

if __name__ == "__main__":
    main(sys.argv)
