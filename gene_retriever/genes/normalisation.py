import re
# --------------------------------------------------------------------------
# Normalization / filename helpers
# --------------------------------------------------------------------------

def normalise_hyphens(s: str) -> str:
    """Replace various unicode hyphen-like characters with ASCII '-'."""
    if not isinstance(s, str):
        s = str(s)
    for ch in ["\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2212"]:
        s = s.replace(ch, "-")
    return s

def normalise_apostrophes(s: str) -> str:
    """Replace various apostrophe-like characters with ASCII apostrophe '."""
    if not isinstance(s, str):
        s = str(s)
    for ch in ["\u2018", "\u2019", "\u02BC", "\u2032", "’", "‘", "ʹ", "ʼ", "ʾ"]:
        s = s.replace(ch, "'")
    return s

def strip_wrapping_quotes(s: str) -> str:
    """Remove a single pair of wrapping quotes if present (leading+trailing)."""
    if not isinstance(s, str):
        s = str(s)
    s = s.strip()
    if len(s) >= 2 and ((s[0] == '"' and s[-1] == '"') or (s[0] == "'" and s[-1] == "'")):
        s = s[1:-1].strip()
    if s[-1] == ".":
        s = s[:-1].strip()
    return s
    
def canonical_gene_name(g: str) -> str:
    """Convert a raw gene string into its canonical biological form."""
    g = strip_wrapping_quotes(g)
    g = normalise_hyphens(g)
    g = normalise_apostrophes(g)
    g = g.strip()
    return g

def clean_gene_term(gene):
    """Light canonicalisation used before generating patterns."""
    if gene is None:
        return ""
    s = canonical_gene_name(str(gene).strip())
    # remove trailing ' like' or ' type'
    for pattern in [r'\s*[-\s]*like\s*$', r'\s*[-\s]*type\s*$']:
        s = re.sub(pattern, '', s, flags=re.IGNORECASE)
    s = s.strip()
    # keep only first token if extra words are present: 23S rRNA --> 23S
    parts = s.split()
    if len(parts) > 1:
        s = parts[0]
    return s

# --------------------------------------------------------------------------
# Pattern generation for gene names
# --------------------------------------------------------------------------

def generate_gene_patterns(query_gene):
    """
    Generate a list of normalized patterns (lowercase) to try for a gene name.
    Handles generic separators and special 'bla' forms.
    """
    q = clean_gene_term(query_gene)
    q_lower = q.lower()
    patterns = []
    # special handling for bla prefix
    bla_match = re.match(r'^bla[_-]?(.*)$', q_lower, re.IGNORECASE)
    if bla_match:
        core = bla_match.group(1) or ""
        core = normalise_hyphens(core)
        # core variations
        patterns.extend([
            core,
            core.replace('-', '_'),
            core.replace('_', '-'),
            core.replace('-', '').replace('_', '')
        ])
        core_us = core.replace('-', '_')
        core_hy = core.replace('_', '-')
        core_ns = core.replace('-', '').replace('_', '')
        patterns.extend([
            f"bla{core}", f"bla-{core}", f"bla_{core}",
            f"bla{core_us}", f"bla{core_hy}", f"bla{core_ns}",
            f"bla-{core_us}", f"bla-{core_hy}",
            f"bla_{core_us}", f"bla_{core_hy}",
        ])
    else:
        core = normalise_hyphens(q_lower)
        patterns.extend([
            core,
            core.replace('-', '_'),
            core.replace('_', '-'),
            core.replace('-', '').replace('_', ''),
        ])
    # deduplicate preserving order
    seen = set()
    out = []
    for p in patterns:
        if not p:
            continue
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out
