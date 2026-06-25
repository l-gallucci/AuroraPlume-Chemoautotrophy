# Sulfurimonas 16S Phylogenetic Placement — PS137 Aurora Vent Field

**Version:** 5-final | **Author:** lgallucc

---

## Overview

Builds a full-length 16S rRNA reference phylogeny for *Sulfurimonas* from hydrothermal vent and type strain sequences, then places short V3–V4 ASVs from the PS137 expedition onto the reference tree using SEPP.

---

## Software Requirements

| Tool | Version | Role |
|------|---------|------|
| Python | 3.13.x | scripting |
| BioPython | 1.84 | sequence I/O |
| NumPy | 2.0.x | alignment array ops |
| SciPy | 1.14.x | hierarchical clustering |
| barrnap | ≥ 1.0 (HRGV fork) | 16S rRNA prediction (Infernal CMs) |
| SINA | 1.7.2 (stable) | structure-guided alignment QC |
| CD-HIT | 4.8.1 | dereplication at 98.6% |
| seqkit | 2.13.0 | FASTA filtering / extraction |
| MAFFT | 7.526 | `--add` profile alignment |
| IQ-TREE3 | 3.1.2 | ML tree + ModelFinder + UFBoot |
| bedtools | 2.31.1 | genome 16S extraction |
| RAxML | 8.2.12 (`raxmlHPC`) | GTRGAMMA optimisation for SEPP |
| SEPP | 4.5.1 | phylogenetic placement (ensemble HMM) |
| pplacer | 1.1.alpha19+ | backend for SEPP |
| DendroPy | 4.6.1 | RF distance + branch-length analysis |

> **SEPP environment:** `conda activate SEPP`. Replace bundled pplacer if needed:
> `cp $(which pplacer) ~/.sepp/bundled-v*/pplacer`

---

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| Reference set: hydrothermal vent + type strains only | Focused on ecologically relevant diversity |
| Outgroup: *Sulfuricurvum* only (~10 seqs) | Caminibacter caused LBA (confirmed by RF analysis; B vs C RF = 0% with *Sulfuricurvum*) |
| SILVA `full_align_trunc` pre-alignment | Structure-guided; genomes/outgroup added via `mafft --add` (preserves profile) |
| No trimAl | SILVA alignment already curated; trimming reduces accuracy for single-marker genes (Tan et al. 2015) |
| Dereplication at 98.6% | Species-level threshold for 16S rRNA (Kim et al. 2014) |
| Duplicate accessions deduplicated before IQ-TREE3 | Keep longest copy per genome |
| IQ-TREE3: MFP + UFBoot 1000 + `--bnni` + SH-aLRT 1000 | Best-fit model: TVM+R4 (ModelFinder BIC) |
| Long branch flagging at mean + 3 SD | Conservative threshold |
| SEPP for placement | Ensemble HMM; better than EPA-ng for V3–V4 (~250 bp) on full-length ref tree (~1400 bp) |
| Bootstrap labels restored via SEPP rename script | Auto-generated per run |

---

## Path Variables

```bash
BASE="/bioinf/home/lgallucc/04_DualPlume/03_16S/Sulfurimonas_v3Aurora"
SILVA_DB="/bioinf/home/lgallucc/databases/silva_db/138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
SILVA_ALIGNED="/bioinf/home/lgallucc/databases/silva_db/138.2/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
SILVA_ARB="/bioinf/home/lgallucc/databases/silva_db/138.2/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb"
QUALITY="/bioinf/home/lgallucc/databases/silva_db/138.2/SILVA_138.2_SSURef_Nr99.quality"
SILVA_META="/bioinf/home/lgallucc/databases/silva_db/138.2/silva_metadata_sulfurimonas_categorized.tsv"
GENOME_DIR="/bioinf/home/lgallucc/04_DualPlume/09_GlobDB/SLFMN/SLFMNGenomes"
GLOBDB_META="$BASE/SulfurimonasGlobDBv2.tsv"

mkdir -p $BASE/{01_silva,02_genomes,03_merged,04_alignment,05_outgroup_expansion,06_phylogeny,07_shortreads/derep,07_shortreads/placement}
```

---

## Directory Structure

```
$BASE/
├── 01_silva/           SILVA extraction + outgroup
├── 02_genomes/         barrnap → SINA QC → filtered genome 16S
├── 03_merged/          merged + dereplicated sequences + metadata
├── 04_alignment/       SILVA pre-alignment + mafft --add
├── 05_outgroup_expansion/  LBA detection trees (A/B/C)
├── 06_phylogeny/       IQ-TREE3 ML tree + RAxML SEPP prep
└── 07_shortreads/
    ├── derep/          PS137 ASV input
    └── placement/      SEPP output + iTOL files
```

---

## Step 1 — Extract *Sulfurimonas* from SILVA 138.2

**Inclusion criteria:**
1. Genus = *Sulfurimonas* (penultimate SILVA taxonomy field)
2. `environment == 'Hydrothermal vent'` **OR** `Subtype` in known TypeStrain set **OR** SILVA taxonomy last field starts with a formally described species name
3. Named type strains accepted even without a metadata entry
4. Pintail ≥ 75, alignment quality ≥ 75, length ≥ 800 bp
5. Not in blacklist (`GQ349112 JN662179 JN662232 HG513655`)

**Output:** `01_silva/silva_sulfurimonas_qc.fasta`, `01_silva/pintail_report.tsv`

```python
python3 -c "
import csv
from Bio import SeqIO

SILVA   = '$SILVA_DB'
QUALITY = '$QUALITY'
META    = '$SILVA_META'
OUTDIR  = '$BASE/01_silva'

SILVA_BLACKLIST = {'GQ349112', 'JN662179', 'JN662232', 'HG513655'}

VALID_SPECIES_PREFIXES = {
    'Sulfurimonas denitrificans':  'S. denitrificans',
    'Sulfurimonas autotrophica':   'S. autotrophica',
    'Sulfurimonas gotlandica':     'S. gotlandica',
    'Sulfurimonas marisnigri':     'S. marisnigri',
    'Sulfurimonas hongkongensis':  'S. hongkongensis',
    'Sulfurimonas paralvinellae':  'S. paralvinellae',
    'Sulfurimonas indica':         'S. indica',
}

TYPESTRAIN_SUBTYPES = set(VALID_SPECIES_PREFIXES.values()) | {'TypeStrain'}

def match_named_species(species_str):
    for prefix, short in VALID_SPECIES_PREFIXES.items():
        if species_str.startswith(prefix): return short
    return None

def is_keep(acc, subtype, environment, taxonomy_str):
    if environment == 'Hydrothermal vent': return True, 'HV'
    if subtype in TYPESTRAIN_SUBTYPES:     return True, f'TypeStrain({subtype})'
    tax_fields = [x.strip() for x in taxonomy_str.strip(';').split(';') if x.strip()]
    species = tax_fields[-1] if tax_fields else ''
    short = match_named_species(species)
    if short: return True, f'NamedSpecies({short})'
    return False, 'excluded'

print('Loading quality scores...')
quality = {}
with open(QUALITY) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        acc = row['Primary Accession'].strip()
        try: quality[acc] = (float(row['Pintail Quality']), float(row['Alignment Quality']))
        except ValueError: pass
print(f'  Loaded: {len(quality)} entries')

print('Loading metadata...')
metadata = {}
with open(META) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        metadata[row['acc'].strip()] = row
print(f'  Loaded: {len(metadata)} entries')

print('Extracting Sulfurimonas...')
kept = []
removed = {'blacklist':[], 'wrong_tax':[], 'not_in_meta':[],
           'not_hv_or_strain':[], 'pintail':[], 'align_q':[], 'short':[], 'no_quality':[]}

for rec in SeqIO.parse(SILVA, 'fasta'):
    if 'Sulfurimonas' not in rec.description: continue
    parts = rec.description.split(' ', 1)
    taxonomy_str = parts[1] if len(parts) > 1 else ''
    tax_fields = [x.strip() for x in taxonomy_str.strip(';').split(';')]
    genus = tax_fields[-2] if len(tax_fields) >= 2 else tax_fields[-1] if tax_fields else ''
    if genus != 'Sulfurimonas':
        removed['wrong_tax'].append((rec.id.split('.')[0], genus)); continue
    acc = rec.id.split('.')[0]
    if acc in SILVA_BLACKLIST:
        removed['blacklist'].append(acc); continue

    # Named species allowed even without a metadata entry
    if acc not in metadata:
        tax_fields_check = [x.strip() for x in taxonomy_str.strip(';').split(';') if x.strip()]
        species_check = tax_fields_check[-1] if tax_fields_check else ''
        short = match_named_species(species_check)
        if not short:
            removed['not_in_meta'].append(acc); continue
        m = {'Subtype': short, 'environment': ''}
    else:
        m = metadata[acc]

    subtype     = m.get('Subtype','').strip()
    environment = m.get('environment','').strip()
    keep, reason = is_keep(acc, subtype, environment, taxonomy_str)
    if not keep:
        removed['not_hv_or_strain'].append((acc, environment, subtype)); continue

    if len(rec.seq) < 800: removed['short'].append(acc); continue

    # Named type strains may not be in quality file — assign perfect scores
    if acc not in quality:
        if match_named_species(tax_fields[-1] if tax_fields else ''):
            pt, aq = 100.0, 100.0
        else:
            removed['no_quality'].append(acc); continue
    else:
        pt, aq = quality[acc]

    if pt < 75: removed['pintail'].append((acc, pt)); continue
    if aq < 75: removed['align_q'].append((acc, aq)); continue

    rec.id = f'SILVA_{acc}'; rec.description = ''; kept.append(rec)

SeqIO.write(kept, f'{OUTDIR}/silva_sulfurimonas_qc.fasta', 'fasta')
print(f'  Kept:                    {len(kept)}')
print(f'  Removed blacklisted:     {len(removed[\"blacklist\"])} {removed[\"blacklist\"]}')
print(f'  Removed wrong taxonomy:  {len(removed[\"wrong_tax\"])}')
print(f'  Removed not in metadata: {len(removed[\"not_in_meta\"])}')
print(f'  Removed not HV/strain:   {len(removed[\"not_hv_or_strain\"])}')
print(f'  Removed pintail<75:      {len(removed[\"pintail\"])}')
print(f'  Removed align<75:        {len(removed[\"align_q\"])}')
print(f'  Removed short<800bp:     {len(removed[\"short\"])}')

with open(f'{OUTDIR}/pintail_report.tsv','w') as f:
    f.write('acc\tpintail\talign_quality\tenvironment\tsubtype\tstatus\n')
    for rec in kept:
        acc = rec.id.replace('SILVA_','')
        pt, aq = quality.get(acc, (100.0, 100.0))
        m_row = metadata.get(acc, {'environment':'','Subtype':''})
        f.write(f'{acc}\t{pt}\t{aq}\t{m_row.get(\"environment\",\"\")}\t{m_row.get(\"Subtype\",\"\")}\tKEPT\n')
    for acc in removed['not_in_meta']:
        f.write(f'{acc}\t\t\t\t\tREMOVED_NOT_IN_META\n')
    for acc, pt in removed['pintail']:
        f.write(f'{acc}\t{pt}\t\t\t\tREMOVED_PINTAIL\n')
    for acc, aq in removed['align_q']:
        f.write(f'{acc}\t\t{aq}\t\t\tREMOVED_ALIGN\n')
    for acc, env, sub in removed['not_hv_or_strain']:
        f.write(f'{acc}\t\t\t{env}\t{sub}\tREMOVED_NOT_HV_OR_STRAIN\n')
print(f'Report: {OUTDIR}/pintail_report.tsv')
print('Step 1 done.')
"
```

---

## Step 1b — Extended Metadata for Selected SILVA Sequences

Combines manual metadata, SILVA full metadata, and Pintail scores. Flags each sequence as TypeStrain, Named species, or Environmental HV. Used for review and publication.

**Output:** `01_silva/silva_selection_metadata.csv`

```python
python3 -c "
import csv
from Bio import SeqIO
from collections import Counter

SILVA_QC   = '$BASE/01_silva/silva_sulfurimonas_qc.fasta'
SILVA_META = '$SILVA_META'
SILVA_FULL = '/bioinf/home/lgallucc/databases/silva_db/138.2/SILVA_138.2_SSURef_Nr99.full_metadata'
QUALITY    = '$QUALITY'
SILVA_DB   = '$SILVA_DB'
OUT        = '$BASE/01_silva/silva_selection_metadata.csv'

VALID_SPECIES_PREFIXES = {
    'Sulfurimonas denitrificans':  'S. denitrificans',
    'Sulfurimonas autotrophica':   'S. autotrophica',
    'Sulfurimonas gotlandica':     'S. gotlandica',
    'Sulfurimonas marisnigri':     'S. marisnigri',
    'Sulfurimonas hongkongensis':  'S. hongkongensis',
    'Sulfurimonas paralvinellae':  'S. paralvinellae',
    'Sulfurimonas indica':         'S. indica',
}
TYPESTRAIN_SUBTYPES = set(VALID_SPECIES_PREFIXES.values()) | {'TypeStrain'}

def match_species(sp):
    for prefix, short in VALID_SPECIES_PREFIXES.items():
        if sp.startswith(prefix): return short
    return None

meta_old = {}
with open(SILVA_META) as f:
    for row in csv.DictReader(f, delimiter='\t'): meta_old[row['acc'].strip()] = row

meta_full = {}
with open(SILVA_FULL, errors='replace') as f:
    for row in csv.DictReader(f, delimiter='\t'): meta_full[row['acc'].strip()] = row

quality = {}
with open(QUALITY) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        acc = row['Primary Accession'].strip()
        try: quality[acc] = (float(row['Pintail Quality']), float(row['Alignment Quality']))
        except ValueError: pass

selected_accs = set()
for rec in SeqIO.parse(SILVA_QC, 'fasta'):
    selected_accs.add(rec.id.replace('SILVA_',''))

tax_map = {}
for rec in SeqIO.parse(SILVA_DB, 'fasta'):
    acc = rec.id.split('.')[0]
    if acc not in selected_accs: continue
    parts = rec.description.split(' ', 1)
    if acc not in tax_map: tax_map[acc] = parts[1] if len(parts) > 1 else ''

rows_out = []
for acc in selected_accs:
    old = meta_old.get(acc, {}); full = meta_full.get(acc, {})
    pt, aq = quality.get(acc, (100.0, 100.0))
    tax_str = tax_map.get(acc, '')
    tax_fields = [x.strip() for x in tax_str.strip(';').split(';') if x.strip()]
    species = tax_fields[-1] if tax_fields else ''
    subtype = old.get('Subtype','').strip()
    env     = old.get('environment','').strip()
    sp_short = match_species(species)
    if subtype in TYPESTRAIN_SUBTYPES: flag = f'TypeStrain ({subtype})'
    elif sp_short:                     flag = f'Named species ({sp_short})'
    else:                              flag = 'Environmental HV'
    rows_out.append({'acc':acc,'silva_id':f'SILVA_{acc}',
        'species_silva':species,'species_short':sp_short or '',
        'subtype_manual':subtype,'environment':env,'flag':flag,
        'isolation_source':full.get('isolation_source',''),
        'strain':full.get('strain',''),
        'culture_collection':full.get('culture_collection',''),
        'lat_lon':full.get('lat_lon',''),
        'depth_slv':full.get('depth_slv',''),
        'pintail':pt,'align_quality':aq,'taxonomy_full':tax_str.strip()})

rows_out.sort(key=lambda x: (x['flag'], x['acc']))
with open(OUT,'w',newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
    writer.writeheader(); writer.writerows(rows_out)
flags = Counter(r['flag'] for r in rows_out)
print('Selection summary:')
for k,v in sorted(flags.items(), key=lambda x:-x[1]): print(f'  {k:<45} {v}')
print(f'Total: {len(rows_out)} | Output: {OUT}')
"
```

---

## Step 1c — Extract *Sulfuricurvum* Outgroup from SILVA 138.2

Target: ~10 sequences. Specific accessions prioritised (`AB053951`, `CP002355`, `MIBO01000018`); remainder selected by highest Pintail.

**QC:** Pintail ≥ 75, alignment quality ≥ 75, length ≥ 1000 bp.

**Output:** `01_silva/outgroups.fasta`

```python
python3 -c "
import csv
from Bio import SeqIO

SILVA   = '$SILVA_DB'
QUALITY = '$QUALITY'
OUTDIR  = '$BASE/01_silva'

quality = {}
with open(QUALITY) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        acc = row['Primary Accession'].strip()
        try: quality[acc] = (float(row['Pintail Quality']), float(row['Alignment Quality']))
        except ValueError: pass

SPECIFIC = {'AB053951', 'CP002355', 'MIBO01000018'}
N_TARGET = 10
hits = []; specific = {}

for rec in SeqIO.parse(SILVA, 'fasta'):
    if 'Sulfuricurvum' not in rec.description: continue
    parts = rec.description.split(' ', 1)
    taxonomy = parts[1] if len(parts) > 1 else ''
    tax_fields = [x.strip() for x in taxonomy.strip(';').split(';')]
    rec_genus = tax_fields[-2] if len(tax_fields) >= 2 else tax_fields[-1] if tax_fields else ''
    if rec_genus != 'Sulfuricurvum': continue
    acc = rec.id.split('.')[0]
    if len(rec.seq) < 1000: continue
    pt, aq = quality.get(acc, (0, 0))
    if pt < 75 or aq < 75: continue
    if acc in SPECIFIC: specific[acc] = (rec, pt)
    else: hits.append((rec, pt, acc))

all_og = []
print('Sulfuricurvum outgroup:')
for acc, (rec, pt) in specific.items():
    rec.id = f'OUTGROUP_Sulfuricurvum_{acc}'; rec.description = ''
    all_og.append(rec)
    print(f'  [SPECIFIC] {acc}  pintail={pt:.0f}  len={len(rec.seq)}')
missing = SPECIFIC - set(specific.keys())
if missing: print(f'  WARNING: not in SILVA NR99: {missing}')
hits.sort(key=lambda x: -x[1])
added = 0
for rec, pt, acc in hits:
    if added >= N_TARGET - len(specific): break
    rec.id = f'OUTGROUP_Sulfuricurvum_{acc}'; rec.description = ''
    all_og.append(rec)
    print(f'  [AUTO]     {acc}  pintail={pt:.0f}  len={len(rec.seq)}')
    added += 1

SeqIO.write(all_og, f'{OUTDIR}/outgroups.fasta', 'fasta')
print(f'Total outgroups: {len(all_og)}')
"
```

---

## Step 2a — Extract 16S rRNA from Genomes (barrnap ≥ 1.0)

barrnap ≥ 1.0 uses Infernal (Covariance Models). GFF on stdout is filtered for `16S_rRNA`; `bedtools getfasta` extracts sequences. Only the longest 16S copy per genome is retained. Sequences < 800 bp are excluded.

**Output:** `02_genomes/genome_16S_all.fasta`

```bash
for fa in $GENOME_DIR/*.fa $GENOME_DIR/*.fna $GENOME_DIR/*.fasta; do
    [ -f "$fa" ] || continue
    base=$(basename $fa | sed 's/\.[^.]*$//')
    barrnap --kingdom bac --threads 4 --quiet "$fa" \
        > $BASE/02_genomes/${base}.gff
    grep "16S_rRNA" $BASE/02_genomes/${base}.gff \
        | bedtools getfasta -fi "$fa" -bed - -s -name \
        > $BASE/02_genomes/${base}_16S_raw.fasta
done

python3 -c "
from Bio import SeqIO
import os, glob
OUTDIR = '$BASE/02_genomes'
kept = []; skipped = []
for f in sorted(glob.glob(f'{OUTDIR}/*_16S_raw.fasta')):
    base = os.path.basename(f).replace('_16S_raw.fasta','')
    seqs = list(SeqIO.parse(f, 'fasta'))
    if not seqs: skipped.append(base); continue
    best = max(seqs, key=lambda r: len(r.seq))
    best.id = f'GENOME_{base}'; best.description = ''
    if len(best.seq) < 800: skipped.append(f'{base} (len={len(best.seq)})'); continue
    kept.append(best)
SeqIO.write(kept, f'{OUTDIR}/genome_16S_all.fasta', 'fasta')
print(f'Full-length genomes (reference tree): {len(kept)}')
print(f'Skipped (no 16S or too short):        {len(skipped)}')
"
```

---

## Step 2b — SINA Quality Check on Genome 16S Sequences

> **⚠ PAUSE:** Review `sina_qc_report.tsv` before running Step 2c.

Thresholds: identity ≥ 0.92, alignment quality ≥ 75, taxonomy must include *Sulfurimonas*. Environment filter: hydrothermal vent or isolate/type strain only.

**Output:** `02_genomes/genome_16S_sina.fasta`, `02_genomes/sina_qc_report.tsv`

```bash
sina \
    -i $BASE/02_genomes/genome_16S_all.fasta \
    -r $SILVA_ARB \
    -o $BASE/02_genomes/genome_16S_sina.fasta \
    --meta-fmt header \
    --search --search-db $SILVA_ARB \
    --search-max-result 1 \
    --lca-fields tax_slv \
    --search-min-sim 0.7 \
    -p 8 \
    2> $BASE/02_genomes/sina.log

echo "SINA done: $(grep -c '^>' $BASE/02_genomes/genome_16S_sina.fasta) sequences"
```

```python
python3 -c "
import re, csv
from Bio import SeqIO

FASTA      = '$BASE/02_genomes/genome_16S_sina.fasta'
GLOBDB_META= '$GLOBDB_META'
OUTDIR     = '$BASE/02_genomes'
BLACKLIST  = {'GENOME_BCRBG_25302', 'GENOME_GCA_015487435'}
MIN_IDENT  = 0.92
MIN_ALIGN  = 75

globdb = {}
with open(GLOBDB_META) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gid = row['ID'].strip()
        globdb[gid] = row
        globdb[re.sub(r'\.\d+$','',gid)] = row

def get_globdb(label):
    gid = label.replace('GENOME_','')
    return globdb.get(gid, globdb.get(re.sub(r'\.\d+$','',gid), {}))

def is_hv_or_isolate(label):
    r = get_globdb(label)
    env = r.get('EnvType','').strip()
    sub = r.get('Subtype','').strip()
    if env == 'Hydrothermal vent': return True, 'HV'
    kw = {'isolate','culture','type strain','typestrain','pure culture','laboratory'}
    if any(k in env.lower() for k in kw): return True, f'Isolate({env})'
    if any(k in sub.lower() for k in kw): return True, f'Isolate({sub})'
    if not env and r: return True, 'Isolate(no_env)'
    return False, f'excluded({env})'

def parse_field(desc, field):
    m = re.search(rf'\[{field}=([^\]]+)\]', desc)
    return m.group(1).strip() if m else ''

def parse_nearest(s):
    if not s: return 0.0
    try: return float(s.strip().split()[0].split('~')[1])
    except: return 0.0

def is_sulfurimonas(tax):
    return 'Sulfurimonas' in [x.strip() for x in tax.strip(';').split(';')] if tax else False

def get_deepest(tax):
    fields = [x.strip() for x in tax.strip(';').split(';') if x.strip()]
    return fields[-1] if fields else ''

rows = []
for rec in SeqIO.parse(FASTA, 'fasta'):
    name = rec.id; desc = rec.description
    aln_q = float(parse_field(desc, 'align_quality_slv') or 0)
    taxon = parse_field(desc, 'lca_tax_slv')
    ident = parse_nearest(parse_field(desc, 'nearest_slv'))
    hv_ok, hv_reason = is_hv_or_isolate(name)
    rows.append((name, aln_q, ident, get_deepest(taxon), taxon, hv_ok, hv_reason))

rows.sort(key=lambda x: x[2])
print('=' * 90)
print('GENOME 16S QC — SINA 1.7.2')
print('=' * 90)
print(f'{\"Genome\":<35} {\"AlnQ\":>6} {\"Ident\":>7} {\"Deepest tax\":<25} Flag')
print('-' * 90)

passed = []; failed = []; blacklisted = []
for name, aln_q, ident, deepest, taxon, hv_ok, hv_reason in rows:
    flags = []
    if name in BLACKLIST:          flags.append('BLACKLISTED')
    if aln_q < MIN_ALIGN:          flags.append(f'LOW_ALN<{MIN_ALIGN}')
    if ident < MIN_IDENT:          flags.append(f'LOW_ID<{MIN_IDENT}')
    if not is_sulfurimonas(taxon): flags.append(f'NOT_SLFMN({deepest})')
    if not hv_ok:                  flags.append(hv_reason)
    flag_str = ' | '.join(flags) if flags else f'OK ({hv_reason})'
    print(f'{name:<35} {aln_q:>6.1f} {ident:>7.3f} {deepest:<25} {flag_str}')
    if 'BLACKLISTED' in flag_str: blacklisted.append(name)
    elif flags: failed.append((name, flags))
    else: passed.append(name)

print('-' * 90)
print(f'PASSED={len(passed)} | FAILED={len(failed)} | BLACKLISTED={len(blacklisted)}')

with open(f'{OUTDIR}/sina_qc_report.tsv','w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['genome_id','align_quality','identity_nearest','deepest_tax',
                     'full_taxonomy','hv_or_isolate','flags','decision'])
    for name, aln_q, ident, deepest, taxon, hv_ok, hv_reason in rows:
        flags = []
        if name in BLACKLIST:          flags.append('BLACKLISTED')
        if aln_q < MIN_ALIGN:          flags.append('LOW_ALN_QUAL')
        if ident < MIN_IDENT:          flags.append('LOW_IDENT')
        if not is_sulfurimonas(taxon): flags.append('NOT_SULFURIMONAS')
        if not hv_ok:                  flags.append('NOT_HV_OR_ISOLATE')
        decision = 'REMOVE' if flags else 'KEEP'
        writer.writerow([name, aln_q, ident, deepest, taxon, hv_reason,
                         '|'.join(flags), decision])
print('>>> PAUSE: review sina_qc_report.tsv, override if needed, then run Step 2c')
"
```

---

## Step 2c — Apply SINA QC Decisions

> **⚠ PAUSE BEFORE THIS STEP — review `sina_qc_report.tsv` first.**

Reads `decision` column; removes any row with `REMOVE`. Allows manual overrides by editing the TSV.

**Output:** `02_genomes/genome_16S.fasta`

```python
python3 -c "
import csv
from Bio import SeqIO
FASTA  = '$BASE/02_genomes/genome_16S_all.fasta'
REPORT = '$BASE/02_genomes/sina_qc_report.tsv'
OUT    = '$BASE/02_genomes/genome_16S.fasta'
remove = set()
with open(REPORT) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['decision'].strip() == 'REMOVE': remove.add(row['genome_id'].strip())
print(f'Removing: {sorted(remove)}')
kept = [rec for rec in SeqIO.parse(FASTA,'fasta') if rec.id not in remove]
SeqIO.write(kept, OUT, 'fasta')
print(f'Genome 16S kept: {len(kept)}')
"
```

---

## Step 3 — Merge + CD-HIT Dereplication at 98.6%

Merges SILVA QC-passed (includes named type strains) + genomes + outgroups. CD-HIT at 98.6% collapses identical sequences including multiple 16S copies from the same genome (Kim et al. 2014).

**Output:** `03_merged/all_sequences_derep986.fasta`

```bash
cat $BASE/01_silva/silva_sulfurimonas_qc.fasta \
    $BASE/02_genomes/genome_16S.fasta \
    $BASE/01_silva/outgroups.fasta \
    > $BASE/03_merged/all_sequences_raw.fasta

echo "Total before CD-HIT: $(grep -c '^>' $BASE/03_merged/all_sequences_raw.fasta)"
echo "Outgroups:"; grep "^>OUTGROUP" $BASE/03_merged/all_sequences_raw.fasta | sed 's/^>/  /'

cd-hit-est \
    -i $BASE/03_merged/all_sequences_raw.fasta \
    -o $BASE/03_merged/all_sequences_derep986.fasta \
    -c 0.986 -n 8 -T 8 -M 8000 -d 0

echo "After CD-HIT 98.6%: $(grep -c '^>' $BASE/03_merged/all_sequences_derep986.fasta)"
echo "Type strains surviving:"
grep "^>SILVA_CP\|^>SILVA_AFRZ\|^>SILVA_AB\|^>SILVA_L40808" \
    $BASE/03_merged/all_sequences_derep986.fasta | sed 's/^>/  /'
echo "Outgroups surviving:"; grep "^>OUTGROUP" $BASE/03_merged/all_sequences_derep986.fasta | sed 's/^>/  /'
```

---

## Step 3b — Unified Metadata for Manual Review

> **⚠ PAUSE:** Fill `REMOVE=1` for any sequences to exclude manually, then proceed to Step 4.

Combines SILVA old/new metadata + GlobDB for genomes. Produces a single TSV for editorial review.

**Output:** `03_merged/unified_metadata_review.tsv`

```python
python3 -c "
import csv, re
from Bio import SeqIO

DEREP          = '$BASE/03_merged/all_sequences_derep986.fasta'
SILVA_META     = '$SILVA_META'
SILVA_META_NEW = '/bioinf/home/lgallucc/databases/silva_db/138.2/SILVA_138.2_SSURef_Nr99.full_metadata'
GLOBDB_META    = '$GLOBDB_META'
SILVA_SEL_META = '$BASE/01_silva/silva_selection_metadata.csv'
OUT            = '$BASE/03_merged/unified_metadata_review.tsv'

silva_old = {}
with open(SILVA_META) as f:
    for row in csv.DictReader(f, delimiter='\t'): silva_old[row['acc'].strip()] = row

silva_sel = {}
with open(SILVA_SEL_META) as f:
    for row in csv.DictReader(f): silva_sel[row['acc'].strip()] = row

silva_new = {}
with open(SILVA_META_NEW, errors='replace') as f:
    for row in csv.DictReader(f, delimiter='\t'): silva_new[row['acc'].strip()] = row

globdb = {}
with open(GLOBDB_META) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gid = row['ID'].strip(); globdb[gid] = row; globdb[re.sub(r'\.\d+$','',gid)] = row

kept_labels = []
with open(OUT,'w',newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['id','source','acc','subtype_manual','flag','envtype_globdb',
                     'environment_silva','isolation_source','habitat_slv',
                     'depth_slv','lat_lon','REMOVE','NOTES'])
    for rec in SeqIO.parse(DEREP, 'fasta'):
        label = rec.id; kept_labels.append(label)
        if label.startswith('OUTGROUP_'):
            writer.writerow([label,'OUTGROUP',label,'Outgroup','Outgroup','','','','','','','',''])
        elif label.startswith('SILVA_'):
            acc = label.replace('SILVA_','').split('_')[0]
            old = silva_old.get(acc,{}); new = silva_new.get(acc,{}); sel = silva_sel.get(acc,{})
            writer.writerow([label,'SILVA',acc,
                old.get('Subtype',''),sel.get('flag',''),'',old.get('environment',''),
                new.get('isolation_source',''),new.get('habitat_slv',''),
                new.get('depth_slv',''),new.get('lat_lon',''),'',''])
        elif label.startswith('GENOME_'):
            gid = label.replace('GENOME_','')
            g = globdb.get(gid, globdb.get(re.sub(r'\.\d+$','',gid),{}))
            writer.writerow([label,'GENOME',gid,
                g.get('Subtype',''),'',g.get('EnvType',''),
                '','','','','','',''])

print(f'Total: {len(kept_labels)} | SILVA: {sum(1 for l in kept_labels if l.startswith(\"SILVA_\"))} | GENOME: {sum(1 for l in kept_labels if l.startswith(\"GENOME_\"))} | OUTGROUP: {sum(1 for l in kept_labels if l.startswith(\"OUTGROUP_\"))}')
print('>>> PAUSE: review unified_metadata_review.tsv, fill REMOVE=1 if needed')
"
```

---

## Step 4 — Assign Final Habitat Categories + Apply Manual Exclusions

Reads `REMOVE` column; assigns final habitat categories based on `subtype_manual` and environment keywords.

**Output:** `03_merged/all_sequences_clean.fasta`, `03_merged/unified_metadata_clean.tsv`

**Habitat categories:**

| Category | Criteria |
|----------|---------|
| Outgroup | source == OUTGROUP |
| HV inactive chimney | subtype contains inactive/dead/extinct |
| HV active chimney | subtype contains chimney/active/bubbling pool |
| HV fluids | subtype contains fluids/fluid/diffuse |
| HV pelagic | subtype contains plume/background water |
| HV artificial | subtype contains artificial/colonizer |
| HV benthic | hydrothermal vent, no subtype match |
| TypeStrain | all other |

```python
python3 -c "
import csv
from Bio import SeqIO
from collections import Counter

META     = '$BASE/03_merged/unified_metadata_review.tsv'
DEREP    = '$BASE/03_merged/all_sequences_derep986.fasta'
OUT      = '$BASE/03_merged/all_sequences_clean.fasta'
OUT_META = '$BASE/03_merged/unified_metadata_clean.tsv'

def get_category(row):
    src = row.get('source','')
    sub = row.get('subtype_manual','').strip().lower()
    env = row.get('environment_silva','').strip() if src=='SILVA' else row.get('envtype_globdb','').strip()
    if src == 'OUTGROUP': return 'Outgroup'
    if env == 'Hydrothermal vent' or 'hv' in sub or 'hydrothermal' in sub:
        if any(x in sub for x in ['inactive','dead','extinct']): return 'HV inactive chimney'
        if any(x in sub for x in ['chimney','active','bubbling pool']): return 'HV active chimney'
        if any(x in sub for x in ['fluids','fluid','diffuse']):         return 'HV fluids'
        if any(x in sub for x in ['plume','background water']):         return 'HV pelagic'
        if any(x in sub for x in ['artificial','colonizer']):           return 'HV artificial'
        return 'HV benthic'
    return 'TypeStrain'

meta = {}
with open(META) as f:
    for row in csv.DictReader(f, delimiter='\t'): meta[row['id'].strip()] = row

remove_ids = {l for l, r in meta.items() if r.get('REMOVE','').strip() == '1'}
print(f'Manual removals: {len(remove_ids)}')

kept = [rec for rec in SeqIO.parse(DEREP,'fasta') if rec.id not in remove_ids]
SeqIO.write(kept, OUT, 'fasta')
print(f'Kept: {len(kept)}')

with open(OUT_META,'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['id','source','acc','final_category','subtype_manual',
                     'flag','envtype_globdb','environment_silva',
                     'isolation_source','habitat_slv','depth_slv','lat_lon','NOTES'])
    for rec in kept:
        row = meta.get(rec.id,{})
        writer.writerow([rec.id,row.get('source',''),row.get('acc',''),
            get_category(row),row.get('subtype_manual',''),row.get('flag',''),
            row.get('envtype_globdb',''),row.get('environment_silva',''),
            row.get('isolation_source',''),row.get('habitat_slv',''),
            row.get('depth_slv',''),row.get('lat_lon',''),row.get('NOTES','')])

counts = Counter(get_category(meta.get(rec.id,{})) for rec in kept)
print('\nCategory distribution:')
for k,v in sorted(counts.items(), key=lambda x:-x[1]): print(f'  {k:<30} {v}')
"
```

---

## Step 5 — Structure-Guided Alignment

**Strategy:** SILVA `full_align_trunc` pre-alignment as base profile → `mafft --add` for genomes → `mafft --add` for *Sulfuricurvum* outgroup. Adding outgroup last prevents it from distorting the ingroup profile.

### 5a — Extract Pre-Aligned *Sulfurimonas* from SILVA `full_align_trunc`

**Output:** `04_alignment/sulfurimonas_silva_prealigned.fasta`

```python
python3 -c "
import csv, gzip
from Bio import SeqIO
SILVA_ALIGNED = '$SILVA_ALIGNED'
META_CLEAN    = '$BASE/03_merged/unified_metadata_clean.tsv'
OUTFILE       = '$BASE/04_alignment/sulfurimonas_silva_prealigned.fasta'
wanted = set()
with open(META_CLEAN) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['source'] == 'SILVA': wanted.add(row['acc'].strip())
print(f'Wanted SILVA accessions: {len(wanted)}')
kept = []
open_fn = gzip.open if SILVA_ALIGNED.endswith('.gz') else open
with open_fn(SILVA_ALIGNED,'rt') as f:
    for rec in SeqIO.parse(f,'fasta'):
        acc = rec.id.split('.')[0]
        if acc in wanted:
            rec.id = f'SILVA_{acc}'; rec.description = ''; kept.append(rec)
SeqIO.write(kept, OUTFILE, 'fasta')
print(f'Extracted: {len(kept)} | Raw alignment length: {len(kept[0].seq) if kept else 0}')
"
```

### 5b — Remove All-Gap Columns + Convert `.` to `-`

The SILVA alignment uses `.` for gaps; MAFFT requires `-`.

**Output:** `04_alignment/sulfurimonas_silva_degapped.fasta`

```python
python3 -c "
import numpy as np
from Bio import SeqIO
INFILE  = '$BASE/04_alignment/sulfurimonas_silva_prealigned.fasta'
OUTFILE = '$BASE/04_alignment/sulfurimonas_silva_degapped.fasta'
recs = list(SeqIO.parse(INFILE,'fasta'))
print(f'Sequences: {len(recs)} | Raw length: {len(recs[0].seq)}')
arr = np.array([list(str(r.seq).upper()) for r in recs])
not_all_gap = np.any((arr != '-') & (arr != '.'), axis=0)
print(f'Informative columns: {not_all_gap.sum()}')
for i, rec in enumerate(recs):
    rec.seq = rec.seq.__class__(''.join(arr[i,not_all_gap]).replace('.', '-'))
SeqIO.write(recs, OUTFILE, 'fasta')
print(f'Degapped alignment: {len(recs[0].seq)} columns')
"
```

### 5c — Add Genomes via `mafft --add`

**Output:** `04_alignment/sulfurimonas_with_genomes.fasta`

```bash
seqkit grep -r -p "^GENOME_" $BASE/03_merged/all_sequences_clean.fasta \
    > $BASE/04_alignment/genomes_only.fasta
echo "Genomes to add: $(grep -c '^>' $BASE/04_alignment/genomes_only.fasta)"

mafft --thread 8 \
    --add $BASE/04_alignment/genomes_only.fasta \
    $BASE/04_alignment/sulfurimonas_silva_degapped.fasta \
    > $BASE/04_alignment/sulfurimonas_with_genomes.fasta \
    2> $BASE/04_alignment/mafft_add_genomes.log
echo "Exit: $? | With genomes: $(grep -c '^>' $BASE/04_alignment/sulfurimonas_with_genomes.fasta)"
```

### 5d — Add *Sulfuricurvum* Outgroup via `mafft --add`

**Output:** `04_alignment/sulfurimonas_trimmed.fasta`

```bash
seqkit grep -r -p "^OUTGROUP_Sulfuricurvum" $BASE/01_silva/outgroups.fasta \
    > $BASE/04_alignment/sulfuricurvum_only.fasta

mafft --thread 8 \
    --add $BASE/04_alignment/sulfuricurvum_only.fasta \
    $BASE/04_alignment/sulfurimonas_with_genomes.fasta \
    > $BASE/04_alignment/sulfurimonas_trimmed.fasta \
    2> $BASE/04_alignment/mafft_add_sulfuricurvum.log
echo "Exit: $? | Final: $(grep -c '^>' $BASE/04_alignment/sulfurimonas_trimmed.fasta)"
echo "Outgroups:"; grep "^>OUTGROUP" $BASE/04_alignment/sulfurimonas_trimmed.fasta | sed 's/^>/  /'
```

**Verify outgroup coverage (expect > 85% non-gap):**

```python
python3 -c "
from Bio import SeqIO
for rec in SeqIO.parse('$BASE/04_alignment/sulfurimonas_trimmed.fasta','fasta'):
    if 'OUTGROUP' in rec.id:
        seq = str(rec.seq); gaps = seq.count('-') + seq.count('.')
        total = len(seq)
        print(f'  {rec.id}: {total-gaps}/{total} bp ({100*(total-gaps)/total:.1f}% non-gap)')
"
```

---

## Step 6 — Outgroup Expansion RF Test (LBA Detection)

> **⚠ PAUSE:** If B vs C RF > 10%, remove flagged sequences and rerun from Step 5.

Three tree comparisons:
- **A** — ingroup only (no outgroup)
- **B** — ingroup + *Sulfuricurvum* (with `-o` flag)
- **C** — full alignment, no `-o` flag (naive rooting)

RF(B, C) near 0% → no LBA. RF(B, C) > 10% → long-branch attraction suspected from outgroup placement.

```bash
TRIMMED="$BASE/04_alignment/sulfurimonas_trimmed.fasta"
OG_DIR="$BASE/05_outgroup_expansion"
OG_SULFURICURVUM=$(grep "^>OUTGROUP_Sulfuricurvum" $TRIMMED | sed 's/^>//' | tr -d '\r' | paste -sd ',' -)
echo "Outgroup: $OG_SULFURICURVUM"

seqkit grep -v -r -p "^OUTGROUP_" $TRIMMED > $OG_DIR/stageA.fasta
cp $TRIMMED $OG_DIR/stageB.fasta
cp $TRIMMED $OG_DIR/stageC.fasta

iqtree3 -s $OG_DIR/stageA.fasta -m MFP --prefix $OG_DIR/stageA_tree \
    -B 1000 -nm 1000 --alrt 1000 -nt 8 --redo --quiet
iqtree3 -s $OG_DIR/stageB.fasta -m MFP --prefix $OG_DIR/stageB_tree \
    -B 1000 -nm 1000 --alrt 1000 -nt 8 -o "$OG_SULFURICURVUM" --redo --quiet
iqtree3 -s $OG_DIR/stageC.fasta -m MFP --prefix $OG_DIR/stageC_tree \
    -B 1000 -nm 1000 --alrt 1000 -nt 8 --redo --quiet
```

```python
python3 -c "
import dendropy
OUTDIR = '$OG_DIR'
stages = {'A':'stageA_tree','B':'stageB_tree','C':'stageC_tree'}
taxon_namespace = dendropy.TaxonNamespace()
trees = {}
for name, prefix in stages.items():
    try:
        t = dendropy.Tree.get(path=f'{OUTDIR}/{prefix}.treefile', schema='newick',
                              taxon_namespace=taxon_namespace)
        og = [tx for tx in t.taxon_namespace if 'OUTGROUP' in (tx.label or '')]
        t.prune_taxa(og); t.is_rooted = False; trees[name] = t
    except Exception as e: print(f'Could not load {name}: {e}')
labels = list(trees.keys())
n_taxa = len(taxon_namespace)
max_rf = 2*(n_taxa-3) if n_taxa > 3 else None
print('\n=== RF distances (ingroup only) ===')
print(f'{\"\":6}', end='')
for n in labels: print(f'{n:>8}', end='')
print()
for n1 in labels:
    print(f'{n1:<6}', end='')
    for n2 in labels:
        try: rf = dendropy.calculate.treecompare.symmetric_difference(trees[n1],trees[n2]); print(f'{rf:>8}', end='')
        except: print(f'{\"NA\":>8}', end='')
    print()
if max_rf:
    print(f'Max possible RF: {max_rf}')
    for (n1,n2), label in [(('B','C'),'LBA test'),(('A','B'),'Ingroup vs +Sulfuricurvum')]:
        if n1 in trees and n2 in trees:
            rf = dendropy.calculate.treecompare.symmetric_difference(trees[n1],trees[n2])
            flag = ' <- LBA suspected' if rf/max_rf > 0.10 else ''
            print(f'  {n1} vs {n2} ({label}): RF={rf} ({rf/max_rf*100:.1f}%){flag}')
"
```

> **If LBA detected:**
> ```bash
> seqkit grep -v -n -p "SEQ_ID" $TRIMMED > $BASE/04_alignment/sulfurimonas_trimmed_clean.fasta
> # Then rerun Step 7 with sulfurimonas_trimmed_clean.fasta
> ```

---

## Step 7 — Final ML Tree with IQ-TREE3

> **⚠ PAUSE:** Review `flagged_long_branches.txt` before Step 8.

Deduplicates identical accessions (multiple 16S copies per genome) — keeps longest copy. Runs IQ-TREE3 with MFP model selection + UFBoot 1000 + SH-aLRT 1000. Flags long branches at mean + 3 SD.

**Output:** `06_phylogeny/sulfurimonas_final.treefile`, `06_phylogeny/flagged_long_branches.txt`

```bash
TRIMMED="$BASE/04_alignment/sulfurimonas_trimmed.fasta"
OG_SULFURICURVUM=$(grep "^>OUTGROUP_Sulfuricurvum" $TRIMMED | sed 's/^>//' | tr -d '\r' | paste -sd ',' -)
```

```python
# Deduplicate: keep longest copy per accession
python3 -c "
from Bio import SeqIO
from collections import defaultdict
INFILE  = '$BASE/04_alignment/sulfurimonas_trimmed.fasta'
OUTFILE = '$BASE/04_alignment/sulfurimonas_trimmed_deduped.fasta'
seen = defaultdict(list)
for rec in SeqIO.parse(INFILE,'fasta'):
    seen[rec.id].append(rec)
kept = []; removed = 0
for seq_id, recs in seen.items():
    if len(recs) == 1: kept.append(recs[0])
    else:
        best = max(recs, key=lambda r: len(str(r.seq).replace('-','').replace('.','')))
        kept.append(best); removed += len(recs)-1
        print(f'  Deduped {seq_id}: kept 1 of {len(recs)} copies')
SeqIO.write(kept, OUTFILE, 'fasta')
print(f'Kept: {len(kept)} | Removed duplicates: {removed}')
"
```

```bash
TRIMMED="$BASE/04_alignment/sulfurimonas_trimmed_deduped.fasta"
OG_SULFURICURVUM=$(grep "^>OUTGROUP_Sulfuricurvum" $TRIMMED | sed 's/^>//' | tr -d '\r' | paste -sd ',' -)

iqtree3 -s $TRIMMED -m MFP \
    --prefix $BASE/06_phylogeny/sulfurimonas_final \
    -B 1000 --bnni --alrt 1000 -nt 8 \
    -o "$OG_SULFURICURVUM" --redo --quiet

grep "Best-fit model" $BASE/06_phylogeny/sulfurimonas_final.iqtree
```

```python
# Flag long branches at mean + 3 SD
python3 -c "
import dendropy, statistics
t = dendropy.Tree.get(path='$BASE/06_phylogeny/sulfurimonas_final.treefile', schema='newick')
branches = [(nd.taxon.label, nd.edge.length) for nd in t.leaf_node_iter() if nd.edge.length and nd.taxon]
lengths = [b[1] for b in branches]
mean_l = statistics.mean(lengths); sd_l = statistics.stdev(lengths); thresh3 = mean_l + 3*sd_l
print(f'Mean: {mean_l:.5f} | SD: {sd_l:.5f} | Threshold 3SD: {thresh3:.5f}')
flagged = [(l,v) for l,v in sorted(branches, key=lambda x:-x[1]) if v > thresh3 and 'OUTGROUP' not in l]
print(f'Flagged at 3SD (non-outgroup): {len(flagged)}')
for l,v in flagged: print(f'  {l}\t{v:.6f}')
print('\nOutgroup branches:')
for l,v in sorted([(l,v) for l,v in branches if 'OUTGROUP' in l], key=lambda x:-x[1]):
    print(f'  {l}\t{v:.6f}')
with open('$BASE/06_phylogeny/flagged_long_branches.txt','w') as f:
    [f.write(f'{l}\t{v:.6f}\n') for l,v in flagged]
"
```

> **If long branches require removal:**
> ```bash
> seqkit grep -v -n -p "SEQ_ID" $TRIMMED > $BASE/04_alignment/sulfurimonas_trimmed_final.fasta
> # Rerun Step 7 with sulfurimonas_trimmed_final.fasta
> ```

---

## Step 8 — SEPP Phylogenetic Placement of PS137 ASVs

> **Requires:** `conda activate SEPP`
> **pplacer alpha17+ must be installed in the SEPP env.** Replace bundled binary if needed:
> `cp $(which pplacer) ~/.sepp/bundled-v*/pplacer`

**Strategy:** SEPP uses ensemble HMM decomposition for placing short V3–V4 reads (~250 bp) onto a full-length reference tree (~1400 bp). Better than EPA-ng for this fragment/reference length ratio.

Bootstrap labels are restored via the auto-generated rename script produced by SEPP.

### 8a — U→T Conversion

pplacer requires DNA (not RNA) input.

```bash
sed '/^>/! s/U/T/g; /^>/! s/u/t/g' \
    $BASE/04_alignment/sulfurimonas_trimmed_deduped.fasta \
    > $BASE/04_alignment/sulfurimonas_trimmed_dna.fasta

ALIGN_DNA="$BASE/04_alignment/sulfurimonas_trimmed_dna.fasta"
ALIGN_RED=$ALIGN_DNA; [ -f "${ALIGN_DNA}.reduced" ] && ALIGN_RED="${ALIGN_DNA}.reduced"
```

### 8b — RAxML: Optimise GTRGAMMA on Fixed Tree Topology

RAxML evaluates the fixed IQ-TREE3 topology under GTRGAMMA to produce the info file needed by SEPP/pplacer.

```bash
rm -f $BASE/06_phylogeny/RAxML_info.sepp_eval \
      $BASE/06_phylogeny/RAxML_result.sepp_eval \
      $BASE/06_phylogeny/RAxML_bestTree.sepp_eval \
      $BASE/06_phylogeny/RAxML_parsimonyTree.sepp_eval \
      $BASE/06_phylogeny/RAxML_log.sepp_eval

raxmlHPC \
    -g $BASE/06_phylogeny/sulfurimonas_final.treefile \
    -s $ALIGN_RED -m GTRGAMMA -n sepp_eval -p 42 \
    -w $BASE/06_phylogeny/

grep "alpha" $BASE/06_phylogeny/RAxML_info.sepp_eval

grep -v "^Partition: 0 with name: No Name Provided$" \
    $BASE/06_phylogeny/RAxML_info.sepp_eval \
    > $BASE/06_phylogeny/RAxML_info.sepp_final
```

### 8c — Run SEPP

```bash
rm -f $BASE/07_shortreads/placement/sulfurimonas_sepp*
export LC_ALL=C; export LANG=C

run_sepp.py \
    -t $BASE/06_phylogeny/sulfurimonas_final.treefile \
    -a $ALIGN_DNA \
    -f $BASE/07_shortreads/derep/Sulfurimonas_ASVs_lineplot.fasta \
    -r $BASE/06_phylogeny/RAxML_info.sepp_final \
    -A 25 -P 25 -x 8 \
    -d $BASE/07_shortreads/placement/ \
    -o sulfurimonas_sepp
```

### 8d — Restore Bootstrap Labels + Convert to Newick

```bash
# Restore original bootstrap labels via SEPP auto-generated rename script
cat $BASE/07_shortreads/placement/sulfurimonas_sepp_placement.json \
    | python $BASE/07_shortreads/placement/sulfurimonas_sepp_rename-json.py \
    > $BASE/07_shortreads/placement/sulfurimonas_sepp_placement_final.json

# Convert placement JSON to Newick (tog = toggling tool in pplacer suite)
GUPPY=$(find ~/.sepp/ -name "guppy" 2>/dev/null | head -1)
echo "Using guppy: $GUPPY"
$GUPPY tog \
    $BASE/07_shortreads/placement/sulfurimonas_sepp_placement_final.json \
    -o $BASE/07_shortreads/placement/sulfurimonas_sepp_placement_final.tog.tre

echo "Done: $BASE/07_shortreads/placement/sulfurimonas_sepp_placement_final.tog.tre"
```

---

## Step 9 — Pairwise 16S Identity + Species Cluster Annotation

Full-length sequences only (no ASVs, no outgroups) to avoid length bias. Hierarchical clustering (average linkage) at 98.6% → putative species clusters. iTOL color strip generated.

**Output:** `03_merged/pairwise_16S_identity.tsv`, `03_merged/species_clusters_986.tsv`, `07_shortreads/itol_species_clusters.txt`

```python
python3 -c "
from Bio import SeqIO
import itertools, csv, statistics

ALN = '$BASE/04_alignment/sulfurimonas_trimmed_deduped.fasta'
OUT = '$BASE/03_merged/pairwise_16S_identity.tsv'

aln = [rec for rec in SeqIO.parse(ALN,'fasta')
       if 'OUTGROUP' not in rec.id and 'ASV' not in rec.id]
print(f'Full-length sequences for pairwise: {len(aln)}')

results = []
for r1, r2 in itertools.combinations(aln, 2):
    s1 = str(r1.seq).upper(); s2 = str(r2.seq).upper()
    matches = sum(a==b for a,b in zip(s1,s2) if a not in '-.' and b not in '-.')
    total   = sum(1   for a,b in zip(s1,s2) if a not in '-.' and b not in '-.')
    results.append((r1.id, r2.id, round(matches/total*100, 2) if total else 0))

with open(OUT,'w') as f:
    f.write('seq1\tseq2\tidentity_pct\n')
    for r in results: f.write(f'{r[0]}\t{r[1]}\t{r[2]}\n')

idents = [r[2] for r in results]
print(f'Comparisons: {len(results)} | Min: {min(idents):.2f}% | Max: {max(idents):.2f}% | Mean: {statistics.mean(idents):.2f}%')
print(f'Pairs below 98.6%: {sum(1 for r in results if r[2] < 98.6)}')
"

python3 -c "
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
import csv, colorsys
from Bio import SeqIO
from collections import Counter

PAIRWISE  = '$BASE/03_merged/pairwise_16S_identity.tsv'
ALN       = '$BASE/04_alignment/sulfurimonas_trimmed_deduped.fasta'
OUT_CLUST = '$BASE/03_merged/species_clusters_986.tsv'
OUT_ITOL  = '$BASE/07_shortreads/itol_species_clusters.txt'
THRESHOLD = 98.6

ids = [rec.id for rec in SeqIO.parse(ALN,'fasta')
       if 'OUTGROUP' not in rec.id and 'ASV' not in rec.id]
n = len(ids); id_idx = {s:i for i,s in enumerate(ids)}

dist_mat = np.zeros((n,n))
with open(PAIRWISE) as f:
    next(f)
    for line in f:
        s1, s2, ident = line.strip().split('\t')
        if s1 in id_idx and s2 in id_idx:
            d = 1 - float(ident)/100
            i,j = id_idx[s1], id_idx[s2]
            dist_mat[i,j] = dist_mat[j,i] = d

Z = linkage(squareform(dist_mat), method='average')
clusters = fcluster(Z, t=1-(THRESHOLD/100), criterion='distance')
cluster_map = {ids[i]: int(clusters[i]) for i in range(n)}
n_clusters = len(set(clusters))
print(f'Sequences: {n} | Clusters at {THRESHOLD}%: {n_clusters}')

with open(OUT_CLUST,'w') as f:
    f.write('seq_id\tcluster\n')
    for sid, clust in sorted(cluster_map.items(), key=lambda x: x[1]):
        f.write(f'{sid}\t{clust}\n')

def get_color(i, n):
    r,g,b = __import__('colorsys').hsv_to_rgb(i/n, 0.7, 0.85)
    return '#{:02x}{:02x}{:02x}'.format(int(r*255),int(g*255),int(b*255))

sorted_clusters = sorted(set(clusters.tolist()))
cluster_colors = {c: get_color(i,n_clusters) for i,c in enumerate(sorted_clusters)}

lines = ['DATASET_COLORSTRIP','SEPARATOR TAB',
         'DATASET_LABEL\t16S Species Clusters (98.6%)','COLOR\t#aaaaaa',
         'LEGEND_TITLE\t16S identity clusters (98.6%)',
         'STRIP_WIDTH\t30','MARGIN\t2','DATA']
for sid, clust in cluster_map.items():
    lines.append(f'{sid}\t{cluster_colors[clust]}\tCluster_{clust}')

with open(OUT_ITOL,'w') as f: f.write('\n'.join(lines)+'\n')
print(f'iTOL cluster file: {OUT_ITOL}')
"
```

---

## Step 10 — iTOL Annotation Files

**Output files:**

| File | Content |
|------|---------|
| `itol_colorstrip.txt` | Habitat category color strips |
| `itol_symbols.txt` | PS137 ASV symbols (chimney vs plume) |
| `itol_rename.txt` | Clean labels (accession / species name) |
| `itol_species_clusters.txt` | 16S identity clusters at 98.6% |

**Color scheme:**

| Category | Color |
|----------|-------|
| HV active chimney | `#D62728` (red) |
| HV fluids | `#FF7F0E` (orange) |
| HV pelagic | `#FFBB78` (light orange) |
| HV benthic | `#8C564B` (brown) |
| HV artificial | `#BCBD22` (olive) |
| HV inactive chimney | `#E377C2` (pink) |
| TypeStrain | `#000000` (black) |
| Outgroup | `#888888` (grey) |
| PS137 Chimney ASVs | `#FFD700` (gold) |
| PS137 Plume ASVs | `#7FFF00` (chartreuse) |

**PS137 ASVs:**
- Chimney only: `ASV35`, `ASV66`
- Plume only: `ASV1`, `ASV2`

```python
python3 -c "
import csv

META_CLEAN = '$BASE/03_merged/unified_metadata_clean.tsv'
SILVA_SEL  = '$BASE/01_silva/silva_selection_metadata.csv'
OUTDIR     = '$BASE/07_shortreads'

ref_colors = {
    'HV active chimney':   '#D62728',
    'HV fluids':           '#FF7F0E',
    'HV pelagic':          '#FFBB78',
    'HV benthic':          '#8C564B',
    'HV artificial':       '#BCBD22',
    'HV inactive chimney': '#E377C2',
    'TypeStrain':          '#000000',
    'Outgroup':            '#888888',
}
asv_colors = {
    'PS137 Chimney': '#FFD700',
    'PS137 Plume':   '#7FFF00',
}

chimney_only = {'ASV35','ASV66'}
plume_only   = {'ASV1','ASV2'}

VALID_SPECIES_PREFIXES = {
    'Sulfurimonas denitrificans':  'S. denitrificans',
    'Sulfurimonas autotrophica':   'S. autotrophica',
    'Sulfurimonas gotlandica':     'S. gotlandica',
    'Sulfurimonas marisnigri':     'S. marisnigri',
    'Sulfurimonas hongkongensis':  'S. hongkongensis',
    'Sulfurimonas paralvinellae':  'S. paralvinellae',
    'Sulfurimonas indica':         'S. indica',
}

def match_species(sp):
    for prefix, short in VALID_SPECIES_PREFIXES.items():
        if sp.startswith(prefix): return short
    return None

species_map = {}
with open(SILVA_SEL) as f:
    for row in csv.DictReader(f):
        acc = row['acc'].strip()
        sp  = row['species_silva'].strip()
        short = match_species(sp)
        if short: species_map[f'SILVA_{acc}'] = short

all_cats = list(ref_colors.keys()) + list(asv_colors.keys())
all_colors_map = {**ref_colors, **asv_colors}
legend_colors = '\t'.join(all_colors_map[c] for c in all_cats)
legend_shapes = '\t'.join('1' for _ in all_cats)
legend_labels = '\t'.join(all_cats)

strip = ['DATASET_COLORSTRIP','SEPARATOR TAB','DATASET_LABEL\tHabitat','COLOR\t#aaaaaa',
         'LEGEND_TITLE\tHabitat category',
         f'LEGEND_SHAPES\t{legend_shapes}',
         f'LEGEND_COLORS\t{legend_colors}',
         f'LEGEND_LABELS\t{legend_labels}',
         'STRIP_WIDTH\t40','MARGIN\t5','BORDER_WIDTH\t0','DATA']
sym = ['DATASET_SYMBOL','SEPARATOR TAB','DATASET_LABEL\tPS137 ASVs','COLOR\t#aaaaaa',
       'LEGEND_TITLE\tPS137 ASVs','LEGEND_SHAPES\t4\t4',
       f'LEGEND_COLORS\t{asv_colors[\"PS137 Chimney\"]}\t{asv_colors[\"PS137 Plume\"]}',
       'LEGEND_LABELS\tPS137 Chimney\tPS137 Plume','DATA']
ren = ['LABELS','SEPARATOR TAB','DATA']

with open(META_CLEAN) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        label = row['id'].strip(); cat = row['final_category'].strip()
        acc   = row['acc'].strip(); src = row['source'].strip()
        if src == 'OUTGROUP':
            strip.append(f'{label}\t{ref_colors[\"Outgroup\"]}\tOutgroup')
            ren.append(f'{label}\t{label.replace(\"OUTGROUP_\",\"\").replace(\"_\",\" \")}')
        else:
            col = ref_colors.get(cat,'#EEEEEE')
            strip.append(f'{label}\t{col}\t{cat}')
            name = species_map.get(label, acc)
            ren.append(f'{label}\t{name}')

for asv_id, cat in [*[(a,'PS137 Chimney') for a in chimney_only],
                     *[(a,'PS137 Plume') for a in plume_only]]:
    col = asv_colors[cat]
    strip.append(f'{asv_id}\t{col}\t{cat}')
    sym.append(f'{asv_id}\t4\t14\t{col}\t1\t1')
    ren.append(f'{asv_id}\t★ {asv_id}_PS137')

with open(f'{OUTDIR}/itol_colorstrip.txt','w') as f: f.write('\n'.join(strip)+'\n')
with open(f'{OUTDIR}/itol_symbols.txt','w') as f:    f.write('\n'.join(sym)+'\n')
with open(f'{OUTDIR}/itol_rename.txt','w') as f:     f.write('\n'.join(ren)+'\n')
print(f'iTOL annotation files written to {OUTDIR}/')
print(f'Species names applied to {len(species_map)} sequences')
"
```

---

## References

- Kim et al. (2014) — 98.6% identity threshold for 16S rRNA species-level clustering
- Tan et al. (2015) — automated trimming reduces phylogenetic accuracy for single-marker genes
- Mirarab et al. (2012) — SEPP: SATé-Enabled Phylogenetic Placement
- Quast et al. (2013) — SILVA ribosomal RNA database (138.2)

---

*Database: SILVA 138.2 SSURef NR99 | Expedition: PS137 Aurora Vent Field*
