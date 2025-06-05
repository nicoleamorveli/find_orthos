#  Find orthologues

This tool analyzes Pfam protein families against OMA (Orthologous MAtrix) ortholog groups to find missing orthologues.

---

## Overview

This tool identifies potential orthologs that OMA captures but Pfam families don't by:

- Fetching all proteins in a Pfam family using bulk UniProt queries  
- Extracting OMA fingerprints from those proteins  
- Finding proteins with the same OMA fingerprints that lack the Pfam domain  
- Generating comprehensive reports of the analysis  

---

##  Usage


```bash
python fast_pfam_oma_analyzer.py PF10181
```
### Advanced Usage
Custom minimum count threshold (mininum OMA appearences in the Uniprot)

```bash
python find_orthologues.py PF10181 -c 5
```
Custom report output name 
```bash
python find_orthologues.py PF10181 -o output.txt
```

```bash
python find_orthologues.py PF10181 --verbose
```
## 🧾 Command Line Options

| Option              | Description                             | Default     |
|---------------------|-----------------------------------------|-------------|
| `pfam_id`           | Pfam family ID (e.g., PF10181)          | Required    |
| `--min-count`, `-c` | Minimum occurrences for OMA groups      | 3           |
| `--output`, `-o`    | Output report file                      | report.txt  |
| `--verbose`, `-v`   | Enable detailed logging                 | False       |

## 📄 Example Output

```
SUMMARY:
Pfam PF10181: 1,247 proteins
OMA groups (>=3 occurrences): 23
Proteins unique to OMA: 156
```

---

## Location

`/homes/nflores/find_orthologues.py`

---

##  Requirements

- Python 3.6+  
- `requests` library


---

##  Considerations

- The default `count` is 3.
- Lower count thresholds (e.g., 1) may introduce many false positives.
- However, the report includes full UniProt descriptions to help identify false positives more easily.

**Example**: In `PF10181` with `count = 1`:

- `A0A6P4EDK6` (from *Drosophila*) is found.
- But many other false positives appear. These proteins (not sure about their relevance) also appear:
  - `A0A064BYA0`
  - `A0A174IDN3` (both entries are found by `PF04087`)


 
