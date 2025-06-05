# find_orthos

This script analyzes Pfam families by comparing them with OMA (Orthologous MAtrix) data.  
It identifies potential orthologs that are captured by OMA but missed by Pfam.

**Input**: Pfam folder (containing the `scores` file)  
**Output**: UniProt IDs captured by OMA but not by Pfam

---

## Steps

1. Reads the `scores` file from the Pfam folder and retrieves the UniProt IDs (only hits in the `ALIGN` file).
2. Retrieves the OMA fingerprint for each UniProt ID.
3. Creates a list and filters it based on the number of appearances (default count = 3).
4. The count is user-defined and corresponds to how many times an OMA ID appears in the Pfam UniProt IDs.
5. Orders the list by number of appearances (i.e., relevance).
6. Creates another list by gathering all UniProt IDs associated with each OMA identifier.
7. Filters this list to keep only UniProt IDs that are uniquely associated with one OMA identifier.
8. These are unique hits captured by OMA but not Pfam.
9. Builds a report, including the counts of OMA appearances.

---

## Usage

usage: ortho_counts.py [-h] [--output OUTPUT] [--min-count MIN_COUNT] [--verbose]

Location 
- /homes/nflores/ortho_counts.py
- /homes/nflores/ortho_full_taxonomy.py 

## Considerations
- Only UniProt hits from the Pfam `seq` file are considered, while OMA uses the entire UniProt database.
- This means many hits may appear to be missing from the selected Pfam but can actually be found when checking InterPro.
- The default `count` is 3, but setting `count = 1` retrieves many additional hits not present in Pfam (as seen in InterPro).
- Lower count thresholds (e.g., 1) also introduce many false positives.
- However, the report includes full UniProt descriptions to help identify false positives more easily.
- **Example**: In `PF10181` with `count = 1`:
  - `A0A6P4EDK6` (from *Drosophila*) is found.
  - But many other false positives appear. These proteins (not sure about their relevance) also appear:
    - `A0A064BYA0`
    - `A0A174IDN3` (both entries are found by `PF04087`)

## New version 
this new version ortho_full_taxonomy.py gives the detailed taxonomy divided by genus (I think) but it's very long
 
