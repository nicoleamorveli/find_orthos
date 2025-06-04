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

## Considerations
It only takes uniprot hits from pfam seq, but OMA uses the whole uniprot database so a lot of hits may seem not be found by the selected
Pfam but when cheking Interpro they are there.
The default of counts is 3 but when using counts =1 a lot of hits not present in Pfam (Interpro) are retrieved. 
A lot of false positives are also found, but in the report the uniprot description is available to detect them faster.
Eg: in PF10181 (c=1), A0A6P4EDK6 (Drosophila) is found but a lot of false positives are also found. 
A0A064BYA0 and A0A174IDN3 (hits in PF04087) not sure if these are relevant


