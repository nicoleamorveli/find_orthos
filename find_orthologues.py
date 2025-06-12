"""
Fast Pfam OMA Ortholog Analyzer

This script analyzes Pfam families by comparing them with OMA (Orthologous Matrix) data
using bulk UniProt API queries for much faster performance.

Usage:
    python fast_pfam_oma_analyzer.py PF10181
    python fast_pfam_oma_analyzer.py PF10181 --min-count 5 --output my_report.txt
"""

import sys
import time
import json
import gzip
from collections import Counter, defaultdict
import argparse
import logging
from typing import Dict, List, Set, Tuple, Optional
import re
import csv
from io import StringIO

try:
    import requests
except ImportError:
    print("Error: 'requests' module is required but not installed.")
    print("Please install it using: pip3 install requests")
    sys.exit(1)

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FastPfamOMAAnalyzer:
    def __init__(self):
    
        self.uniprot_base_url = "https://rest.uniprot.org"
        self.request_delay = 0.1  # Small delay to be respectful
    def get_oma_proteins(self, oma_fingerprint: Optional[str] = None, pfam_id: Optional[str] = None, in_pfam: bool = True) -> Tuple[List[Dict], Set[str]]:
        
        """
        Unified function to get proteins based on OMA fingerprint and/or Pfam domain,
        with UniProt status (Swiss-Prot or TrEMBL) included.

        Returns:
            Tuple of (protein_data_list, oma_fingerprints_set)
        """
        query_parts = []

        if oma_fingerprint:
            query_parts.append(f"xref:oma-{oma_fingerprint}")

        if pfam_id:
            if in_pfam:
                query_parts.append(pfam_id)
                operation = "AND"
            else:
                query_parts.append(f"NOT {pfam_id}")
                operation = "AND"

        query = f"({query_parts[0]})" if len(query_parts) == 1 else f"({query_parts[0]}) {operation} {query_parts[1]}"

        url = f"{self.uniprot_base_url}/uniprotkb/stream"
        params = {
            'compressed': 'true',
            'fields': 'accession,id,protein_name,xref_pfam,xref_oma,reviewed',
            'format': 'tsv',
            'query': query
            }

        try:
            response = requests.get(url, params=params, stream=True, timeout=60)
            response.raise_for_status()

            content = None
            if response.headers.get('content-encoding') == 'gzip':
                try:
                    content = gzip.decompress(response.content).decode('utf-8', errors='replace')
                except Exception:
                    params['compressed'] = 'false'
                    response = requests.get(url, params=params, stream=True, timeout=60)
                    response.raise_for_status()
                    content = response.text
            else:
                content = response.text

            content = content.replace('\x00', '').replace('\r', '')
            content = ''.join(c for c in content if c.isprintable() or c in '\t\n')

            if not content.strip():
                return [], set()

            proteins = []
            oma_fingerprints = set()

            try:
                csv_reader = csv.DictReader(StringIO(content), delimiter='\t')

                for row_num, row in enumerate(csv_reader, 1):
                    try:
                        accession = row.get('Entry', '').strip()
                        entry_name = row.get('Entry Name', '').strip()
                        protein_name = row.get('Protein names', '').strip()
                        pfam_refs = row.get('Pfam', '').strip()
                        oma_refs = row.get('OMA', '').strip()
                        reviewed = row.get('Reviewed', '').strip().lower()
                        uniprot_status = "Swiss-Prot" if reviewed == "reviewed" else "TrEMBL"

                        if not accession:
                            continue

                        protein_data = {
                            'accession': accession,
                            'entry_name': entry_name,
                            'protein_name': protein_name,
                            'pfam_refs': pfam_refs,
                            'oma_refs': oma_refs,
                            'uniprot_status': uniprot_status
                        }

                        proteins.append(protein_data)

                        if oma_refs:
                            oma_matches = re.findall(r'[A-Z]{7}', oma_refs)
                            oma_fingerprints.update(oma_matches)

                    except Exception as e:
                        logger.debug(f"Error processing row {row_num}: {e}")
                        continue

            except csv.Error as e:
                logger.error(f"CSV parsing error: {e}")
                lines = content.split('\n')
                if len(lines) > 1:
                    for line_num, line in enumerate(lines[1:], 2):
                        try:
                            if line.strip():
                                parts = line.split('\t')
                                if len(parts) >= 1 and parts[0].strip():
                                    accession = parts[0].strip()
                                    oma_refs = parts[4].strip() if len(parts) > 4 else ''
                                    reviewed = parts[5].strip().lower() if len(parts) > 5 else ''
                                    uniprot_status = "Swiss-Prot" if reviewed == "reviewed" else "TrEMBL"

                                    proteins.append({
                                        'accession': accession,
                                        'entry_name': parts[1].strip() if len(parts) > 1 else '',
                                        'protein_name': parts[2].strip() if len(parts) > 2 else '',
                                        'pfam_refs': parts[3].strip() if len(parts) > 3 else '',
                                        'oma_refs': oma_refs,
                                        'uniprot_status': uniprot_status
                                    })

                                    if oma_refs:
                                        oma_matches = re.findall(r'[A-Z]{7}', oma_refs)
                                        oma_fingerprints.update(oma_matches)
                        except Exception as e:
                            logger.debug(f"Error processing line {line_num}: {e}")
                            continue

            return proteins, oma_fingerprints

        except requests.RequestException as e:
            logger.error(f"Error fetching proteins: {e}")
            return [], set()
        except Exception as e:
            logger.error(f"Error parsing protein data: {e}")
            if 'content' in locals() and content:
                logger.debug(f"Content preview (first 200 chars): {repr(content[:200])}")
            return [], set()

   

    def get_total_oma_group_size(self, oma_fingerprint: str) -> int:
        """
        Get the total number of proteins in an OMA group.
        """
        logger.debug(f"Getting total size for OMA group {oma_fingerprint}")
        
        url = f"{self.uniprot_base_url}/uniprotkb/search"
        params = {
            'query': f'(xref:oma-{oma_fingerprint})',
            'format': 'json',
            'size': '1'  # We only need the count
        }
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            total_count = data.get('count', 0)
            
            logger.debug(f"OMA group {oma_fingerprint} has {total_count} total members")
            return total_count
            
        except requests.RequestException as e:
            logger.error(f"Error getting OMA group size for {oma_fingerprint}: {e}")
            return 0
        except Exception as e:
            logger.error(f"Error parsing OMA group size for {oma_fingerprint}: {e}")
            return 0

    def analyze_pfam_family(self, pfam_id: str, min_count: int = 3) -> Dict:
        """
        Analyze a Pfam family against OMA groups using bulk API queries.
        """
        logger.info(f"Starting fast analysis of Pfam family {pfam_id}")
        
        # Validate Pfam ID format
        if not re.match(r'^PF\d{5}$', pfam_id):
            logger.error(f"Invalid Pfam ID format: {pfam_id}. Expected format: PF#####")
            return {}
        
        # Step 1: Get all proteins with this Pfam domain and their OMA fingerprints
        pfam_proteins, oma_fingerprints = self.get_oma_proteins(pfam_id=pfam_id)
        
        if not pfam_proteins:
            logger.error(f"No proteins found for Pfam {pfam_id}")
            return {}
        
        if not oma_fingerprints:
            logger.warning(f"No OMA fingerprints found for proteins in Pfam {pfam_id}")
            return {
                'pfam_id': pfam_id,
                'pfam_protein_count': len(pfam_proteins),
                'pfam_proteins': pfam_proteins,
                'oma_fingerprints': {},
                'unique_to_oma': {},
                'unique_to_oma_count': 0
            }
        
        # Step 2: Count occurrences of each OMA fingerprint in the Pfam
        oma_counts = Counter()
        pfam_accessions = set()
        
        for protein in pfam_proteins:
            pfam_accessions.add(protein['accession'])
            if protein['oma_refs']:
                oma_matches = re.findall(r'[A-Z]{7}', protein['oma_refs'])
                oma_counts.update(oma_matches)
        
        # Step 3: Filter OMA groups by minimum count
        frequent_omas = {oma: count for oma, count in oma_counts.items() 
                        if count >= min_count}
        
        logger.info(f"Found {len(frequent_omas)} OMA groups with at least {min_count} occurrences")
        
        # Step 4: For each frequent OMA group, find proteins not in the Pfam
        unique_to_oma = {}
        total_unique_count = 0
        
        for oma_fingerprint, pfam_count in frequent_omas.items():
            logger.info(f"Processing OMA group {oma_fingerprint} ({pfam_count} in Pfam)")
            
            # Get proteins with this OMA fingerprint that don't have the Pfam domain
            oma_only_proteins, _ = self.get_oma_proteins(oma_fingerprint=oma_fingerprint, pfam_id=pfam_id, in_pfam=False)
            
            # Get total OMA group size
            total_oma_size = self.get_total_oma_group_size(oma_fingerprint)
            
            unique_to_oma[oma_fingerprint] = {
                'pfam_count': pfam_count,
                'total_oma_size': total_oma_size,
                'unique_proteins': oma_only_proteins,
                'unique_count': len(oma_only_proteins)
            }
            
            total_unique_count += len(oma_only_proteins)
            
            # Small delay to be respectful
            time.sleep(self.request_delay)
        
        results = {
            'pfam_id': pfam_id,
            'min_count': min_count,
            'pfam_protein_count': len(pfam_proteins),
            'pfam_proteins': pfam_proteins,
            'oma_fingerprints': dict(frequent_omas),
            'unique_to_oma': unique_to_oma,
            'unique_to_oma_count': total_unique_count,
            'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        logger.info(f"Analysis complete. Found {total_unique_count} proteins unique to OMA groups")
        return results

    def generate_report(self, results: Dict, output_file: str = "report.txt") -> str:
        """
        Generate a comprehensive report of the analysis and write directly to file.
        """
        if not results:
            return "No results to report."

        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                # Write header
                f.write("="*80 + "\n")
                f.write("FAST PFAM-OMA ORTHOLOG ANALYSIS REPORT\n")
                f.write("="*80 + "\n")
                f.write(f"Pfam Family: {results['pfam_id']}\n")
                f.write(f"Minimum Count Threshold: {results['min_count']}\n")
                f.write(f"Analysis Date: {results.get('analysis_timestamp', 'Unknown')}\n")
                f.write("\n")
                
                # Write summary statistics
                f.write("SUMMARY STATISTICS\n")
                f.write("-" * 50 + "\n")
                f.write(f"Total proteins in Pfam family: {results['pfam_protein_count']}\n")
                f.write(f"OMA groups with >={results['min_count']} occurrences: {len(results['oma_fingerprints'])}\n")
                f.write(f"Proteins unique to OMA groups: {results['unique_to_oma_count']}\n")
                f.write("\n")

                # Write frequent OMA groups section
                if results['oma_fingerprints']:
                    f.write(f"FREQUENT OMA GROUPS (>={results['min_count']} occurrences)\n")
                    f.write("-" * 50 + "\n")
                    
                    # Sort by count in descending order
                    sorted_omas = sorted(results['oma_fingerprints'].items(), 
                                       key=lambda x: x[1], reverse=True)
                    
                    for oma_fingerprint, pfam_count in sorted_omas:
                        oma_data = results['unique_to_oma'].get(oma_fingerprint, {})
                        total_size = oma_data.get('total_oma_size', 0)
                        unique_count = oma_data.get('unique_count', 0)
                        
                        f.write(f"{oma_fingerprint}: {pfam_count} in Pfam, "
                               f"{total_size} total in OMA, {unique_count} unique to OMA\n")
                    f.write("\n")

                # Write unique proteins section
                if results['unique_to_oma_count'] > 0:
                    f.write("PROTEINS UNIQUE TO OMA GROUPS (Not in Pfam)\n")
                    f.write("-" * 50 + "\n")
                    
                    # Sort OMA groups by fingerprint for consistent output
                    for oma_fingerprint in sorted(results['unique_to_oma'].keys()):
                        oma_data = results['unique_to_oma'][oma_fingerprint]
                        unique_proteins = oma_data['unique_proteins']
                        
                        if unique_proteins:
                            f.write(f"OMA GROUP: {oma_fingerprint}\n")
                            f.write("=" * 60 + "\n")
                            
                            # Sort proteins by accession
                            sorted_proteins = sorted(unique_proteins, key=lambda x: x['accession'])
                            
                            for protein in sorted_proteins:
                                accession = protein['accession']
                                protein_name = protein['protein_name'][:80] + "..." if len(protein['protein_name']) > 80 else protein['protein_name']
                                f.write(f"      {accession} | {protein_name}\n")
                            
                            f.write("\n")
                else:
                    f.write("No proteins found that are unique to OMA groups.\n")
                    f.write("\n")

                # Write footer
                f.write("="*80 + "\n")

            logger.info(f"Report saved to {output_file}")
            return f"Report successfully written to {output_file}"
            
        except IOError as e:
            logger.error(f"Failed to write report to {output_file}: {e}")
            return f"Error: Failed to write report to {output_file}: {e}"

def main():
    parser = argparse.ArgumentParser(
        description="Fast analysis of Pfam families vs OMA ortholog groups using bulk UniProt API queries"
    )
    parser.add_argument("pfam_id", help="Pfam family ID (e.g., PF10181)")
    parser.add_argument("--output", "-o", default="report.txt", 
                       help="Output file for the report (default: report.txt)")
    parser.add_argument("--min-count", "-c", type=int, default=3, 
                       help="Minimum count threshold for OMA groups (default: 3)")
    parser.add_argument("--verbose", "-v", action="store_true", 
                       help="Enable verbose logging")

    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate Pfam ID format
    if not re.match(r'^PF\d{5}$', args.pfam_id):
        logger.error(f"Invalid Pfam ID format: {args.pfam_id}. Expected format: PF##### (e.g., PF10181)")
        sys.exit(1)

    try:
        analyzer = FastPfamOMAAnalyzer()
        results = analyzer.analyze_pfam_family(args.pfam_id, min_count=args.min_count)
        
        if not results:
            logger.error("Analysis failed - no results generated")
            sys.exit(1)

        report_status = analyzer.generate_report(results, args.output)
        print(f"\nAnalysis complete! {report_status}")
        
        # Print summary to console
        print(f"\nSUMMARY:")
        print(f"Pfam {args.pfam_id}: {results['pfam_protein_count']} proteins")
        print(f"OMA groups (>={args.min_count} occurrences): {len(results['oma_fingerprints'])}")
        print(f"Proteins unique to OMA: {results['unique_to_oma_count']}")

    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Analysis failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

[nflores@codon-slurm-login-02 ~]$ nano prueba.py 
[nflores@codon-slurm-login-02 ~]$ ls
final_version_ortho    old_orthos  PF09356  PF23463  PF24894    PF25382    PF25459  prueba.py       report.txt
find_orthologues_1.py  ondemand    PF10181  PF23464  PF25203_2  PF25435    PF25803  report_PF20857
find_orthologues.py    PF03142     PF13632  PF23625  PF25203_C  PF25455_N  PF26560  report_PF22212
[nflores@codon-slurm-login-02 ~]$ rm prueba.py
[nflores@codon-slurm-login-02 ~]$ nano prueba.py
[nflores@codon-slurm-login-02 ~]$ python3 prueba.py
Traceback (most recent call last):
  File "prueba.py", line 364, in <module>
    main()
NameError: name 'main' is not defined
[nflores@codon-slurm-login-02 ~]$ rm prueba.py
[nflores@codon-slurm-login-02 ~]$ nano prueba.py
[nflores@codon-slurm-login-02 ~]$ python3 prueba.py
usage: prueba.py [-h] [--output OUTPUT] [--min-count MIN_COUNT] [--verbose]
                 pfam_id
prueba.py: error: the following arguments are required: pfam_id
[nflores@codon-slurm-login-02 ~]$ python3 prueba.py PF10181
2025-06-12 20:36:11,917 - ERROR - Analysis failed with error: name 'FastPfamOMAAnalyzer' is not defined
[nflores@codon-slurm-login-02 ~]$ cat find_orthologues.py
"""
Fast Pfam OMA Ortholog Analyzer

This script analyzes Pfam families by comparing them with OMA (Orthologous Matrix) data
using bulk UniProt API queries for much faster performance.

Usage:
    python fast_pfam_oma_analyzer.py PF10181
    python fast_pfam_oma_analyzer.py PF10181 --min-count 5 --output my_report.txt
"""

import sys
import time
import json
import gzip
from collections import Counter, defaultdict
import argparse
import logging
from typing import Dict, List, Set, Tuple, Optional
import re
import csv
from io import StringIO

try:
    import requests
except ImportError:
    print("Error: 'requests' module is required but not installed.")
    print("Please install it using: pip3 install requests")
    sys.exit(1)

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FastPfamOMAAnalyzer:
    def __init__(self):
    
        self.uniprot_base_url = "https://rest.uniprot.org"
        self.request_delay = 0.1  # Small delay to be respectful
    def get_oma_proteins(self, oma_fingerprint: Optional[str] = None, pfam_id: Optional[str] = None, in_pfam: bool = True) -> Tuple[List[Dict], Set[str]]:
        
        """
        Unified function to get proteins based on OMA fingerprint and/or Pfam domain,
        with UniProt status (Swiss-Prot or TrEMBL) included.

        Returns:
            Tuple of (protein_data_list, oma_fingerprints_set)
        """
        query_parts = []

        if oma_fingerprint:
            query_parts.append(f"xref:oma-{oma_fingerprint}")

        if pfam_id:
            if in_pfam:
                query_parts.append(pfam_id)
                operation = "AND"
            else:
                query_parts.append(f"NOT {pfam_id}")
                operation = "AND"

        query = f"({query_parts[0]})" if len(query_parts) == 1 else f"({query_parts[0]}) {operation} {query_parts[1]}"

        url = f"{self.uniprot_base_url}/uniprotkb/stream"
        params = {
            'compressed': 'true',
            'fields': 'accession,id,protein_name,xref_pfam,xref_oma,reviewed',
            'format': 'tsv',
            'query': query
            }

        try:
            response = requests.get(url, params=params, stream=True, timeout=60)
            response.raise_for_status()

            content = None
            if response.headers.get('content-encoding') == 'gzip':
                try:
                    content = gzip.decompress(response.content).decode('utf-8', errors='replace')
                except Exception:
                    params['compressed'] = 'false'
                    response = requests.get(url, params=params, stream=True, timeout=60)
                    response.raise_for_status()
                    content = response.text
            else:
                content = response.text

            content = content.replace('\x00', '').replace('\r', '')
            content = ''.join(c for c in content if c.isprintable() or c in '\t\n')

            if not content.strip():
                return [], set()

            proteins = []
            oma_fingerprints = set()

            try:
                csv_reader = csv.DictReader(StringIO(content), delimiter='\t')

                for row_num, row in enumerate(csv_reader, 1):
                    try:
                        accession = row.get('Entry', '').strip()
                        entry_name = row.get('Entry Name', '').strip()
                        protein_name = row.get('Protein names', '').strip()
                        pfam_refs = row.get('Pfam', '').strip()
                        oma_refs = row.get('OMA', '').strip()
                        reviewed = row.get('Reviewed', '').strip().lower()
                        uniprot_status = "Swiss-Prot" if reviewed == "reviewed" else "TrEMBL"

                        if not accession:
                            continue

                        protein_data = {
                            'accession': accession,
                            'entry_name': entry_name,
                            'protein_name': protein_name,
                            'pfam_refs': pfam_refs,
                            'oma_refs': oma_refs,
                            'uniprot_status': uniprot_status
                        }

                        proteins.append(protein_data)

                        if oma_refs:
                            oma_matches = re.findall(r'[A-Z]{7}', oma_refs)
                            oma_fingerprints.update(oma_matches)

                    except Exception as e:
                        logger.debug(f"Error processing row {row_num}: {e}")
                        continue

            except csv.Error as e:
                logger.error(f"CSV parsing error: {e}")
                lines = content.split('\n')
                if len(lines) > 1:
                    for line_num, line in enumerate(lines[1:], 2):
                        try:
                            if line.strip():
                                parts = line.split('\t')
                                if len(parts) >= 1 and parts[0].strip():
                                    accession = parts[0].strip()
                                    oma_refs = parts[4].strip() if len(parts) > 4 else ''
                                    reviewed = parts[5].strip().lower() if len(parts) > 5 else ''
                                    uniprot_status = "Swiss-Prot" if reviewed == "reviewed" else "TrEMBL"

                                    proteins.append({
                                        'accession': accession,
                                        'entry_name': parts[1].strip() if len(parts) > 1 else '',
                                        'protein_name': parts[2].strip() if len(parts) > 2 else '',
                                        'pfam_refs': parts[3].strip() if len(parts) > 3 else '',
                                        'oma_refs': oma_refs,
                                        'uniprot_status': uniprot_status
                                    })

                                    if oma_refs:
                                        oma_matches = re.findall(r'[A-Z]{7}', oma_refs)
                                        oma_fingerprints.update(oma_matches)
                        except Exception as e:
                            logger.debug(f"Error processing line {line_num}: {e}")
                            continue

            return proteins, oma_fingerprints

        except requests.RequestException as e:
            logger.error(f"Error fetching proteins: {e}")
            return [], set()
        except Exception as e:
            logger.error(f"Error parsing protein data: {e}")
            if 'content' in locals() and content:
                logger.debug(f"Content preview (first 200 chars): {repr(content[:200])}")
            return [], set()

   

    def get_total_oma_group_size(self, oma_fingerprint: str) -> int:
        """
        Get the total number of proteins in an OMA group.
        """
        logger.debug(f"Getting total size for OMA group {oma_fingerprint}")
        
        url = f"{self.uniprot_base_url}/uniprotkb/search"
        params = {
            'query': f'(xref:oma-{oma_fingerprint})',
            'format': 'json',
            'size': '1'  # We only need the count
        }
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            total_count = data.get('count', 0)
            
            logger.debug(f"OMA group {oma_fingerprint} has {total_count} total members")
            return total_count
            
        except requests.RequestException as e:
            logger.error(f"Error getting OMA group size for {oma_fingerprint}: {e}")
            return 0
        except Exception as e:
            logger.error(f"Error parsing OMA group size for {oma_fingerprint}: {e}")
            return 0

    def analyze_pfam_family(self, pfam_id: str, min_count: int = 3) -> Dict:
        """
        Analyze a Pfam family against OMA groups using bulk API queries.
        """
        logger.info(f"Starting fast analysis of Pfam family {pfam_id}")
        
        # Validate Pfam ID format
        if not re.match(r'^PF\d{5}$', pfam_id):
            logger.error(f"Invalid Pfam ID format: {pfam_id}. Expected format: PF#####")
            return {}
        
        # Step 1: Get all proteins with this Pfam domain and their OMA fingerprints
        pfam_proteins, oma_fingerprints = self.get_oma_proteins(pfam_id=pfam_id)
        
        if not pfam_proteins:
            logger.error(f"No proteins found for Pfam {pfam_id}")
            return {}
        
        if not oma_fingerprints:
            logger.warning(f"No OMA fingerprints found for proteins in Pfam {pfam_id}")
            return {
                'pfam_id': pfam_id,
                'pfam_protein_count': len(pfam_proteins),
                'pfam_proteins': pfam_proteins,
                'oma_fingerprints': {},
                'unique_to_oma': {},
                'unique_to_oma_count': 0
            }
        
        # Step 2: Count occurrences of each OMA fingerprint in the Pfam
        oma_counts = Counter()
        pfam_accessions = set()
        
        for protein in pfam_proteins:
            pfam_accessions.add(protein['accession'])
            if protein['oma_refs']:
                oma_matches = re.findall(r'[A-Z]{7}', protein['oma_refs'])
                oma_counts.update(oma_matches)
        
        # Step 3: Filter OMA groups by minimum count
        frequent_omas = {oma: count for oma, count in oma_counts.items() 
                        if count >= min_count}
        
        logger.info(f"Found {len(frequent_omas)} OMA groups with at least {min_count} occurrences")
        
        # Step 4: For each frequent OMA group, find proteins not in the Pfam
        unique_to_oma = {}
        total_unique_count = 0
        
        for oma_fingerprint, pfam_count in frequent_omas.items():
            logger.info(f"Processing OMA group {oma_fingerprint} ({pfam_count} in Pfam)")
            
            # Get proteins with this OMA fingerprint that don't have the Pfam domain
            oma_only_proteins, _ = self.get_oma_proteins(oma_fingerprint=oma_fingerprint, pfam_id=pfam_id, in_pfam=False)
            
            # Get total OMA group size
            total_oma_size = self.get_total_oma_group_size(oma_fingerprint)
            
            unique_to_oma[oma_fingerprint] = {
                'pfam_count': pfam_count,
                'total_oma_size': total_oma_size,
                'unique_proteins': oma_only_proteins,
                'unique_count': len(oma_only_proteins)
            }
            
            total_unique_count += len(oma_only_proteins)
            
            # Small delay to be respectful
            time.sleep(self.request_delay)
        
        results = {
            'pfam_id': pfam_id,
            'min_count': min_count,
            'pfam_protein_count': len(pfam_proteins),
            'pfam_proteins': pfam_proteins,
            'oma_fingerprints': dict(frequent_omas),
            'unique_to_oma': unique_to_oma,
            'unique_to_oma_count': total_unique_count,
            'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        logger.info(f"Analysis complete. Found {total_unique_count} proteins unique to OMA groups")
        return results

    def generate_report(self, results: Dict, output_file: str = "report.txt") -> str:
        """
        Generate a comprehensive report of the analysis and write directly to file.
        """
        if not results:
            return "No results to report."

        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                # Write header
                f.write("="*80 + "\n")
                f.write("FAST PFAM-OMA ORTHOLOG ANALYSIS REPORT\n")
                f.write("="*80 + "\n")
                f.write(f"Pfam Family: {results['pfam_id']}\n")
                f.write(f"Minimum Count Threshold: {results['min_count']}\n")
                f.write(f"Analysis Date: {results.get('analysis_timestamp', 'Unknown')}\n")
                f.write("\n")
                
                # Write summary statistics
                f.write("SUMMARY STATISTICS\n")
                f.write("-" * 50 + "\n")
                f.write(f"Total proteins in Pfam family: {results['pfam_protein_count']}\n")
                f.write(f"OMA groups with >={results['min_count']} occurrences: {len(results['oma_fingerprints'])}\n")
                f.write(f"Proteins unique to OMA groups: {results['unique_to_oma_count']}\n")
                f.write("\n")

                # Write frequent OMA groups section
                if results['oma_fingerprints']:
                    f.write(f"FREQUENT OMA GROUPS (>={results['min_count']} occurrences)\n")
                    f.write("-" * 50 + "\n")
                    
                    # Sort by count in descending order
                    sorted_omas = sorted(results['oma_fingerprints'].items(), 
                                       key=lambda x: x[1], reverse=True)
                    
                    for oma_fingerprint, pfam_count in sorted_omas:
                        oma_data = results['unique_to_oma'].get(oma_fingerprint, {})
                        total_size = oma_data.get('total_oma_size', 0)
                        unique_count = oma_data.get('unique_count', 0)
                        
                        f.write(f"{oma_fingerprint}: {pfam_count} in Pfam, "
                               f"{total_size} total in OMA, {unique_count} unique to OMA\n")
                    f.write("\n")

                # Write unique proteins section
                if results['unique_to_oma_count'] > 0:
                    f.write("PROTEINS UNIQUE TO OMA GROUPS (Not in Pfam)\n")
                    f.write("-" * 50 + "\n")
                    
                    # Sort OMA groups by fingerprint for consistent output
                    for oma_fingerprint in sorted(results['unique_to_oma'].keys()):
                        oma_data = results['unique_to_oma'][oma_fingerprint]
                        unique_proteins = oma_data['unique_proteins']
                        
                        if unique_proteins:
                            f.write(f"OMA GROUP: {oma_fingerprint}\n")
                            f.write("=" * 60 + "\n")
                            
                            # Sort proteins by accession
                            sorted_proteins = sorted(unique_proteins, key=lambda x: x['accession'])
                            
                            for protein in sorted_proteins:
                                accession = protein['accession']
                                protein_name = protein['protein_name'][:80] + "..." if len(protein['protein_name']) > 80 else protein['protein_name']
                                f.write(f"      {accession} | {protein_name}\n")
                            
                            f.write("\n")
                else:
                    f.write("No proteins found that are unique to OMA groups.\n")
                    f.write("\n")

                # Write footer
                f.write("="*80 + "\n")

            logger.info(f"Report saved to {output_file}")
            return f"Report successfully written to {output_file}"
            
        except IOError as e:
            logger.error(f"Failed to write report to {output_file}: {e}")
            return f"Error: Failed to write report to {output_file}: {e}"

def main():
    parser = argparse.ArgumentParser(
        description="Fast analysis of Pfam families vs OMA ortholog groups using bulk UniProt API queries"
    )
    parser.add_argument("pfam_id", help="Pfam family ID (e.g., PF10181)")
    parser.add_argument("--output", "-o", default="report.txt", 
                       help="Output file for the report (default: report.txt)")
    parser.add_argument("--min-count", "-c", type=int, default=3, 
                       help="Minimum count threshold for OMA groups (default: 3)")
    parser.add_argument("--verbose", "-v", action="store_true", 
                       help="Enable verbose logging")

    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate Pfam ID format
    if not re.match(r'^PF\d{5}$', args.pfam_id):
        logger.error(f"Invalid Pfam ID format: {args.pfam_id}. Expected format: PF##### (e.g., PF10181)")
        sys.exit(1)

    try:
        analyzer = FastPfamOMAAnalyzer()
        results = analyzer.analyze_pfam_family(args.pfam_id, min_count=args.min_count)
        
        if not results:
            logger.error("Analysis failed - no results generated")
            sys.exit(1)

        report_status = analyzer.generate_report(results, args.output)
        print(f"\nAnalysis complete! {report_status}")
        
        # Print summary to console
        print(f"\nSUMMARY:")
        print(f"Pfam {args.pfam_id}: {results['pfam_protein_count']} proteins")
        print(f"OMA groups (>={args.min_count} occurrences): {len(results['oma_fingerprints'])}")
        print(f"Proteins unique to OMA: {results['unique_to_oma_count']}")

    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Analysis failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
