#!/usr/bin/env python3
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

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FastPfamOMAAnalyzer:
    def __init__(self):
        self.uniprot_base_url = "https://rest.uniprot.org"
        self.request_delay = 0.1  # Small delay to be respectful

    def get_pfam_proteins_bulk(self, pfam_id: str) -> Tuple[List[Dict], Set[str]]:
        """
        Get all proteins containing the specified Pfam domain using UniProt API.
        Returns tuple of (protein_data_list, oma_fingerprints_set)
        """
        logger.info(f"Fetching all proteins with Pfam domain {pfam_id}")
        
        # Construct the UniProt query URL for bulk download
        url = f"{self.uniprot_base_url}/uniprotkb/stream"
        params = {
            'compressed': 'true',
            'fields': 'accession,id,protein_name,xref_pfam,xref_oma',
            'format': 'tsv',
            'query': f'({pfam_id})'
        }
        
        try:
            response = requests.get(url, params=params, stream=True, timeout=60)
            response.raise_for_status()
            
            # Handle response content with better error handling
            content = None
            if response.headers.get('content-encoding') == 'gzip':
                try:
                    content = gzip.decompress(response.content).decode('utf-8', errors='replace')
                except Exception as e:
                    logger.warning(f"Error decompressing gzip for {pfam_id}: {e}")
                    # Try without compression
                    params['compressed'] = 'false'
                    response = requests.get(url, params=params, stream=True, timeout=60)
                    response.raise_for_status()
                    content = response.text
            else:
                content = response.text
            
            # Clean content: remove null bytes and other problematic characters
            if content:
                content = content.replace('\x00', '')  # Remove null bytes
                content = content.replace('\r', '')    # Remove carriage returns that might cause issues
                # Remove any other control characters except tabs and newlines
                content = ''.join(char for char in content if char.isprintable() or char in '\t\n')
            
            if not content or not content.strip():
                logger.error(f"Empty response for Pfam {pfam_id}")
                return [], set()
            
            # Parse TSV data
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
                        
                        if not accession:
                            continue
                        
                        protein_data = {
                            'accession': accession,
                            'entry_name': entry_name,
                            'protein_name': protein_name,
                            'pfam_refs': pfam_refs,
                            'oma_refs': oma_refs
                        }
                        
                        proteins.append(protein_data)
                        
                        # Extract OMA fingerprints
                        if oma_refs:
                            # OMA references are typically in format like "CYLTGQH"
                            oma_matches = re.findall(r'[A-Z]{7}', oma_refs)
                            oma_fingerprints.update(oma_matches)
                            
                    except Exception as e:
                        logger.debug(f"Error processing row {row_num} for {pfam_id}: {e}")
                        continue
                        
            except csv.Error as e:
                logger.error(f"CSV parsing error for {pfam_id}: {e}")
                # Try to salvage what we can by processing line by line
                lines = content.split('\n')
                if len(lines) > 1:  # Skip header
                    for line_num, line in enumerate(lines[1:], 2):
                        try:
                            if line.strip():
                                parts = line.split('\t')
                                if len(parts) >= 1 and parts[0].strip():
                                    accession = parts[0].strip()
                                    oma_refs = parts[4].strip() if len(parts) > 4 else ''
                                    
                                    proteins.append({
                                        'accession': accession,
                                        'entry_name': parts[1].strip() if len(parts) > 1 else '',
                                        'protein_name': parts[2].strip() if len(parts) > 2 else '',
                                        'pfam_refs': parts[3].strip() if len(parts) > 3 else '',
                                        'oma_refs': oma_refs
                                    })
                                    
                                    # Extract OMA fingerprints
                                    if oma_refs:
                                        oma_matches = re.findall(r'[A-Z]{7}', oma_refs)
                                        oma_fingerprints.update(oma_matches)
                        except Exception as e:
                            logger.debug(f"Error processing line {line_num}: {e}")
                            continue
            
            logger.info(f"Found {len(proteins)} proteins with {pfam_id}")
            logger.info(f"Found {len(oma_fingerprints)} unique OMA fingerprints")
            
            return proteins, oma_fingerprints
            
        except requests.RequestException as e:
            logger.error(f"Error fetching Pfam proteins: {e}")
            return [], set()
        except Exception as e:
            logger.error(f"Error parsing Pfam protein data: {e}")
            # Log the first few characters of content for debugging
            if 'content' in locals() and content:
                logger.debug(f"Content preview (first 200 chars): {repr(content[:200])}")
            return [], set()

    def get_oma_proteins_not_in_pfam(self, oma_fingerprint: str, pfam_id: str) -> List[Dict]:
        """
        Get all proteins with the OMA fingerprint that DON'T have the specified Pfam domain.
        """
        logger.info(f"Fetching proteins with OMA fingerprint {oma_fingerprint} NOT in {pfam_id}")
        
        url = f"{self.uniprot_base_url}/uniprotkb/stream"
        params = {
            'compressed': 'true',
            'fields': 'accession,id,protein_name,xref_pfam,xref_oma',
            'format': 'tsv',
            'query': f'(xref:oma-{oma_fingerprint}) NOT {pfam_id}'
        }
        
        try:
            response = requests.get(url, params=params, stream=True, timeout=60)
            response.raise_for_status()
            
            # Handle response content with better error handling
            content = None
            if response.headers.get('content-encoding') == 'gzip':
                try:
                    content = gzip.decompress(response.content).decode('utf-8', errors='replace')
                except Exception as e:
                    logger.warning(f"Error decompressing gzip for {oma_fingerprint}: {e}")
                    # Try without compression
                    params['compressed'] = 'false'
                    response = requests.get(url, params=params, stream=True, timeout=60)
                    response.raise_for_status()
                    content = response.text
            else:
                content = response.text
            
            # Clean content: remove null bytes and other problematic characters
            if content:
                content = content.replace('\x00', '')  # Remove null bytes
                content = content.replace('\r', '')    # Remove carriage returns that might cause issues
                # Remove any other control characters except tabs and newlines
                content = ''.join(char for char in content if char.isprintable() or char in '\t\n')
            
            if not content or not content.strip():
                logger.warning(f"Empty response for OMA fingerprint {oma_fingerprint}")
                return []
            
            proteins = []
            
            try:
                csv_reader = csv.DictReader(StringIO(content), delimiter='\t')
                
                for row_num, row in enumerate(csv_reader, 1):
                    try:
                        accession = row.get('Entry', '').strip()
                        entry_name = row.get('Entry Name', '').strip()
                        protein_name = row.get('Protein names', '').strip()
                        pfam_refs = row.get('Pfam', '').strip()  
                        oma_refs = row.get('OMA', '').strip()
                        
                        if not accession:
                            continue
                        
                        protein_data = {
                            'accession': accession,
                            'entry_name': entry_name,
                            'protein_name': protein_name,
                            'pfam_refs': pfam_refs,
                            'oma_refs': oma_refs
                        }
                        
                        proteins.append(protein_data)
                        
                    except Exception as e:
                        logger.debug(f"Error processing row {row_num} for {oma_fingerprint}: {e}")
                        continue
                        
            except csv.Error as e:
                logger.error(f"CSV parsing error for {oma_fingerprint}: {e}")
                # Try to salvage what we can by processing line by line
                lines = content.split('\n')
                if len(lines) > 1:  # Skip header
                    for line_num, line in enumerate(lines[1:], 2):
                        try:
                            if line.strip():
                                parts = line.split('\t')
                                if len(parts) >= 1 and parts[0].strip():
                                    proteins.append({
                                        'accession': parts[0].strip(),
                                        'entry_name': parts[1].strip() if len(parts) > 1 else '',
                                        'protein_name': parts[2].strip() if len(parts) > 2 else '',
                                        'pfam_refs': parts[3].strip() if len(parts) > 3 else '',
                                        'oma_refs': parts[4].strip() if len(parts) > 4 else ''
                                    })
                        except Exception as e:
                            logger.debug(f"Error processing line {line_num}: {e}")
                            continue
            
            logger.info(f"Found {len(proteins)} proteins with OMA {oma_fingerprint} not in {pfam_id}")
            return proteins
            
        except requests.RequestException as e:
            logger.error(f"Error fetching OMA proteins for {oma_fingerprint}: {e}")
            return []
        except Exception as e:
            logger.error(f"Error parsing OMA protein data for {oma_fingerprint}: {e}")
            # Log the first few characters of content for debugging
            if 'content' in locals() and content:
                logger.debug(f"Content preview (first 200 chars): {repr(content[:200])}")
            return []

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
        pfam_proteins, oma_fingerprints = self.get_pfam_proteins_bulk(pfam_id)
        
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
            oma_only_proteins = self.get_oma_proteins_not_in_pfam(oma_fingerprint, pfam_id)
            
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
        Generate a comprehensive report of the analysis.
        """
        if not results:
            return "No results to report."

        lines = []
        lines.append("="*80)
        lines.append("FAST PFAM-OMA ORTHOLOG ANALYSIS REPORT")
        lines.append("="*80)
        lines.append(f"Pfam Family: {results['pfam_id']}")
        lines.append(f"Minimum Count Threshold: {results['min_count']}")
        lines.append(f"Analysis Date: {results.get('analysis_timestamp', 'Unknown')}")
        lines.append("")
        
        lines.append("SUMMARY STATISTICS")
        lines.append("-" * 50)
        lines.append(f"Total proteins in Pfam family: {results['pfam_protein_count']}")
        lines.append(f"OMA groups with >={results['min_count']} occurrences: {len(results['oma_fingerprints'])}")
        lines.append(f"Proteins unique to OMA groups: {results['unique_to_oma_count']}")
        lines.append("")

        if results['oma_fingerprints']:
            lines.append(f"FREQUENT OMA GROUPS (>={results['min_count']} occurrences)")
            lines.append("-" * 50)
            
            # Sort by count in descending order
            sorted_omas = sorted(results['oma_fingerprints'].items(), 
                               key=lambda x: x[1], reverse=True)
            
            for oma_fingerprint, pfam_count in sorted_omas:
                oma_data = results['unique_to_oma'].get(oma_fingerprint, {})
                total_size = oma_data.get('total_oma_size', 0)
                unique_count = oma_data.get('unique_count', 0)
                
                lines.append(f"{oma_fingerprint}: {pfam_count} in Pfam, "
                           f"{total_size} total in OMA, {unique_count} unique to OMA")
            lines.append("")

        if results['unique_to_oma_count'] > 0:
            lines.append("PROTEINS UNIQUE TO OMA GROUPS (Not in Pfam)")
            lines.append("-" * 50)
            
            # Sort OMA groups by fingerprint for consistent output
            for oma_fingerprint in sorted(results['unique_to_oma'].keys()):
                oma_data = results['unique_to_oma'][oma_fingerprint]
                unique_proteins = oma_data['unique_proteins']
                
                if unique_proteins:
                    lines.append(f"OMA GROUP: {oma_fingerprint}")
                    lines.append("=" * 60)
                    
                    # Sort proteins by accession
                    sorted_proteins = sorted(unique_proteins, key=lambda x: x['accession'])
                    
                    for protein in sorted_proteins:
                        accession = protein['accession']
                        protein_name = protein['protein_name'][:80] + "..." if len(protein['protein_name']) > 80 else protein['protein_name']
                        lines.append(f"      {accession} | {protein_name}")
                    
                    lines.append("")
        else:
            lines.append("No proteins found that are unique to OMA groups.")
            lines.append("")

        lines.append("="*80)

        report_text = "\n".join(lines)
        
        # Save report to file
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(report_text)
            logger.info(f"Report saved to {output_file}")
        except IOError as e:
            logger.error(f"Failed to write report to {output_file}: {e}")
        
        return report_text

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

        report = analyzer.generate_report(results, args.output)
        print(f"\nAnalysis complete! Report saved to {args.output}")
        
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
