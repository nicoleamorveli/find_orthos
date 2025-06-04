#!/usr/bin/env python3
"""
Pfam OMA Ortholog Analyzer

This script analyzes Pfam families by comparing them with OMA (Orthologous MAtrix) data.
It identifies potential orthologs that OMA captures but Pfam doesn't.

Usage:
    python check_orthologs.py /path/to/PF00001
"""

import sys
import time
import json
from collections import Counter
import argparse
import logging
from typing import Dict, List, Set, Tuple, Optional
import re
import os

try:
    import requests
except ImportError:
    print("Error: 'requests' module is required but not installed.")
    print("Please install it using: pip3 install requests")
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PfamOMAAnalyzer:
    def __init__(self):
        self.oma_base_url = "https://omabrowser.org/api/"
        self.uniprot_base_url = "https://rest.uniprot.org"
        self.request_delay = 0.2  # Increased delay to be more respectful to APIs

    def read_pfam_scores_file(self, pfam_folder: str) -> List[str]:
        scores_file = os.path.join(pfam_folder, "scores")

        if not os.path.exists(scores_file):
            logger.error(f"Scores file not found: {scores_file}")
            return []

        logger.info(f"Reading UniProt IDs from {scores_file}")
        uniprot_ids = []

        try:
            with open(scores_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split()
                    if len(parts) < 2:
                        continue

                    potential_uniprot = parts[1]
                    uniprot_id = potential_uniprot.split('.')[0] if '.' in potential_uniprot else potential_uniprot

                    if re.match(r'^[A-Z0-9]{6,10}$', uniprot_id):
                        uniprot_ids.append(uniprot_id)
                    else:
                        logger.debug(f"Skipping invalid UniProt ID format on line {line_num}: {uniprot_id}")
        except Exception as e:
            logger.error(f"Error reading scores file: {e}")
            return []

        logger.info(f"Extracted {len(uniprot_ids)} UniProt IDs from scores file")
        return uniprot_ids

    def get_pfam_id_from_folder(self, pfam_folder: str) -> str:
        folder_name = os.path.basename(pfam_folder.rstrip('/'))
        if re.match(r'^PF\d{5}$', folder_name):
            return folder_name

        desc_file = os.path.join(pfam_folder, "DESC")
        if os.path.exists(desc_file):
            try:
                with open(desc_file, 'r') as f:
                    for line in f:
                        if line.startswith('AC '):
                            pfam_id = line.split()[1].rstrip(';')
                            if re.match(r'^PF\d{5}$', pfam_id):
                                return pfam_id
            except Exception as e:
                logger.debug(f"Error reading DESC file: {e}")

        return folder_name

    def get_uniprot_description(self, uniprot_id: str) -> str:
        """
        Get the protein description from UniProt.
        """
        try:
            url = f"{self.uniprot_base_url}/uniprotkb/{uniprot_id}"
            params = {'format': 'json'}
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                # Extract protein name/description
                if 'proteinDescription' in data:
                    desc = data['proteinDescription']
                    if 'recommendedName' in desc and 'fullName' in desc['recommendedName']:
                        return desc['recommendedName']['fullName']['value']
                    elif 'submissionNames' in desc and desc['submissionNames']:
                        return desc['submissionNames'][0]['fullName']['value']
                
                # Fallback to entry name
                if 'uniProtkbId' in data:
                    return data['uniProtkbId']
                    
            time.sleep(0.1)  # Small delay to respect API limits
            return "Unknown"
            
        except requests.RequestException as e:
            logger.debug(f"Error getting description for {uniprot_id}: {e}")
            return "Unknown"
        except Exception as e:
            logger.debug(f"Error parsing description for {uniprot_id}: {e}")
            return "Unknown"

    def get_oma_fingerprint(self, uniprot_id: str) -> Optional[str]:
        try:
            # First try the OMA API directly
            url = f"{self.oma_base_url}protein/{uniprot_id}/"
            response = requests.get(url, timeout=10)

            if response.status_code == 200:
                data = response.json()
                if 'oma_group' in data and data['oma_group']:
                    return data['oma_group']

            time.sleep(self.request_delay)
            
            # Fallback: scrape UniProt page for OMA links
            uniprot_url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}"
            response = requests.get(uniprot_url, timeout=10)

            if response.status_code == 200:
                oma_pattern = r'omabrowser\.org/oma/omagroup/([^/\s"]+)'
                match = re.search(oma_pattern, response.text)
                if match:
                    return match.group(1)

        except requests.RequestException as e:
            logger.debug(f"Error getting OMA fingerprint for {uniprot_id}: {e}")

        return None

    def batch_get_oma_fingerprints(self, uniprot_ids: List[str]) -> Dict[str, str]:
        logger.info(f"Fetching OMA fingerprints for {len(uniprot_ids)} UniProt IDs")
        oma_mapping = {}

        for i, uniprot_id in enumerate(uniprot_ids, 1):
            if i % 50 == 0:
                logger.info(f"Processed {i}/{len(uniprot_ids)} UniProt IDs")

            oma_id = self.get_oma_fingerprint(uniprot_id)
            if oma_id:
                oma_mapping[uniprot_id] = oma_id

            time.sleep(self.request_delay)

        logger.info(f"Found OMA fingerprints for {len(oma_mapping)} out of {len(uniprot_ids)} UniProt IDs")
        return oma_mapping

    def filter_frequent_omas(self, oma_mapping: Dict[str, str], min_count: int = 3) -> Dict[str, int]:
        """
        Filter OMA groups based on minimum count threshold instead of percentage.
        """
        oma_counts = Counter(oma_mapping.values())

        frequent_omas = {}
        for oma_id, count in oma_counts.items():
            if count >= min_count:
                frequent_omas[oma_id] = count

        frequent_omas = dict(sorted(frequent_omas.items(), key=lambda x: x[1], reverse=True))
        logger.info(f"Found {len(frequent_omas)} OMA groups with at least {min_count} occurrences")
        return frequent_omas

    def get_oma_fingerprint_from_group(self, oma_id: str) -> Optional[str]:
        """
        Get the OMA fingerprint for a given OMA group ID.
        """
        try:
            url = f"{self.oma_base_url}group/{oma_id}/"
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            # The fingerprint is typically the OMA group ID itself
            if 'fingerprint' in data:
                return data['fingerprint']
            else:
                # If no explicit fingerprint field, use the group ID
                return oma_id
                
        except requests.RequestException as e:
            logger.error(f"Error fetching OMA group {oma_id}: {e}")
            return None

    def get_uniprot_ids_with_oma_fingerprint(self, fingerprint: str) -> List[str]:
        """
        Search UniProt for all proteins with a specific OMA fingerprint.
        Uses UniProt's search API with the query: (xref:oma-FINGERPRINT)
        """
        try:
            # Correct UniProt search query for OMA cross-references
            query = f'(xref:oma-{fingerprint})'
            url = f"{self.uniprot_base_url}/uniprotkb/search"
            
            params = {
                'query': query,
                'format': 'json',
                'size': 500  # Adjust size as needed
            }
            
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            uniprot_ids = []
            
            if 'results' in data:
                for result in data['results']:
                    if 'primaryAccession' in result:
                        uniprot_ids.append(result['primaryAccession'])
            
            logger.info(f"Found {len(uniprot_ids)} UniProt IDs with OMA fingerprint {fingerprint}")
            return uniprot_ids
            
        except requests.RequestException as e:
            logger.error(f"Error searching UniProt for OMA fingerprint {fingerprint}: {e}")
            return []

    def get_oma_group_members_with_valid_uniprot(self, oma_id: str) -> List[Dict]:
        """
        Get OMA group members but only keep those with valid UniProt cross-references.
        """
        try:
            url = f"{self.oma_base_url}group/{oma_id}/"
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            valid_members = []
            
            if 'members' in data:
                for member in data['members']:
                    # Extract potential UniProt ID
                    uniprot_id = member.get('canonicalid', '')
                    
                    # Validate UniProt ID format and check if it exists
                    if uniprot_id and re.match(r'^[A-Z0-9]{6,10}$', uniprot_id):
                        # Verify the UniProt ID exists by making a quick API call
                        if self._validate_uniprot_id(uniprot_id):
                            member_info = {
                                'uniprot_id': uniprot_id,
                                'species': member.get('species', {}).get('name', ''),
                                'taxon_id': member.get('species', {}).get('taxon_id', ''),
                                'kingdom': self._get_kingdom_from_taxon(member.get('species', {}).get('taxon_id', ''))
                            }
                            valid_members.append(member_info)
            
            logger.info(f"Found {len(valid_members)} members with valid UniProt IDs for OMA group {oma_id}")
            return valid_members

        except requests.RequestException as e:
            logger.error(f"Error fetching OMA group {oma_id}: {e}")
            return []

    def _validate_uniprot_id(self, uniprot_id: str) -> bool:
        """
        Validate that a UniProt ID exists by checking the UniProt API.
        """
        try:
            url = f"{self.uniprot_base_url}/uniprotkb/{uniprot_id}"
            response = requests.head(url, timeout=5)  # Use HEAD to avoid downloading full data
            time.sleep(0.05)  # Small delay to be respectful to the API
            return response.status_code == 200
        except requests.RequestException:
            return False

    def _get_kingdom_from_lineage(self, lineage: List[str]) -> str:
        """
        Determine kingdom from UniProt lineage information.
        """
        if not lineage:
            return "Unknown"
        
        lineage_str = ' '.join(lineage).lower()
        
        if 'bacteria' in lineage_str:
            return "Bacteria"
        elif 'archaea' in lineage_str:
            return "Archaea"
        elif any(term in lineage_str for term in ['metazoa', 'animals', 'vertebrata']):
            return "Metazoa"
        elif any(term in lineage_str for term in ['viridiplantae', 'plants', 'plantae']):
            return "Viridiplantae"
        elif 'fungi' in lineage_str:
            return "Fungi"
        elif 'eukaryota' in lineage_str:
            return "Eukaryota"
        else:
            return "Unknown"

    def _get_kingdom_from_taxon(self, taxon_id: str) -> str:
        """
        Determine kingdom from taxon ID patterns with better separation.
        """
        if not taxon_id:
            return "Unknown"
        
        taxon_id = str(taxon_id)
        
        # More specific bacterial taxonomy
        if taxon_id.startswith('2'):
            if taxon_id.startswith('2157'):
                return "Archaea"
            else:
                return "Bacteria"
        
        # Eukaryotic kingdoms with better separation
        elif taxon_id.startswith('33208'):  # Metazoa
            return "Metazoa"
        elif taxon_id.startswith('33090'):  # Viridiplantae
            return "Viridiplantae"
        elif taxon_id.startswith('4751'):   # Fungi
            return "Fungi"
        elif taxon_id.startswith('554915'): # Amoebozoa
            return "Amoebozoa"
        elif taxon_id.startswith('2763'):   # Rhodophyta
            return "Rhodophyta"
        elif taxon_id.startswith('3027'):   # Cryptophyta
            return "Cryptophyta"
        elif taxon_id.startswith('5878'):   # Ciliophora
            return "Ciliophora"
        elif taxon_id.startswith('5747'):   # Apicomplexa
            return "Apicomplexa"
        elif taxon_id.startswith('207245'): # Diplomonadida
            return "Diplomonadida"
        elif taxon_id.startswith('2759'):   # Other Eukaryota
            return "Other Eukaryota"
        else:
            return "Unknown"

    def analyze_pfam_family(self, pfam_folder: str, min_count: int = 3) -> Dict:
        pfam_id = self.get_pfam_id_from_folder(pfam_folder)
        logger.info(f"Starting analysis of Pfam family {pfam_id} from folder {pfam_folder}")

        pfam_uniprot_ids = self.read_pfam_scores_file(pfam_folder)
        if not pfam_uniprot_ids:
            logger.error(f"No UniProt IDs found in scores file from {pfam_folder}")
            return {}

        pfam_uniprot_set = set(pfam_uniprot_ids)
        oma_mapping = self.batch_get_oma_fingerprints(pfam_uniprot_ids)
        if not oma_mapping:
            logger.error("No OMA fingerprints found")
            return {}

        frequent_omas = self.filter_frequent_omas(oma_mapping, min_count=min_count)

        all_oma_uniprot_ids = set()
        oma_details = {}

        for oma_id, count in frequent_omas.items():
            logger.info(f"Processing OMA group {oma_id} (count: {count})")
            
            # Get the OMA fingerprint for this group
            fingerprint = self.get_oma_fingerprint_from_group(oma_id)
            if not fingerprint:
                logger.warning(f"Could not get fingerprint for OMA group {oma_id}, using group ID as fingerprint")
                fingerprint = oma_id
            
            logger.info(f"Using fingerprint '{fingerprint}' for OMA group {oma_id}")
            # Get all UniProt IDs with this fingerprint
            oma_uniprot_ids_list = self.get_uniprot_ids_with_oma_fingerprint(fingerprint)
            
            # Filter out UniProt IDs that are already in the Pfam family
            oma_uniprot_ids = set(oma_uniprot_ids_list) - pfam_uniprot_set

            all_oma_uniprot_ids.update(oma_uniprot_ids)

            oma_details[oma_id] = {
                'count': count,
                'total_members': len(oma_uniprot_ids_list),
                'fingerprint': fingerprint,
                'uniprot_ids': oma_uniprot_ids
            }

            time.sleep(self.request_delay)

        unique_to_oma = all_oma_uniprot_ids
        results = {
            'pfam_id': pfam_id,
            'pfam_folder': pfam_folder,
            'min_count': min_count,
            'pfam_uniprot_count': len(pfam_uniprot_ids),
            'pfam_uniprot_ids': pfam_uniprot_ids,
            'oma_mapping': oma_mapping,
            'frequent_omas': frequent_omas,
            'oma_details': oma_details,
            'unique_to_oma_count': len(unique_to_oma),
            'unique_to_oma_ids': list(unique_to_oma)
        }

        logger.info(f"Analysis complete. Found {len(unique_to_oma)} UniProt IDs unique to OMA")
        return results

    def generate_report(self, results: Dict, output_file: str = "report.txt") -> str:
        if not results:
            return "No results to report."

        lines = []
        lines.append("="*80)
        lines.append(f"PFAM-OMA ORTHOLOG ANALYSIS REPORT")
        lines.append("="*80)
        lines.append(f"Pfam Family: {results['pfam_id']}")
        lines.append(f"Minimum Count Threshold: {results['min_count']}")
        lines.append(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")
        lines.append("SUMMARY STATISTICS")
        lines.append("-" * 50)
        lines.append(f"Total UniProt IDs in Pfam family: {results['pfam_uniprot_count']}")
        lines.append(f"UniProt IDs with OMA fingerprints: {len(results['oma_mapping'])}")
        lines.append(f"Frequent OMA groups (>={results['min_count']} occurrences): {len(results['frequent_omas'])}")
        lines.append(f"UniProt IDs unique to OMA: {results['unique_to_oma_count']}")
        lines.append("")

        lines.append(f"FREQUENT OMA GROUPS (>={results['min_count']} occurrences)")
        lines.append("-" * 50)
        for oma_id, count in results['frequent_omas'].items():
            total_members = results['oma_details'][oma_id]['total_members']
            lines.append(f"{oma_id}: {count} occurrences in Pfam, {total_members} total valid members")
        lines.append("")

        if results['unique_to_oma_ids']:
            lines.append("UNIPROT IDs UNIQUE TO OMA (Not in Pfam)")
            lines.append("-" * 50)
            
            # Group by kingdom for better organization
            kingdom_groups = {}
            
            logger.info("Fetching UniProt descriptions and kingdom information...")
            for oma_id, details in results['oma_details'].items():
                if details['uniprot_ids']:
                    for uniprot_id in details['uniprot_ids']:
                        # Get description
                        description = self.get_uniprot_description(uniprot_id)
                        
                        # For now, we'll use a simplified kingdom classification
                        # You might want to enhance this by fetching organism info from UniProt
                        kingdom = "Eukaryota"  # Default assumption
                        
                        if kingdom not in kingdom_groups:
                            kingdom_groups[kingdom] = []
                        
                        kingdom_groups[kingdom].append({
                            'uniprot_id': uniprot_id,
                            'description': description,
                            'oma_id': oma_id
                        })
            
            # Sort kingdoms for consistent output
            for kingdom in sorted(kingdom_groups.keys()):
                lines.append(f"KINGDOM: {kingdom}")
                lines.append("=" * 40)
                
                # Sort by UniProt ID
                sorted_entries = sorted(kingdom_groups[kingdom], key=lambda x: x['uniprot_id'])
                
                for entry in sorted_entries:
                    lines.append(f"      {entry['uniprot_id']} | {entry['description']} | OMA: {entry['oma_id']}")
                
                lines.append("")
        else:
            lines.append("No UniProt IDs found that are unique to OMA groups.")

        lines.append("="*80)

        report_text = "\n".join(lines)
        
        # Always write to report file
        try:
            with open(output_file, 'w') as f:
                f.write(report_text)
            logger.info(f"Report saved to {output_file}")
        except IOError as e:
            logger.error(f"Failed to write report to {output_file}: {e}")
        
        return report_text

def main():
    parser = argparse.ArgumentParser(description="Analyze Pfam families vs OMA ortholog groups using local Pfam folder")
    parser.add_argument("pfam_folder", help="Path to Pfam family folder containing the 'scores' file")
    parser.add_argument("--output", "-o", default="report.txt", help="Output file for the report (default: report.txt)")
    parser.add_argument("--min-count", "-c", type=int, default=3, help="Minimum count threshold for OMA groups (default: 3)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not os.path.exists(args.pfam_folder) or not os.path.isdir(args.pfam_folder):
        logger.error(f"Pfam folder not found or is not a directory: {args.pfam_folder}")
        sys.exit(1)

    scores_file = os.path.join(args.pfam_folder, "scores")
    if not os.path.exists(scores_file):
        logger.error(f"Scores file not found in Pfam folder: {scores_file}")
        sys.exit(1)

    try:
        analyzer = PfamOMAAnalyzer()
        results = analyzer.analyze_pfam_family(args.pfam_folder, min_count=args.min_count)
        if not results:
            sys.exit(1)

        report = analyzer.generate_report(results, args.output)
        print(f"Analysis complete. Report saved to {args.output}")

    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Analysis failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
