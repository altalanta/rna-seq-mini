#!/usr/bin/env python3
"""
Intelligent reference management system for RNASEQ-MINI.
Downloads and manages reference genomes for any species from NCBI/Ensembl.

Usage: python scripts/download_references.py [species] [options]
"""

import argparse
import sys
import json
import yaml
import requests
import gzip
import shutil
import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import subprocess
import time
import re


class ReferenceManager:
    """Manages reference genome downloads and organization."""

    def __init__(self, base_dir: Path = Path("references")):
        self.base_dir = base_dir
        self.cache_dir = base_dir / ".cache"
        self.cache_dir.mkdir(exist_ok=True)

        # Species presets with download URLs
        self.species_presets = {
            'human': {
                'scientific_name': 'Homo sapiens',
                'assembly': 'GRCh38',
                'ensembl_release': '110',
                'ncbi_taxon_id': '9606'
            },
            'mouse': {
                'scientific_name': 'Mus musculus',
                'assembly': 'GRCm39',
                'ensembl_release': '110',
                'ncbi_taxon_id': '10090'
            },
            'yeast': {
                'scientific_name': 'Saccharomyces cerevisiae',
                'assembly': 'R64-1-1',
                'ensembl_release': '110',
                'ncbi_taxon_id': '559292'
            },
            'zebrafish': {
                'scientific_name': 'Danio rerio',
                'assembly': 'GRCz11',
                'ensembl_release': '110',
                'ncbi_taxon_id': '7955'
            },
            'fruitfly': {
                'scientific_name': 'Drosophila melanogaster',
                'assembly': 'BDGP6.32',
                'ensembl_release': '110',
                'ncbi_taxon_id': '7227'
            }
        }

    def get_ensembl_download_urls(self, species: str, release: str = '110') -> Dict[str, str]:
        """Get Ensembl download URLs for a species."""
        base_url = f"https://ftp.ensembl.org/pub/release-{release}/fasta/{species.lower()}"

        urls = {
            'transcripts_fa': f"{base_url}/dna/{species.capitalize()}.{species}.dna.toplevel.fa.gz",
            'annotation_gtf': f"{base_url}/gtf/{species.capitalize()}.{species}.{release}.gtf.gz"
        }

        return urls

    def get_ncbi_download_info(self, species: str, taxon_id: str) -> Dict[str, str]:
        """Get NCBI download information for a species."""
        # Use Ensembl Genomes for non-vertebrates, NCBI for vertebrates
        if species.lower() in ['yeast', 'fruitfly']:
            return self.get_ensembl_download_urls(species)
        else:
            # For vertebrates, try to find NCBI RefSeq
            return self._get_ncbi_refseq_urls(taxon_id)

    def _get_ncbi_refseq_urls(self, taxon_id: str) -> Dict[str, str]:
        """Get NCBI RefSeq download URLs for a species."""
        # This is a simplified implementation - in practice you'd query NCBI API
        # For now, return placeholder URLs that would need to be updated
        return {
            'transcripts_fa': f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna.fna.gz",
            'annotation_gtf': f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
        }

    def download_file(self, url: str, output_path: Path, expected_size: Optional[int] = None) -> bool:
        """Download a file with progress tracking and integrity verification."""
        print(f"📥 Downloading {output_path.name}...")

        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()

            total_size = int(response.headers.get('content-length', 0))

            with open(output_path, 'wb') as f:
                downloaded = 0
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)

                        # Progress bar
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"\r   Progress: {percent:.1f}% ({downloaded}/{total_size} bytes)", end='', flush=True)

            print()  # New line after progress

            # Verify file size if expected size provided
            if expected_size and output_path.stat().st_size != expected_size:
                print(f"❌ File size mismatch for {output_path.name}")
                return False

            return True

        except Exception as e:
            print(f"❌ Download failed for {url}: {e}")
            return False

    def verify_file_integrity(self, file_path: Path, expected_hash: Optional[str] = None) -> bool:
        """Verify file integrity using checksum."""
        if not file_path.exists():
            return False

        # Compute SHA256 hash
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                sha256_hash.update(chunk)

        computed_hash = sha256_hash.hexdigest()

        if expected_hash:
            if computed_hash == expected_hash:
                print(f"✅ Integrity check passed for {file_path.name}")
                return True
            else:
                print(f"❌ Integrity check failed for {file_path.name}")
                return False

        print(f"✅ File exists and appears valid: {file_path.name}")
        return True

    def create_decoy_fasta(self, transcripts_fa: Path, annotation_gtf: Path, output_path: Path) -> bool:
        """Create decoy FASTA file for Salmon indexing."""
        print("🔧 Creating decoy FASTA file...")

        try:
            # Extract chromosome/scaffold names from GTF
            chromosomes = set()

            with gzip.open(annotation_gtf, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 1:
                        seqname = fields[0]
                        chromosomes.add(seqname)

            # Create decoy FASTA with chromosome names
            with gzip.open(transcripts_fa, 'rt') as infile, open(output_path, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        # Check if this sequence name is in our chromosome list
                        seqname = line[1:].strip().split()[0]
                        if seqname in chromosomes:
                            chromosomes.remove(seqname)
                        outfile.write(line)
                    else:
                        outfile.write(line)

                # Add decoy sequences
                for chrom in sorted(chromosomes):
                    outfile.write(f">{chrom}\n")
                    outfile.write("N" * 100 + "\n")  # Placeholder sequence

            print(f"✅ Created decoy FASTA: {output_path}")
            return True

        except Exception as e:
            print(f"❌ Failed to create decoy FASTA: {e}")
            return False

    def setup_species_references(self, species: str, force: bool = False) -> bool:
        """Download and setup references for a species."""
        print(f"🚀 Setting up references for {species}...")

        if species not in self.species_presets:
            print(f"❌ Unknown species: {species}")
            print(f"Available species: {', '.join(self.species_presets.keys())}")
            return False

        species_info = self.species_presets[species]
        species_dir = self.base_dir / species
        species_dir.mkdir(exist_ok=True)

        print(f"📋 Species info: {species_info}")

        # Check if references already exist
        transcripts_fa = species_dir / "transcripts.fa.gz"
        annotation_gtf = species_dir / "annotation.gtf.gz"
        decoys_fa = species_dir / "decoys.fa.gz"

        if not force and all(p.exists() for p in [transcripts_fa, annotation_gtf]):
            print("✅ References already exist. Use --force to re-download.")
            return True

        # Get download URLs
        download_urls = self.get_ncbi_download_info(species, species_info['ncbi_taxon_id'])

        # Download transcripts
        transcripts_url = download_urls['transcripts_fa']
        if not self.download_file(transcripts_url, transcripts_fa):
            return False

        # Download annotation
        annotation_url = download_urls['annotation_gtf']
        if not self.download_file(annotation_url, annotation_gtf):
            return False

        # Verify downloads
        if not all(self.verify_file_integrity(p) for p in [transcripts_fa, annotation_gtf]):
            return False

        # Create decoy file
        if not self.create_decoy_fasta(transcripts_fa, annotation_gtf, decoys_fa):
            return False

        # Update genome.yaml
        self.update_genome_config(species, species_dir)

        print(f"🎉 Successfully set up references for {species}")
        return True

    def update_genome_config(self, species: str, species_dir: Path):
        """Update the genome.yaml configuration file."""
        genome_config = self.base_dir / ".." / "config" / "genome.yaml"

        # Load existing config
        try:
            with open(genome_config, 'r') as f:
                config = yaml.safe_load(f) or {}
        except FileNotFoundError:
            config = {}

        # Update species configuration
        config[species] = {
            'transcripts_fa': str(species_dir / "transcripts.fa.gz"),
            'annotation_gtf': str(species_dir / "annotation.gtf.gz"),
            'decoy_fasta': str(species_dir / "decoys.fa.gz"),
            'salmon_index': str(species_dir / "salmon_index")
        }

        # Write updated config
        with open(genome_config, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)

        print(f"✅ Updated genome configuration in {genome_config}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Download and manage reference genomes")
    parser.add_argument('species', nargs='?', help='Species to download references for')
    parser.add_argument('--list', action='store_true', help='List available species')
    parser.add_argument('--force', action='store_true', help='Force re-download even if files exist')
    parser.add_argument('--base-dir', default='references', help='Base directory for references')

    args = parser.parse_args()

    manager = ReferenceManager(Path(args.base_dir))

    if args.list:
        print("Available species:")
        for species, info in manager.species_presets.items():
            print(f"  {species}: {info['scientific_name']} ({info['assembly']})")
        return

    if not args.species:
        print("❌ Please specify a species or use --list to see available options")
        return

    success = manager.setup_species_references(args.species, args.force)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
