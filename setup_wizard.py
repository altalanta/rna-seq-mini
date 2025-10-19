#!/usr/bin/env python3
"""
Interactive Setup Wizard for RNASEQ-MINI

This script guides users through configuring their RNA-seq analysis pipeline
with intelligent defaults and validation.
"""

import argparse
import os
import sys
import yaml
from pathlib import Path
import subprocess
import shutil


class SetupWizard:
    """Interactive setup wizard for RNASEQ-MINI configuration."""

    def __init__(self):
        self.project_root = Path(__file__).parent.parent
        self.config_dir = self.project_root / "config"
        self.config_dir.mkdir(exist_ok=True)

        # Default configurations
        self.default_params = {
            "project": "my-rnaseq-analysis",
            "engine": "snakemake",
            "threads": min(os.cpu_count() or 4, 8),
            "memory_gb": min(os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024.**3), 16),
            "organism": "human",
            "se": False,
            "fastqc": {"extra": ""},
            "multiqc": {"title": "RNA-seq QC"},
            "salmon": {"libtype": "A", "extra": "--validateMappings --gcBias", "threads": 4},
            "r": {
                "design": "~ condition",
                "contrast_variable": "condition",
                "alpha": 0.05,
                "lfc_shrink": True,
                "contrasts_file": "config/contrasts.tsv",
                "gene_id_column": "gene_id",
                "pvalue_adjust": "BH"
            },
            "fgsea": {"genesets": None, "min_size": 15, "max_size": 500, "nperm": 1000, "padj_cutoff": 0.05},
            "report": {
                "author": self._get_git_user(),
                "title": "RNA-seq Analysis Report",
                "output_html": "results/report.html"
            },
            "cache": {"enabled": True, "dir": ".cache", "max_age_days": 30}
        }

        self.organism_presets = {
            "human": {"transcripts": "references/human/transcripts.fa.gz", "annotation": "references/human/annotation.gtf.gz"},
            "mouse": {"transcripts": "references/mouse/transcripts.fa.gz", "annotation": "references/mouse/annotation.gtf.gz"},
            "yeast": {"transcripts": "references/yeast/transcripts.fa.gz", "annotation": "references/yeast/annotation.gtf.gz"},
            "zebrafish": {"transcripts": "references/zebrafish/transcripts.fa.gz", "annotation": "references/zebrafish/annotation.gtf.gz"}
        }

    def _get_git_user(self):
        """Get git username for default author."""
        try:
            result = subprocess.run(['git', 'config', 'user.name'],
                                  capture_output=True, text=True, cwd=self.project_root)
            return result.stdout.strip() if result.returncode == 0 else "Your Name"
        except:
            return "Your Name"

    def _print_header(self, title):
        """Print a formatted header."""
        print(f"\n{'='*60}")
        print(f"🔬 {title}")
        print(f"{'='*60}")

    def _ask_question(self, question, default=None, options=None):
        """Ask user a question with optional validation."""
        if default:
            prompt = f"{question} [{default}]: "
        else:
            prompt = f"{question}: "

        while True:
            answer = input(prompt).strip()

            if not answer and default:
                return default

            if options and answer not in options:
                print(f"Please choose from: {', '.join(options)}")
                continue

            return answer or default

    def _ask_yes_no(self, question, default=True):
        """Ask a yes/no question."""
        default_str = "(Y/n)" if default else "(y/N)"
        answer = input(f"{question} {default_str}: ").strip().lower()

        if not answer:
            return default
        return answer in ['y', 'yes', 'true', '1']

    def _detect_fastq_files(self):
        """Auto-detect FASTQ files in common locations."""
        common_dirs = ['data', 'fastq', 'raw_data', 'reads']

        for dir_name in common_dirs:
            dir_path = self.project_root / dir_name
            if dir_path.exists():
                fastq_files = list(dir_path.glob("*R1*.fastq*")) + list(dir_path.glob("*_1.fastq*"))
                if fastq_files:
                    return dir_path, fastq_files[:3]  # Return first 3 as examples

        return None, []

    def configure_project(self):
        """Step 1: Basic project configuration."""
        self._print_header("Project Configuration")

        print("Let's start by configuring your analysis project.")

        self.default_params["project"] = self._ask_question(
            "Project name", self.default_params["project"]
        )

        self.default_params["author"] = self._ask_question(
            "Your name (for reports)", self._get_git_user()
        )

        # Engine selection
        engine = self._ask_question(
            "Workflow engine",
            self.default_params["engine"],
            ["snakemake", "nextflow"]
        )
        self.default_params["engine"] = engine

        # Resource configuration
        print(f"\nDetected {os.cpu_count() or 'unknown'} CPU cores")
        threads = int(self._ask_question(
            "Number of threads to use",
            str(self.default_params["threads"])
        ))
        self.default_params["threads"] = threads

        # Auto-detect available memory
        try:
            import psutil
            memory_gb = psutil.virtual_memory().total / (1024**3)
            print(f"Detected {memory_gb:.1f}GB of RAM")
        except ImportError:
            memory_gb = 16

        memory = float(self._ask_question(
            "Memory to allocate (GB)",
            str(min(memory_gb * 0.8, 16))
        ))
        self.default_params["memory_gb"] = memory

    def configure_samples(self):
        """Step 2: Sample configuration."""
        self._print_header("Sample Configuration")

        print("Now let's configure your samples. I can help auto-detect FASTQ files.")

        # Try to auto-detect FASTQ files
        data_dir, sample_files = self._detect_fastq_files()

        if sample_files:
            print(f"\n✅ Found potential FASTQ files in {data_dir}/")
            for f in sample_files:
                print(f"   {f.name}")

        # Ask about sample file location
        samples_file = self._ask_question(
            "Path to samples.tsv file (or I'll help you create one)",
            "config/samples.tsv"
        )

        if not Path(samples_file).exists():
            print("\n📋 Let's create a samples file. I need to know about your samples.")
            print("Expected format: sample<tab>condition<tab>fastq_1<tab>fastq_2")

            samples_data = []
            print("\nExample entry:")
            print("sample1	treated	data/sample1_R1.fastq.gz	data/sample1_R2.fastq.gz")

            while True:
                sample = self._ask_question("Sample name (or 'done' to finish)")
                if sample.lower() == 'done':
                    break

                condition = self._ask_question(f"Condition for {sample}")
                fastq_1 = self._ask_question(f"Forward reads file for {sample}")

                # Check if paired-end
                is_paired = self._ask_yes_no(f"Is {sample} paired-end?", True)
                fastq_2 = ""
                if is_paired:
                    fastq_2 = self._ask_question(f"Reverse reads file for {sample}")

                batch = self._ask_question(f"Batch for {sample} (optional)", "")

                samples_data.append([sample, condition, fastq_1, fastq_2, batch])

            # Write samples file
            with open(samples_file, 'w') as f:
                f.write("sample\tcondition\tfastq_1\tfastq_2\tbatch\n")
                for sample_data in samples_data:
                    f.write("\t".join(sample_data) + "\n")

            print(f"✅ Created {samples_file}")

        self.default_params["paths"] = {"samples": samples_file}

    def configure_contrasts(self):
        """Step 3: Contrast configuration."""
        self._print_header("Contrast Configuration")

        print("Now let's define the comparisons you want to make.")
        print("Expected format: groupA<tab>groupB")

        contrasts_file = self._ask_question(
            "Path to contrasts.tsv file (or I'll help you create one)",
            "config/contrasts.tsv"
        )

        if not Path(contrasts_file).exists():
            print("\n📋 Let's create a contrasts file.")
            print("Example entry:")
            print("treated	control")

            contrasts_data = []
            while True:
                group_a = self._ask_question("First condition (or 'done' to finish)")
                if group_a.lower() == 'done':
                    break

                group_b = self._ask_question(f"Second condition to compare {group_a} against")
                contrasts_data.append([group_a, group_b])

            # Write contrasts file
            with open(contrasts_file, 'w') as f:
                f.write("groupA\tgroupB\n")
                for contrast in contrasts_data:
                    f.write("\t".join(contrast) + "\n")

            print(f"✅ Created {contrasts_file}")

        self.default_params["r"]["contrasts_file"] = contrasts_file

    def configure_organism(self):
        """Step 4: Organism and reference configuration."""
        self._print_header("Organism & Reference Configuration")

        print("Choose your organism (this determines which reference genome to use):")
        organisms = list(self.organism_presets.keys())

        for i, org in enumerate(organisms, 1):
            print(f"{i}. {org}")

        while True:
            choice = self._ask_question(
                "Select organism (number or name)",
                "1"
            )

            try:
                if choice.isdigit():
                    idx = int(choice) - 1
                    if 0 <= idx < len(organisms):
                        organism = organisms[idx]
                        break
                elif choice in organisms:
                    organism = choice
                    break
                else:
                    print(f"Please choose from: {', '.join([str(i) for i in range(1, len(organisms)+1)] + organisms)}")
            except:
                print("Invalid choice. Please try again.")

        self.default_params["organism"] = organism

        # Set reference paths
        if organism in self.organism_presets:
            refs = self.organism_presets[organism]
            self.default_params["reference"] = {
                "transcripts_fa": refs["transcripts"],
                "annotation_gtf": refs["annotation"],
                "decoy_fasta": refs.get("decoy_fasta", ""),
                "salmon_index": refs.get("salmon_index", f"references/{organism}/salmon_index"),
                "auto_download": True
            }

        print(f"✅ Configured for {organism}")

    def configure_analysis_options(self):
        """Step 5: Advanced analysis options."""
        self._print_header("Analysis Options")

        print("Configure advanced analysis options:")

        # Single-end vs paired-end
        se = not self._ask_yes_no(
            "Do you have paired-end reads?",
            True
        )
        self.default_params["se"] = se

        # Salmon library type
        print("\nSalmon library type detection:")
        print("A = Auto-detect (recommended for most cases)")
        print("ISR = Inward, Stranded, Reverse")
        print("ISF = Inward, Stranded, Forward")
        print("SR = Stranded, Reverse")
        print("SF = Stranded, Forward")

        libtype = self._ask_question(
            "Salmon library type",
            self.default_params["salmon"]["libtype"]
        )
        self.default_params["salmon"]["libtype"] = libtype

        # Statistical design
        print("\nStatistical design formula for DESeq2:")
        print("~ condition                    # Simple two-group comparison")
        print("~ condition + batch           # With batch correction")
        print("~ condition + batch + covariate  # Multiple factors")

        design = self._ask_question(
            "Design formula",
            self.default_params["r"]["design"]
        )
        self.default_params["r"]["design"] = design

        # Significance threshold
        alpha = float(self._ask_question(
            "Significance threshold (alpha)",
            str(self.default_params["r"]["alpha"])
        ))
        self.default_params["r"]["alpha"] = alpha

    def configure_advanced_features(self):
        """Step 6: Advanced features."""
        self._print_header("Advanced Features")

        print("Configure optional advanced features:")

        # Caching
        enable_cache = self._ask_yes_no(
            "Enable intelligent caching? (recommended for development)",
            True
        )
        self.default_params["cache"]["enabled"] = enable_cache

        # MultiQC title
        title = self._ask_question(
            "MultiQC report title",
            f"RNA-seq QC - {self.default_params['project']}"
        )
        self.default_params["multiqc"]["title"] = title

        # Report title
        report_title = self._ask_question(
            "Analysis report title",
            f"{self.default_params['project']} - RNA-seq Analysis"
        )
        self.default_params["report"]["title"] = report_title

    def generate_configuration(self):
        """Generate final configuration files."""
        self._print_header("Generating Configuration")

        print("Generating configuration files...")

        # Create params.yaml
        params_file = self.config_dir / "params.yaml"
        with open(params_file, 'w') as f:
            yaml.dump(self.default_params, f, default_flow_style=False, sort_keys=False)

        print(f"✅ Created {params_file}")

        # Copy sample and contrast files if they exist
        for filename in ["samples.tsv", "contrasts.tsv"]:
            src = self.config_dir / filename
            dst = self.config_dir / f"params_{filename}"
            if src.exists():
                shutil.copy2(src, dst)
                print(f"✅ Backed up existing {src} to {dst}")

        print("\n🎉 Configuration complete!")
        print(f"\nProject: {self.default_params['project']}")
        print(f"Organism: {self.default_params['organism']}")
        print(f"Engine: {self.default_params['engine']}")
        print(f"Samples file: {self.default_params['paths']['samples']}")
        print(f"Contrasts file: {self.default_params['r']['contrasts_file']}")

        print("\n🚀 Next steps:")
        print("1. Run 'make run' to start the analysis")
        print("2. Monitor progress with 'make monitor'")
        print("3. View results with 'make serve-results'")
        print("4. Check 'docs/workflow.md' for detailed workflow information")

    def run(self):
        """Run the complete setup wizard."""
        print("🔬 RNASEQ-MINI Setup Wizard")
        print("=" * 60)
        print("This wizard will help you configure your RNA-seq analysis pipeline.")
        print("You can accept defaults by pressing Enter, or customize any setting.")

        if not self._ask_yes_no("Ready to begin?", True):
            print("Setup cancelled.")
            return

        try:
            self.configure_project()
            self.configure_samples()
            self.configure_contrasts()
            self.configure_organism()
            self.configure_analysis_options()
            self.configure_advanced_features()
            self.generate_configuration()

        except KeyboardInterrupt:
            print("\n\nSetup interrupted by user.")
            return
        except Exception as e:
            print(f"\n❌ Setup failed: {e}")
            return


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="RNASEQ-MINI Setup Wizard")
    parser.add_argument("--reset", action="store_true",
                       help="Reset configuration and start fresh")

    args = parser.parse_args()

    wizard = SetupWizard()

    if args.reset:
        print("🔄 Resetting configuration...")
        config_files = ["params.yaml", "samples.tsv", "contrasts.tsv"]
        for filename in config_files:
            filepath = wizard.config_dir / filename
            if filepath.exists():
                backup = wizard.config_dir / f"{filename}.backup"
                shutil.move(filepath, backup)
                print(f"  Backed up {filepath} to {backup}")

    wizard.run()


if __name__ == "__main__":
    main()
