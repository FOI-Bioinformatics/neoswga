"""
Interactive workflow selector for neoswga.

Provides a user-friendly menu to discover and launch different
neoswga features:
1. Standard pipeline setup and execution
2. Multi-genome design
3. Analysis tools
4. Simulation
"""

import sys
import subprocess
import shlex
from typing import Optional, List, Tuple


def clear_screen():
    """Clear terminal screen."""
    print("\033[H\033[J", end="")


def execute_command(cmd: List[str]) -> bool:
    """
    Execute a neoswga command.

    Args:
        cmd: Command as a list of arguments

    Returns:
        True if command succeeded, False otherwise
    """
    display_cmd = " ".join(shlex.quote(c) for c in cmd)
    print(f"\nExecuting: {display_cmd}\n")
    print("-" * 60)
    try:
        result = subprocess.run(cmd)
        print("-" * 60)
        if result.returncode == 0:
            print("\nCommand completed successfully.")
            return True
        else:
            print(f"\nCommand exited with code {result.returncode}")
            return False
    except FileNotFoundError:
        print("Error: neoswga command not found. Is it installed?")
        return False
    except Exception as e:
        print(f"Error executing command: {e}")
        return False


def prompt_execute(cmd: List[str], description: str = "") -> bool:
    """
    Show a command and ask user if they want to execute it.

    Args:
        cmd: Command as a list of arguments
        description: What the command does

    Returns:
        True if command was executed and succeeded
    """
    display_cmd = " ".join(shlex.quote(c) for c in cmd)
    print()
    if description:
        print(f"  {description}")
        print()
    print(f"  Command: {display_cmd}")
    print()

    while True:
        choice = input("Run this command? [y/n]: ").strip().lower()
        if choice == 'y' or choice == 'yes':
            return execute_command(cmd)
        elif choice == 'n' or choice == 'no':
            print("Skipped.")
            return False
        else:
            print("Please enter 'y' or 'n'")


def print_menu(title: str, options: List[Tuple[str, str]], show_quit: bool = True):
    """
    Print a numbered menu.

    Args:
        title: Menu title
        options: List of (label, description) tuples
        show_quit: Whether to show quit option
    """
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)
    print()

    for i, (label, description) in enumerate(options, 1):
        print(f"  {i}. {label}")
        if description:
            print(f"     {description}")
        print()

    if show_quit:
        print(f"  q. Quit")
        print()


def get_choice(max_option: int) -> Optional[int]:
    """
    Get user's menu choice.

    Args:
        max_option: Maximum valid option number

    Returns:
        Selected option number, or None for quit
    """
    while True:
        choice = input("Select option: ").strip().lower()

        if choice == 'q':
            return None

        try:
            num = int(choice)
            if 1 <= num <= max_option:
                return num
            print(f"Please enter a number between 1 and {max_option}")
        except ValueError:
            print("Please enter a valid number or 'q' to quit")


def run_workflow_selector():
    """
    Run the interactive workflow selector.

    Displays menus and launches appropriate commands.
    """
    main_options = [
        ("Set up new project", "Create params.json with guided configuration"),
        ("Run primer design pipeline", "Execute count-kmers, filter, score, optimize"),
        ("Validate configuration", "Check params.json for errors before running"),
        ("Interpret results", "Get quality assessment of primer design output"),
        ("Advanced features", "Multi-genome design, simulation, analysis tools"),
    ]

    advanced_options = [
        ("Multi-genome design", "Design primers for multiple targets with differential penalties"),
        ("Simulate amplification", "Run agent-based replication simulation"),
        ("Analyze genome", "Check genome suitability for SWGA"),
        ("Analyze primer set", "Evaluate existing primer sequences"),
        ("Analyze dimers", "Visualize primer-dimer interaction network"),
        ("Active learning", "Iterative optimization with experimental feedback"),
        ("Back to main menu", ""),
    ]

    pipeline_options = [
        ("Run full pipeline", "Execute all steps: count-kmers -> filter -> score -> optimize"),
        ("Step 1: Count k-mers", "Generate k-mer counts from genome"),
        ("Step 2: Filter primers", "Apply frequency and thermodynamic filters"),
        ("Step 3: Score candidates", "Predict amplification efficacy"),
        ("Step 4: Optimize set", "Select optimal primer combination"),
        ("Back to main menu", ""),
    ]

    while True:
        print_menu("NEOSWGA WORKFLOW SELECTOR", main_options)
        choice = get_choice(len(main_options))

        if choice is None:
            print("\nGoodbye!")
            sys.exit(0)

        elif choice == 1:  # Set up new project
            print("\n" + "-" * 40)
            print("SET UP NEW PROJECT")
            print("-" * 40)
            print()
            print("This will analyze your genome and create a params.json")
            print("configuration file with recommended settings.")
            print()

            # Prompt for genome path
            genome = input("Enter path to target genome FASTA (or 'skip' to exit): ").strip()
            if genome.lower() != 'skip' and genome:
                background = input("Enter path to background genome FASTA (or Enter to skip): ").strip()
                cmd = ["neoswga", "init", "--genome", genome]
                if background:
                    cmd.extend(["--background", background])
                prompt_execute(cmd, "Create params.json with recommended settings")
            else:
                print("\nTo set up later, run:")
                print("  neoswga init --genome <target.fasta> [--background <bg.fasta>]")

            input("\nPress Enter to continue...")

        elif choice == 2:  # Run pipeline
            while True:
                print_menu("PRIMER DESIGN PIPELINE", pipeline_options)
                pipe_choice = get_choice(len(pipeline_options))

                if pipe_choice is None or pipe_choice == len(pipeline_options):
                    break

                # Get params file path (used by all pipeline steps)
                params_file = input("Enter params.json path [params.json]: ").strip()
                if not params_file:
                    params_file = "params.json"

                if pipe_choice == 1:  # Full pipeline
                    print("\nRunning full pipeline (4 steps)...")
                    steps = [
                        (["neoswga", "count-kmers", "-j"], "Step 1: Count k-mers"),
                        (["neoswga", "filter", "-j"], "Step 2: Filter primers"),
                        (["neoswga", "score", "-j"], "Step 3: Score candidates"),
                        (["neoswga", "optimize", "-j"], "Step 4: Optimize set"),
                    ]
                    for base_cmd, desc in steps:
                        cmd = base_cmd + [params_file]
                        if not prompt_execute(cmd, desc):
                            print("\nPipeline stopped due to error or user choice.")
                            break
                elif pipe_choice == 2:
                    prompt_execute(["neoswga", "count-kmers", "-j", params_file], "Generate k-mer counts from genome")
                elif pipe_choice == 3:
                    prompt_execute(["neoswga", "filter", "-j", params_file], "Apply frequency and thermodynamic filters")
                elif pipe_choice == 4:
                    prompt_execute(["neoswga", "score", "-j", params_file], "Predict amplification efficacy")
                elif pipe_choice == 5:
                    print("\nOptimization methods available:")
                    print("  1. hybrid (default)")
                    print("  2. dominating-set (8x faster)")
                    print("  3. background-aware (clinical, 10-20x bg reduction)")
                    method = input("Select method [1]: ").strip()
                    methods = {'1': 'hybrid', '2': 'dominating-set', '3': 'background-aware'}
                    opt_method = methods.get(method, 'hybrid')
                    cmd = ["neoswga", "optimize", "-j", params_file, f"--optimization-method={opt_method}"]
                    prompt_execute(cmd, "Select optimal primer combination")

                input("\nPress Enter to continue...")

        elif choice == 3:  # Validate configuration
            print("\n" + "-" * 40)
            print("VALIDATE CONFIGURATION")
            print("-" * 40)
            print()
            print("This checks for:")
            print("  - Missing required parameters")
            print("  - Invalid parameter values")
            print("  - Incompatible settings (e.g., phi29 + high temperature)")
            print("  - Missing genome files")
            print()

            params_file = input("Enter params.json path [params.json]: ").strip()
            if not params_file:
                params_file = "params.json"
            prompt_execute(["neoswga", "validate-params", "-j", params_file], "Check configuration for errors")
            input("\nPress Enter to continue...")

        elif choice == 4:  # Interpret results
            print("\n" + "-" * 40)
            print("INTERPRET RESULTS")
            print("-" * 40)
            print()
            print("This provides:")
            print("  - Quality ratings for coverage, enrichment, uniformity")
            print("  - Go/no-go recommendation for synthesis")
            print("  - Suggested next steps")
            print()

            results_dir = input("Enter results directory [results/]: ").strip()
            if not results_dir:
                results_dir = "results/"
            prompt_execute(["neoswga", "interpret", "-d", results_dir], "Quality assessment of primer design output")
            input("\nPress Enter to continue...")

        elif choice == 5:  # Advanced features
            while True:
                print_menu("ADVANCED FEATURES", advanced_options)
                adv_choice = get_choice(len(advanced_options))

                if adv_choice is None or adv_choice == len(advanced_options):
                    break

                if adv_choice == 1:  # Multi-genome
                    print("\nMulti-genome primer design:")
                    print()
                    print("  neoswga multi-genome \\")
                    print("    --genomes target1.fasta target2.fasta \\")
                    print("    --background host.fasta \\")
                    print("    --blacklist pathogen2.fasta \\")
                    print("    --output results/")
                    print()
                    print("Genome roles:")
                    print("  --genomes    Targets to amplify (maximize binding)")
                    print("  --background Off-target to tolerate (e.g., host DNA)")
                    print("  --blacklist  Off-target to avoid (e.g., competing pathogens)")

                elif adv_choice == 2:  # Simulate
                    print("\nReplication simulation:")
                    print()
                    print("  neoswga simulate \\")
                    print("    --primers SEQ1 SEQ2 SEQ3 \\")
                    print("    --genome target.fasta \\")
                    print("    --output sim/")
                    print()
                    print("Validates primer set performance before synthesis")

                elif adv_choice == 3:  # Analyze genome
                    print("\nGenome suitability analysis:")
                    print()
                    print("  neoswga analyze-genome \\")
                    print("    --genome target.fasta \\")
                    print("    --output analysis/")
                    print()
                    print("Assesses GC content, complexity, and SWGA suitability")

                elif adv_choice == 4:  # Analyze primer set
                    print("\nPrimer set analysis:")
                    print()
                    print("  neoswga analyze-set \\")
                    print("    --primers SEQ1 SEQ2 \\")
                    print("    --fg target.fasta \\")
                    print("    --fg-kmers data/target \\")
                    print("    --output analysis/")

                elif adv_choice == 5:  # Analyze dimers
                    print("\nDimer network analysis:")
                    print()
                    print("  neoswga analyze-dimers \\")
                    print("    --primers SEQ1 SEQ2 SEQ3 \\")
                    print("    --output dimers/ \\")
                    print("    --visualize")
                    print()
                    print("Creates network visualization of primer-dimer interactions")

                elif adv_choice == 6:  # Active learning
                    print("\nActive learning (experimental):")
                    print()
                    print("  neoswga active-learn \\")
                    print("    -j params.json \\")
                    print("    --output active_learn/ \\")
                    print("    --num-candidates 10")
                    print()
                    print("Iterative optimization with experimental feedback")

                input("\nPress Enter to continue...")


if __name__ == '__main__':
    run_workflow_selector()
