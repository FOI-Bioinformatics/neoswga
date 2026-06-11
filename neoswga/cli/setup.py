"""Setup / validation / guidance CLI handlers: validate, doctor, init,
validate-params, schema, validate-model, interpret. Extracted from cli_unified.
"""

import json
import logging
import os
import sys

logger = logging.getLogger(__name__)


def run_validate(args):
    """Dispatcher for `neoswga validate [install|params|model]`.

    - No subcommand or `install`: installation / environment check
      (historical behavior of `neoswga validate`).
    - `params`: delegate to :func:`run_validate_params`.
    - `model`: delegate to :func:`run_validate_model`.
    """
    mode = getattr(args, "validate_mode", None)
    if mode == "params":
        return run_validate_params(args)
    if mode == "model":
        return run_validate_model(args)
    # mode is None (bare `validate`) or 'install': installation check.
    return _run_validate_installation(args)


def _run_validate_installation(args):
    """Run installation / environment validation tests."""
    from neoswga.core.gpu_acceleration import get_gpu_info, is_gpu_available
    from neoswga.core.kmer_counter import check_jellyfish_available
    from neoswga.core.validation import ValidationSuite, quick_validation

    # Show system capabilities
    logger.info("System Capabilities:")
    logger.info(f"  Jellyfish: {'Available' if check_jellyfish_available() else 'NOT FOUND'}")
    if is_gpu_available():
        gpu_info = get_gpu_info()
        logger.info(
            f"  GPU: {gpu_info.get('device_name', 'Available')} ({gpu_info.get('memory_total_gb', 0):.1f} GB)"
        )
    else:
        logger.info("  GPU: Not available (install CuPy for GPU acceleration)")

    if args.quick:
        logger.info("Running quick validation...")
        success = quick_validation()
    elif args.all:
        logger.info("Running full validation + benchmarks...")
        suite = ValidationSuite(verbose=True)
        success = suite.run_all_tests()

        # Also run benchmarks if available
        try:
            import importlib.util

            bench_path = os.path.join(
                os.path.dirname(__file__),
                "..",
                "scripts",
                "benchmarking",
                "benchmark_improvements.py",
            )
            spec = importlib.util.spec_from_file_location("benchmark_improvements", bench_path)
            if spec and spec.loader:
                logger.info("\nRunning benchmarks...")
                benchmark_improvements = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(benchmark_improvements)
                benchmark_improvements.main()
            else:
                logger.info("\nBenchmark script not found, skipping.")
        except Exception as e:
            logger.warning(f"Benchmarks unavailable: {e}")
    else:
        logger.info("Running standard validation...")
        suite = ValidationSuite(verbose=True)
        success = suite.run_all_tests()

    sys.exit(0 if success else 1)


def run_doctor(args):
    """Diagnostic health check.

    Reports version, jellyfish location, registered optimizers with
    additive-awareness, params.json summary, and file-existence status.
    Exits non-zero if any blocking issue is found.
    """
    import json as _json
    import os as _os
    import platform as _platform
    import shutil as _shutil
    import subprocess as _subprocess
    import sys as _sys

    blocking: list = []
    warnings: list = []
    summary: dict = {
        "platform": {
            "system": _platform.system(),
            "python": _sys.version.split()[0],
            "neoswga": _import_neoswga_version(),
        },
        "tools": {},
        "optimizers": [],
        "params": None,
    }

    # Jellyfish
    jelly = _shutil.which("jellyfish")
    if jelly:
        try:
            version = _subprocess.run(
                [jelly, "--version"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            first = (version.stdout or version.stderr or "").splitlines()
            summary["tools"]["jellyfish"] = {
                "path": jelly,
                "version": first[0] if first else "unknown",
            }
        except Exception as e:
            summary["tools"]["jellyfish"] = {"path": jelly, "error": str(e)}
            warnings.append(f"jellyfish --version failed: {e}")
    else:
        blocking.append(
            "jellyfish is not on PATH. Install with: conda install -c bioconda jellyfish"
        )

    # Optimizer capability matrix. Each entry declares whether the optimizer
    # is registered with the factory and whether its scoring path meaningfully
    # references `self.conditions`. The awareness flag is computed by grepping
    # the source file for ReactionConditions usage plus `self.conditions`.
    try:
        from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry

        # Force registration of built-in optimizers
        try:
            from neoswga.core import unified_optimizer as _uo

            _uo._ensure_optimizers_registered()
        except Exception:
            pass
        registered = OptimizerFactory.list_optimizers()
        for name in sorted(registered.keys()):
            cls = OptimizerRegistry.get(name)
            # Phase 15D: structural class-level flag replaces the old
            # source-grep heuristic. Subclasses that actually use
            # ReactionConditions in selection opt-in by setting
            # ADDITIVE_AWARE = True on the class.
            aware = bool(getattr(cls, "ADDITIVE_AWARE", False))
            summary["optimizers"].append(
                {
                    "name": name,
                    "class": f"{cls.__module__}.{cls.__name__}",
                    "additive_aware": aware,
                    "description": registered[name],
                }
            )
    except Exception as e:
        warnings.append(f"Optimizer discovery failed: {e}")

    # params.json summary (optional)
    if args.json_file:
        if not _os.path.isfile(args.json_file):
            blocking.append(f"params.json not found at {args.json_file}")
        else:
            try:
                with open(args.json_file) as fh:
                    raw = _json.load(fh)
                params_summary = {
                    "path": args.json_file,
                    "schema_version": raw.get("schema_version"),
                    "polymerase": raw.get("polymerase"),
                    "reaction_temp": raw.get("reaction_temp"),
                    "mg_conc": raw.get("mg_conc"),
                    "additives": {
                        k: raw.get(k, 0.0)
                        for k in (
                            "dmso_percent",
                            "betaine_m",
                            "trehalose_m",
                            "formamide_percent",
                            "ethanol_percent",
                            "urea_m",
                            "tmac_m",
                        )
                    },
                    "min_k": raw.get("min_k"),
                    "max_k": raw.get("max_k"),
                    "num_targets": len(raw.get("fg_genomes", []) or []),
                    "num_backgrounds": len(raw.get("bg_genomes", []) or []),
                    "num_blacklists": len(raw.get("bl_genomes", []) or []),
                    "adaptive_gc": raw.get("adaptive_gc", True),
                    "genome_gc": raw.get("genome_gc"),
                }
                # File existence checks for each prefix
                missing_files: list = []
                for prefix in raw.get("fg_prefixes", []) or []:
                    mn = raw.get("min_k", 6)
                    mx = raw.get("max_k", 12)
                    for k in range(mn, mx + 1):
                        candidate = f"{prefix}_{k}mer_all.txt"
                        if not _os.path.isfile(candidate):
                            missing_files.append(candidate)
                            if len(missing_files) > 10:
                                break
                    if len(missing_files) > 10:
                        break
                if missing_files:
                    warnings.append(
                        "Missing k-mer count files (run count-kmers?): "
                        + ", ".join(missing_files[:5])
                        + ("…" if len(missing_files) > 5 else "")
                    )
                    params_summary["missing_kmer_files"] = missing_files[:20]
                summary["params"] = params_summary
            except Exception as e:
                blocking.append(f"Could not parse params.json: {e}")

    summary["ok"] = not blocking
    summary["blocking"] = blocking
    summary["warnings"] = warnings

    if args.output_json:
        print(_json.dumps(summary, indent=2))
    else:
        # Human-readable report
        print("=" * 60)
        print("neoswga doctor")
        print("=" * 60)
        print(f"  platform: {summary['platform']['system']}")
        print(f"  python: {summary['platform']['python']}")
        print(f"  neoswga: {summary['platform']['neoswga']}")
        j = summary["tools"].get("jellyfish")
        if j:
            print(f"  jellyfish: {j.get('version', j.get('error', '?'))} at {j.get('path', '')}")
        else:
            print("  jellyfish: NOT FOUND")
        print()
        print(f"Registered optimizers: {len(summary['optimizers'])}")
        aware_count = sum(1 for o in summary["optimizers"] if o["additive_aware"])
        print(f"  additive-aware: {aware_count}")
        print(f"  additive-blind: {len(summary['optimizers']) - aware_count}")
        for opt in summary["optimizers"]:
            marker = "+" if opt["additive_aware"] else "-"
            print(f"  [{marker}] {opt['name']}")
        print()
        if summary["params"]:
            ps = summary["params"]
            print("params.json summary:")
            print(f"  path: {ps['path']}")
            print(f"  schema_version: {ps['schema_version']}")
            print(f"  polymerase: {ps['polymerase']}  reaction_temp: {ps['reaction_temp']}")
            print(f"  min_k: {ps['min_k']}  max_k: {ps['max_k']}")
            print(
                f"  targets: {ps['num_targets']}  backgrounds: {ps['num_backgrounds']}  blacklists: {ps['num_blacklists']}"
            )
            nz_add = {k: v for k, v in ps["additives"].items() if v}
            if nz_add:
                print(f"  additives: {nz_add}")
            print()
        if warnings:
            print("Warnings:")
            for w in warnings:
                print(f"  ! {w}")
            print()
        if blocking:
            print("BLOCKING:")
            for b in blocking:
                print(f"  X {b}")
            print()
        if summary["ok"]:
            print("Status: READY")
        else:
            print("Status: BLOCKED — fix the items above before running the pipeline")

    sys.exit(0 if summary["ok"] else 1)


def _import_neoswga_version() -> str:
    try:
        import neoswga

        return getattr(neoswga, "__version__", "unknown")
    except Exception:
        return "unknown"


def run_init(args):
    """Run the setup wizard to create params.json"""
    import inspect

    from neoswga.core.wizard import run_wizard

    try:
        # Build kwargs, only pass blacklist_paths if wizard supports it
        wizard_kwargs = dict(
            genome_path=args.genome,
            background_path=args.background,
            output_path=args.output,
            output_dir=args.output_dir,
            interactive=not args.non_interactive,
            advanced=args.advanced,
            auto_approve=args.yes,
        )
        if getattr(args, "blacklist", None):
            sig = inspect.signature(run_wizard)
            if "blacklist_paths" in sig.parameters:
                wizard_kwargs["blacklist_paths"] = args.blacklist
            else:
                logger.warning("Blacklist support for init requires wizard update")
        run_wizard(**wizard_kwargs)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("\nSetup cancelled")
        sys.exit(1)


def run_validate_params(args):
    """Validate params.json configuration.

    Invoked as `neoswga validate params -j params.json` (preferred) or
    `neoswga validate-params -j params.json` (deprecated alias).
    """
    from neoswga.core.param_validator import validate_params_file

    if getattr(args, "command", None) == "validate-params":
        logger.warning(
            "`neoswga validate-params` is deprecated; " "use `neoswga validate params` instead."
        )

    success, _ = validate_params_file(args.json_file, verbose=True)
    sys.exit(0 if success else 1)


def run_schema(args):
    """Dump the canonical params.json JSON Schema.

    Either prints to stdout (default with --dump) or writes to --output path.
    Useful for IDE integration (VS Code can autocomplete against the schema)
    and for human-readable review.
    """
    import json

    try:
        from neoswga.core.schema import load_schema

        schema = load_schema()
    except Exception as e:
        logger.error(f"Could not load params schema: {e}")
        sys.exit(1)

    output = json.dumps(schema, indent=2)
    if getattr(args, "output", None):
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(output + "\n")
        logger.info(f"Schema written to {args.output}")
    else:
        print(output)


def run_validate_model(args):
    """Validate mechanistic model against expected behavior.

    Invoked as `neoswga validate model` (preferred) or
    `neoswga validate-model` (deprecated alias).
    """
    import json

    from neoswga.core.model_validation import (
        format_validation_report,
        validate_mechanistic_model,
    )

    if getattr(args, "command", None) == "validate-model":
        logger.warning(
            "`neoswga validate-model` is deprecated; " "use `neoswga validate model` instead."
        )

    logger.info("Validating mechanistic model...")
    results = validate_mechanistic_model()

    if getattr(args, "output_json", False):
        # Output as JSON
        # Convert results to JSON-serializable format
        json_results = []
        for r in results:
            jr = {
                "test": r.get("test", "Unknown"),
                "passed": r.get("passed", False),
                "summary": r.get("summary", ""),
            }
            if "error" in r:
                jr["error"] = r["error"]
            json_results.append(jr)
        print(json.dumps(json_results, indent=2))
    else:
        # Output as formatted report
        report = format_validation_report(results)
        print(report)

    # Exit with appropriate code
    all_passed = all(r.get("passed", False) for r in results)
    sys.exit(0 if all_passed else 1)


def run_interpret(args):
    """Interpret pipeline results"""
    from neoswga.core.results_interpreter import interpret_results

    try:
        interpret_results(args.dir, verbose=True)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)
