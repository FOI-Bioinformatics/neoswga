.PHONY: help install install-dev install-gpu install-all test clean build upload docs format lint type-check

# Default target
.DEFAULT_GOAL := help

# Python executable
PYTHON := python3
PIP := $(PYTHON) -m pip

# Package name
PACKAGE := neoswga

help: ## Show this help message
	@echo "NeoSWGA Development Makefile"
	@echo ""
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install package with core dependencies
	$(PIP) install -e .
	@echo ""
	@echo "✓ NeoSWGA installed successfully!"
	@echo "Run: neoswga --help"

install-dev: ## Install package with development dependencies
	$(PIP) install -e ".[dev]"
	@echo ""
	@echo "✓ NeoSWGA installed with dev dependencies!"

install-gpu: ## Install package with GPU support (CUDA 11.x)
	$(PIP) install -e ".[gpu]"
	@echo ""
	@echo "✓ NeoSWGA installed with GPU support!"

install-gpu-cuda12: ## Install package with GPU support (CUDA 12.x)
	$(PIP) install -e ".[gpu-cuda12]"
	@echo ""
	@echo "✓ NeoSWGA installed with GPU support (CUDA 12)!"

install-torch: ## Install package with PyTorch deep learning
	$(PIP) install -e ".[deep-learning-torch]"
	@echo ""
	@echo "✓ NeoSWGA installed with PyTorch!"

install-tf: ## Install package with TensorFlow deep learning
	$(PIP) install -e ".[deep-learning-tf]"
	@echo ""
	@echo "✓ NeoSWGA installed with TensorFlow!"

install-all: ## Install package with all optional dependencies
	$(PIP) install -e ".[all]"
	@echo ""
	@echo "✓ NeoSWGA installed with all dependencies!"

test: ## Run test suite
	$(PYTHON) -m pytest tests/ -v --cov=$(PACKAGE) --cov-report=html --cov-report=term-missing
	@echo ""
	@echo "✓ Tests complete! Coverage report: htmlcov/index.html"

test-quick: ## Run fast tests only
	$(PYTHON) -m pytest tests/ -v -k "not slow"

test-watch: ## Run tests in watch mode
	$(PYTHON) -m pytest tests/ -v --looponfail

clean: ## Remove build artifacts and cache files
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf htmlcov/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	@echo "✓ Cleaned build artifacts!"

build: clean ## Build distribution packages
	$(PYTHON) -m build
	@echo ""
	@echo "✓ Distribution packages built in dist/"

upload: build ## Upload package to PyPI
	$(PYTHON) -m twine upload dist/*
	@echo ""
	@echo "✓ Package uploaded to PyPI!"

upload-test: build ## Upload package to TestPyPI
	$(PYTHON) -m twine upload --repository testpypi dist/*
	@echo ""
	@echo "✓ Package uploaded to TestPyPI!"

format: ## Format code with black and isort
	$(PYTHON) -m black neoswga/ tests/ --line-length 100
	$(PYTHON) -m isort neoswga/ tests/ --profile black
	@echo "✓ Code formatted!"

lint: ## Run linter (flake8)
	$(PYTHON) -m flake8 neoswga/ tests/ --max-line-length=100 --ignore=E203,W503
	@echo "✓ Linting complete!"

type-check: ## Run type checker (mypy)
	$(PYTHON) -m mypy neoswga/ --ignore-missing-imports
	@echo "✓ Type checking complete!"

check: format lint type-check test ## Run all checks (format, lint, type-check, test)
	@echo ""
	@echo "✓ All checks passed!"

docs: ## Generate documentation (requires sphinx and docs/conf.py)
	@if [ ! -f docs/conf.py ]; then echo "Error: docs/conf.py not found. Sphinx is not yet configured."; exit 1; fi
	$(PYTHON) -m sphinx -b html docs/ docs/_build/html
	@echo ""
	@echo "Documentation built in docs/_build/html/index.html"

serve-docs: docs ## Build and serve documentation locally
	cd docs/_build/html && $(PYTHON) -m http.server 8000

version: ## Show package version
	@$(PYTHON) -c "import neoswga; print(f'NeoSWGA version {neoswga.__version__}')"

info: ## Show package information
	@$(PYTHON) -c "import neoswga; neoswga.print_info()"

example: ## Run plasmid example pipeline
	@echo "Running plasmid example (see examples/plasmid_example/)..."
	@echo "  cd examples/plasmid_example"
	@echo "  neoswga count-kmers -j params.json"
	@echo "  neoswga filter -j params.json"
	@echo "  neoswga score -j params.json"
	@echo "  neoswga optimize -j params.json"

benchmark: ## Run performance benchmarks
	@echo "Running benchmarks..."
	$(PYTHON) -m pytest tests/benchmark_*.py -v --benchmark-only

profile: ## Profile code performance
	@echo "Profiling code..."
	$(PYTHON) -m cProfile -o profile.stats -m pytest tests/
	$(PYTHON) -c "import pstats; p = pstats.Stats('profile.stats'); p.sort_stats('cumulative'); p.print_stats(20)"

requirements: ## Update requirements.txt from pyproject.toml
	$(PIP) install pip-tools
	pip-compile pyproject.toml -o requirements-dev.txt --extra dev
	@echo "✓ Requirements updated!"

freeze: ## Freeze current environment
	$(PIP) freeze > requirements-frozen.txt
	@echo "✓ Environment frozen to requirements-frozen.txt"

check-deps: ## Check for outdated dependencies
	$(PIP) list --outdated
	@echo "Run 'pip install --upgrade <package>' to update"

docker-build: ## Build Docker container
	docker build -t neoswga:latest .
	@echo "✓ Docker container built!"

docker-run: ## Run Docker container
	docker run -it --rm -v $(PWD)/data:/data neoswga:latest bash

init-dev: install-dev ## Initialize development environment
	pre-commit install
	@echo ""
	@echo "✓ Development environment initialized!"
	@echo "Git hooks installed. Code will be checked before commits."

uninstall: ## Uninstall package
	$(PIP) uninstall -y $(PACKAGE)
	@echo "✓ NeoSWGA uninstalled!"

reinstall: uninstall install ## Reinstall package
	@echo "✓ NeoSWGA reinstalled!"

# Quick shortcuts
i: install ## Shortcut for install
t: test ## Shortcut for test
c: clean ## Shortcut for clean
f: format ## Shortcut for format
l: lint ## Shortcut for lint
b: build ## Shortcut for build

.PHONY: all
all: check build docs ## Run all checks and build
	@echo ""
	@echo "========================================="
	@echo "All tasks complete! Package ready."
	@echo "========================================="
