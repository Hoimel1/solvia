# SOLVIA Makefile
# Convenient commands for development and deployment

.PHONY: help setup test clean docker-build docker-up run-pipeline

help:
	@echo "SOLVIA - Hemolytic Toxicity Prediction"
	@echo ""
	@echo "Available commands:"
	@echo "  make setup          - Set up local environment"
	@echo "  make test           - Run tests"
	@echo "  make lint           - Run code linting"
	@echo "  make format         - Format code with black"
	@echo "  make clean          - Clean intermediate files"
	@echo "  make docker-build   - Build Docker containers"
	@echo "  make docker-up      - Start Docker services"
	@echo "  make run-pipeline   - Run full pipeline"
	@echo "  make api            - Start API server"
	@echo "  make docs           - Build documentation"

setup:
	@echo "Setting up SOLVIA..."
	chmod +x setup_local.sh
	./setup_local.sh

test:
	@echo "Running tests..."
	python test_installation.py
	pytest tests/ -v --cov=src

lint:
	@echo "Running linters..."
	flake8 src --max-line-length=127 --extend-ignore=E203
	mypy src --ignore-missing-imports

format:
	@echo "Formatting code..."
	black src tests
	isort src tests

clean:
	@echo "Cleaning intermediate files..."
	rm -rf data/intermediate/*
	rm -rf data/output/*
	rm -rf logs/*
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

docker-build:
	@echo "Building Docker containers..."
	docker-compose build

docker-up:
	@echo "Starting Docker services..."
	docker-compose up -d

docker-down:
	@echo "Stopping Docker services..."
	docker-compose down

run-pipeline:
	@echo "Running SOLVIA pipeline..."
	snakemake --cores all --use-conda

run-pipeline-docker:
	@echo "Running SOLVIA pipeline in Docker..."
	docker-compose run --rm solvia-base snakemake --cores all

api:
	@echo "Starting API server..."
	uvicorn src.api.main:app --reload --host 0.0.0.0 --port 8000

api-docker:
	@echo "Starting API server in Docker..."
	docker-compose up api

jupyter:
	@echo "Starting Jupyter Lab..."
	jupyter lab --ip=0.0.0.0 --port=8888 --no-browser

docs:
	@echo "Building documentation..."
	cd docs && make html
	@echo "Documentation available at: docs/_build/html/index.html"

validate:
	@echo "Validating setup..."
	snakemake --lint
	snakemake -n

benchmark:
	@echo "Running benchmarks..."
	python tests/benchmark.py

deploy:
	@echo "Deploying to production..."
	@echo "Not implemented yet"

# Development helpers
dev-install:
	pip install -e .
	pip install -r requirements-dev.txt

update-deps:
	pip-compile requirements.in
	pip-compile requirements-dev.in

