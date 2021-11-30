.PHONY: install init init-dev all shell test test-all test-fast test-slow test-time coverage mutate \
		clean clean-mutate clean-coverage clean-test

all: test-all

install:
	sudo apt install python3.8-dev ppl-dev python3-tk python3-pip
	pip3 install pipenv

install-mac:
	sudo brew install python3.8-dev ppl-dev python3-tk python3-pip
	pip3 install pipenv

init-dev-mac: install-mac
	pipenv install -d

init-dev: install
	pipenv install -d

init: install
	pipenv install

init-mac: install-mac
	pipenv install

shell:
	pipenv shell

test: test-all

test-all:
	pytest tests/

test-fast:
	pytest -m "not slow" tests/

test-slow:
	pytest -m "slow" tests/

test-time:
	@pytest --durations=0 tests/ > .pytest-durations

coverage:
	pytest --cov=pysymbrobustness --cov-branch --cov-report \
		term-missing:skip-covered --cov-report html:cov_html tests/

mutate:
	@-mutmut --paths-to-mutate=pysymbrobustness run || echo ""
	@mutmut results > .mutate-results

clean: clean-mutate clean-test clean-coverage

clean-mutate:
	@rm -f .mutate-results .mutmut-cache

clean-test:
	@rm -rf .pytest_cache .pytest-durations

clean-coverage:
	@rm -rf htmlcov .coverage