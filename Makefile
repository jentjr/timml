.PHONY: docs

black:
	black timml

docs:
	cd docs && make html
	@echo "\033[95m\n\nBuild successful! View the docs homepage at docs/builddocs/html/index.html.\n\033[0m"
