# Development and Testing

## Local Installation

First install the Python packaging and dependency management tool [`poetry`](https://python-poetry.org/docs/#installation).
You must have Python 3.12 or greater available on your system path, which could be managed by [`pyenv`](https://github.com/pyenv/pyenv) or another package manager. 
Finally, install the dependencies of the project with:

```console
poetry install
```

## Running Tests

To test the codebase, run the following command:

```console
poetry run pytest
```

The command will:

- Execute unit tests with [`pytest`](https://docs.pytest.org/)
- Test the language typing with [`mypy`](https://mypy-lang.org/)
- Test for linting and styling errors with [`ruff`](https://docs.astral.sh/ruff/)
- Emit a testing coverage report with [`coverage`](https://coverage.readthedocs.io/)

To format the code in the library, run the following commands:

```console
poetry run ruff check --select I --fix
poetry run ruff format bedspec tests
```

To generate a code coverage report locally, run the following command:

```console
poetry run coverage html
```
