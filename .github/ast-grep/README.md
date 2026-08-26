# ast-grep rules

This directory contains [ast-grep](https://ast-grep.github.io/) rules that run in CI
and optionally as pre-commit hook using [prek](https://prek.j178.dev/).
`ruff check` does not support custom rules so we use ast-grep to enforce custom rules.

Current rules:

* [no-forbidden-monkeypatch](rules/no-forbidden-monkeypatch.yml): Disallow the use of `monkeypatch.setattr` in tests.

Run the test suite with:

```shell
uvx --from ast-grep-cli ast-grep test --config .github/ast-grep/sgconfig.yml
```
