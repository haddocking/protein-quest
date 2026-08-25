#!/usr/bin/env python3
"""Pre-commit/prek hook: reject monkeypatch usage except for cwd/env changes.

When executed as a script, it lints the given paths and exits non-zero on
violations. When loaded by pytest, the ``test_*`` functions below are collected
and run as unit tests.
"""

import ast
import sys
from pathlib import Path

# Methods of pytest's ``monkeypatch`` fixture that we allow in the codebase.
ALLOWED_MONKEYPATCH_METHODS: frozenset[str] = frozenset(
    {
        "chdir",
        "delenv",
        "setenv",
        "unsetenv",
        "undo",
    }
)


def _is_forbidden_monkeypatch_call(node: ast.AST) -> bool:
    """Return True for calls such as ``monkeypatch.setattr(...)``."""
    return (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "monkeypatch"
        and node.func.attr not in ALLOWED_MONKEYPATCH_METHODS
    )


def find_violations(source: str, filename: str = "<unknown>") -> list[str]:
    """Return human-readable violation messages for a single source string."""
    try:
        tree = ast.parse(source, filename=filename)
    except SyntaxError as exc:  # pragma: no cover - covered by test_syntax_error
        return [f"{filename}: syntax error: {exc}"]

    return [
        f"{filename}:{node.lineno}: forbidden monkeypatch.{node.func.attr}()"  # pyrefly: ignore[missing-attribute]
        for node in ast.walk(tree)
        if _is_forbidden_monkeypatch_call(node)
    ]


def check_file(path: Path) -> list[str]:
    """Check a single Python file."""
    return find_violations(path.read_text(encoding="utf-8"), filename=str(path))


def check_paths(paths: list[Path]) -> list[str]:
    """Check one or more paths (files or directories)."""
    violations: list[str] = []
    for path in paths:
        if path.is_dir():
            for py_file in sorted(path.rglob("*.py")):
                violations.extend(check_file(py_file))
        elif path.suffix == ".py":
            violations.extend(check_file(path))
    return violations


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for use as a pre-commit/prek hook."""
    args = argv if argv is not None else sys.argv[1:]
    if not args:
        # Default to the test suite; src/ is not monkeypatch checked.
        args = ["tests/"]

    violations = check_paths([Path(arg) for arg in args])
    if violations:
        sys.stderr.write("Forbidden monkeypatch usage detected:\n")
        for violation in violations:
            sys.stderr.write(f"{violation}\n")
        return 1
    return 0


# ---------------------------------------------------------------------------
# Internal pytest tests. These are collected when pytest loads this file as a
# test module; they do not run when the file is executed as a hook script.
# ---------------------------------------------------------------------------


def test_allowed_chdir_is_not_flagged() -> None:
    assert find_violations("monkeypatch.chdir('/tmp')") == []


def test_allowed_env_methods_are_not_flagged() -> None:
    assert find_violations("monkeypatch.setenv('X', 'y')") == []
    assert find_violations("monkeypatch.delenv('X')") == []
    assert find_violations("monkeypatch.unsetenv('X')") == []
    assert find_violations("monkeypatch.undo()") == []


def test_forbidden_setattr_is_flagged() -> None:
    violations = find_violations("monkeypatch.setattr('requests.get', mock)")
    assert len(violations) == 1
    assert "forbidden monkeypatch.setattr" in violations[0]


def test_forbidden_delattr_is_flagged() -> None:
    violations = find_violations("monkeypatch.delattr(obj, 'name')")
    assert len(violations) == 1
    assert "forbidden monkeypatch.delattr" in violations[0]


def test_forbidden_setitem_is_flagged() -> None:
    violations = find_violations("monkeypatch.setitem(d, 'k', 'v')")
    assert len(violations) == 1
    assert "forbidden monkeypatch.setitem" in violations[0]


def test_forbidden_delitem_is_flagged() -> None:
    violations = find_violations("monkeypatch.delitem(d, 'k')")
    assert len(violations) == 1
    assert "forbidden monkeypatch.delitem" in violations[0]


def test_forbidden_syspath_prepend_is_flagged() -> None:
    violations = find_violations("monkeypatch.syspath_prepend('/tmp')")
    assert len(violations) == 1
    assert "forbidden monkeypatch.syspath_prepend" in violations[0]


def test_non_monkeypatch_calls_are_ignored() -> None:
    assert find_violations("other.setattr('x', 'y')") == []
    assert find_violations("mp.setattr('x', 'y')") == []


def test_multiline_source_reports_correct_line() -> None:
    source = """
def test_foo(monkeypatch):
    monkeypatch.chdir('/tmp')
    monkeypatch.setenv('A', 'b')
    monkeypatch.setattr('requests.get', mock)
"""
    violations = find_violations(source)
    assert len(violations) == 1
    assert ":5: forbidden monkeypatch.setattr" in violations[0]


def test_syntax_error_is_reported() -> None:
    violations = find_violations("def foo(", filename="bad.py")
    assert len(violations) == 1
    assert "bad.py: syntax error" in violations[0]


if __name__ == "__main__":
    raise SystemExit(main())
