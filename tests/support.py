"""Small helpers shared by the repository's standard-library tests."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_script_module(filename: str, module_name: str):
    """Load one script directly from ``scripts/`` without installing a package."""
    path = REPO_ROOT / "scripts" / filename
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load script module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    return module


def write_text(path: Path, text: str) -> Path:
    """Write UTF-8 text, creating the requested parent directories."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def write_csv(path: Path, header: Sequence[str], rows: Sequence[Sequence[object]]) -> Path:
    """Write a small CSV file with the standard CSV writer."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def run_python(
    args: Sequence[str],
    cwd: Path | str = REPO_ROOT,
) -> subprocess.CompletedProcess[str]:
    """Run the current Python interpreter and capture text output."""
    return subprocess.run(
        [sys.executable, *map(str, args)],
        cwd=Path(cwd),
        text=True,
        capture_output=True,
        check=False,
    )


def normalize_input_line(report: str, replacement: str = "<ABSOLUTE_PATH>") -> str:
    """Replace only the path value on each report line beginning with ``Input:``."""
    normalized = []
    for line in report.splitlines(keepends=True):
        body = line.rstrip("\r\n")
        ending = line[len(body) :]
        normalized.append(
            f"Input: {replacement}{ending}" if body.startswith("Input:") else line
        )
    return "".join(normalized)
