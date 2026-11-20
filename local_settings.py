"""Machine-local paths for the FASTQ Reader scripts and notebooks."""

import json
import os
from pathlib import Path
import shutil

CONFIG_PATH = Path(__file__).resolve().with_name(".fastq-reader.local.json")


def read_settings() -> dict:
    try:
        return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}


def settings_path(value: str) -> Path:
    return CONFIG_PATH.parent / Path(value).expanduser()


def local_path(section: str, name: str) -> Path | None:
    value = read_settings().get(section, {}).get(name)
    return settings_path(value) if value else None


def find_tool(name: str) -> str | None:
    """Find a native executable on PATH or in explicitly configured folders."""
    search_dirs = os.get_exec_path() + [
        str(settings_path(directory)) for directory in read_settings().get("tool_dirs", [])
    ]
    for directory in search_dirs:
        folder = Path(directory)
        if not folder.is_absolute():
            continue
        executable = folder / (f"{name}.exe" if os.name == "nt" else name)
        match = shutil.which(str(executable))
        if match:
            return match
    return None
