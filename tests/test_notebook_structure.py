"""Prevent accidentally saving serialized notebook JSON inside a single cell."""

import ast
import json
from pathlib import Path
import unittest


class NotebookStructureTests(unittest.TestCase):
    def test_batch_notebook_has_real_settings_and_processing_cells(self):
        path = Path(__file__).resolve().parents[1] / "batch_fastq_reader.ipynb"
        notebook = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(notebook["nbformat"], 4)
        code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
        self.assertGreaterEqual(len(code_cells), 2, "Notebook JSON must not be wrapped inside a raw/markdown cell")
        sources = ["".join(cell["source"]) for cell in code_cells]
        for source in sources:
            ast.parse(source)
        self.assertTrue(any("FASTQ_FILE =" in source for source in sources))
        self.assertTrue(any("fastq_file=FASTQ_FILE" in source for source in sources))
        for cell in notebook["cells"]:
            source = "".join(cell["source"]).strip()
            if source.startswith("{"):
                try:
                    value = json.loads(source)
                except ValueError:
                    continue
                self.assertFalse(isinstance(value, dict) and "cells" in value and "nbformat" in value,
                                 "A complete notebook has been embedded as text instead of real cells")


if __name__ == "__main__":
    unittest.main()
