"""Integration check using FlyBase r6.69 roX2 exon records."""

import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GTF = ROOT / "tests" / "data" / "rox2-exons.gtf"


class FlyBaseDatabaseTest(unittest.TestCase):
    def test_rox2_transcripts_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            command = [
                sys.executable, str(ROOT / "build_DB_FromGTF.py"),
                "-g", str(GTF), "-d", "rox2", "--sqlite_db_dir", temporary_directory,
            ]
            subprocess.run([*command, "--initialize"], check=True, capture_output=True, text=True)
            subprocess.run(command, check=True, capture_output=True, text=True)
            with sqlite3.connect(Path(temporary_directory) / "rox2") as database:
                rows = database.execute(
                    "SELECT transcript_id, gene_name, start, end FROM exon "
                    "ORDER BY transcript_id, start"
                ).fetchall()

        self.assertEqual(
            rows,
            [
                ("FBtr0073514", "FBgn0019660", 11579065, 11579124),
                ("FBtr0073514", "FBgn0019660", 11579670, 11580171),
                ("FBtr0308846", "FBgn0019660", 11579065, 11579124),
                ("FBtr0308846", "FBgn0019660", 11579386, 11579399),
                ("FBtr0308846", "FBgn0019660", 11579512, 11579528),
                ("FBtr0308846", "FBgn0019660", 11579670, 11580100),
            ],
        )


if __name__ == "__main__":
    unittest.main()
