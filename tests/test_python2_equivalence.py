"""Compare both implementations on public FlyBase and Pasilla records.

Set JUNCBASE_PY2_ROOT and JUNCBASE_PY2_PYTHON to an upstream checkout and its
Python 2 interpreter. The Python 2 checkout requires its own dependencies.
"""

import gzip
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
PY2_ROOT = os.environ.get("JUNCBASE_PY2_ROOT")
PY2_PYTHON = os.environ.get("JUNCBASE_PY2_PYTHON")


@unittest.skipUnless(PY2_ROOT and PY2_PYTHON, "Set the Python 2 reference paths")
class Python2EquivalenceTest(unittest.TestCase):
    def test_real_rox2_inputs(self) -> None:
        assert PY2_ROOT is not None and PY2_PYTHON is not None
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            # The upstream GTF parser expects transcript_id second and the BAM
            # reader prefixes contigs with chr. Both versions receive the same
            # transformed annotation, with the original coordinates unchanged.
            annotation = work / "annotation.gtf"
            lines = []
            annotation_lines = (DATA / "rox2-exons.gtf").read_text().splitlines()
            annotation_lines += (DATA / "dre4-exons.gtf").read_text().splitlines()
            for line in annotation_lines:
                fields = line.split("\t")
                attributes = [part.strip() for part in fields[8].split(";") if part.strip()]
                transcript = next(part for part in attributes if part.split(maxsplit=1)[0] == "transcript_id")
                attributes.remove(transcript)
                attributes.insert(1, transcript)
                fields[0] = "chr" + fields[0]
                fields[8] = "; ".join(attributes) + ";"
                lines.append("\t".join(fields))
            annotation.write_text("\n".join(lines) + "\n")

            outputs = []
            for label, root, python in (
                ("py2", Path(PY2_ROOT), Path(PY2_PYTHON)),
                ("py3", ROOT, Path(sys.executable)),
            ):
                env = os.environ.copy()
                env["PATH"] = str(python.parent) + os.pathsep + env["PATH"]
                env["R_HOME"] = str(python.parent.parent / "lib" / "R")
                env["PYTHONDONTWRITEBYTECODE"] = "1"

                def run(script: str, *args: str) -> subprocess.CompletedProcess[str]:
                    return subprocess.run([str(python), str(root / script), *args], cwd=root, env=env,
                                          check=True, capture_output=True, text=True)

                database = "annotation_" + label
                db_args = ("-g", str(annotation), "-d", database, "--sqlite_db_dir", str(work))
                run("build_DB_FromGTF.py", *db_args, "--initialize")
                run("build_DB_FromGTF.py", *db_args)
                with sqlite3.connect(work / database) as connection:
                    annotation_rows = tuple(
                        connection.execute(
                            "SELECT transcript_id, gene_name, chr, start, end, strand, classification "
                            f"FROM {table} ORDER BY 1, 2, 3, 4, 5, 6, 7"
                        ).fetchall() for table in ("exon", "intron")
                    )

                for run_id in ("SRR031716", "SRR031708"):
                    run("preProcess_getASEventReadCounts.py", "-s", str(DATA / f"{run_id}-loci.bam"),
                        "-o", str(work), "-n", f"{label}-{run_id}", "-l", "37")

                first = work / f"{label}-SRR031716"
                second = work / f"{label}-SRR031708"
                for sample in (first, second):
                    run("intron_exon_jctn_counts.py", "-b", str(first) + "_junctions.bed",
                        "-a", str(sample) + "_genome_reads.txt.gz", "--out_dir", str(work),
                        "--out_prefix", sample.name, "-l", "37", "-f", "31", "-o", "3")
                prefix = work / (label + "-events")
                run("getASEventReadCounts.py", "--jcn1", str(first) + "_junctions.bed",
                    "--jcn2", str(second) + "_junctions.bed",
                    "--genome_reads1", str(first) + "_genome_reads.txt.gz",
                    "--genome_reads2", str(second) + "_genome_reads.txt.gz",
                    "--ie1", str(first) + "_intron_exon_junction_counts.txt",
                    "--ie2", str(second) + "_intron_exon_junction_counts.txt",
                    "--sqlite_db_dir", str(work), "--txt_db1", database, "--txt_db2", database,
                    "--jcn_seq_len", "62", "--keep_intermediate", "-p", str(prefix))
                if label == "py3":
                    missing_reads = str(work / "missing-genome-reads.txt.gz")
                    with self.assertRaises(subprocess.CalledProcessError) as failure:
                        run("getASEventReadCounts.py", "--jcn1", str(first) + "_junctions.bed",
                            "--jcn2", str(second) + "_junctions.bed", "--genome_reads1", missing_reads,
                            "--genome_reads2", str(second) + "_genome_reads.txt.gz",
                            "--ie1", os.devnull, "--ie2", os.devnull, "--sqlite_db_dir", str(work),
                            "--txt_db1", database, "--txt_db2", database, "--jcn_seq_len", "62",
                            "--keep_intermediate", "-p", str(work / "expected-failure"))
                    self.assertIn("FileNotFoundError", failure.exception.stderr)
                    self.assertFalse((work / "expected-failure_finished.txt").exists())
                junctions = [sorted(Path(str(sample) + "_junctions.bed").read_text().splitlines()[1:])
                             for sample in (first, second)]
                genome_reads = []
                for sample in (first, second):
                    with gzip.open(str(sample) + "_genome_reads.txt.gz", "rt") as handle:
                        genome_reads.append(handle.readlines())
                events = sorted(Path(str(prefix) + "_all_AS_event_info.txt").read_text().splitlines())
                statistics = sorted(run("fisherTest_exclusionInclusion.py", "-f",
                                        str(prefix) + "_all_AS_event_info.txt", "--method", "BH")
                                    .stdout.splitlines())
                sample_root = work / (label + "-sample-tables")
                for sample in ("sample1", "sample2"):
                    sample_dir = sample_root / sample
                    sample_dir.mkdir(parents=True)
                    shutil.copyfile(str(prefix) + "_all_AS_event_info.txt",
                                    sample_dir / (sample + "_all_AS_event_info.txt"))
                tables_prefix = work / (label + "-tables")
                run("createAS_CountTables.py", "-d", str(sample_root), "-s", "sample1,sample2",
                    "-o", str(tables_prefix), "--jcn_seq_len", "62")
                table = sorted(Path(str(tables_prefix) + "_AS_exclusion_inclusion_counts.txt")
                               .read_text().splitlines())
                if label == "py3":
                    clustered = work / "clustered.txt"
                    run("clusterASExons2.py", "-d", str(sample_root), "-s", "sample1,sample2",
                        "-o", str(clustered), "--left_intron", str(work / "left.txt"),
                        "--right_intron", str(work / "right.txt"))
                    self.assertEqual(table, sorted(clustered.read_text().splitlines()))
                    normalized = work / "clustered-normalized.txt"
                    run("clusterASExons2.py", "-d", str(sample_root), "-s", "sample1,sample2",
                        "-o", str(normalized), "--left_intron", str(work / "left-normalized.txt"),
                        "--right_intron", str(work / "right-normalized.txt"),
                        "--lengthNorm", "--jcn_seq_len", "62")
                    expected_normalized = Path(str(tables_prefix) + "_AS_exclusion_inclusion_counts_lenNorm.txt")
                    self.assertEqual(sorted(expected_normalized.read_text().splitlines()),
                                     sorted(normalized.read_text().splitlines()))
                boundary_counts = [sorted(Path(str(sample) + "_intron_exon_junction_counts.txt")
                                          .read_text().splitlines()) for sample in (first, second)]
                outputs.append((annotation_rows, junctions, genome_reads, events, statistics, table,
                                boundary_counts))

            self.assertTrue(outputs[0][3], "The real-data event comparison must not be empty")
            event_types = {line.split("\t")[1] for line in outputs[0][3]}
            self.assertTrue({"cassette", "intron_retention"}.issubset(event_types))
            self.assertGreater(len(outputs[0][5]), 1, "The aggregated count table must contain events")
            self.assertTrue(outputs[0][6][0], "The boundary-count comparison must not be empty")
            self.assertEqual(outputs[0], outputs[1])


if __name__ == "__main__":
    unittest.main()
