import json
import tempfile
from pathlib import Path
import importlib.util


def load_pipeline_module():
    module_path = Path(__file__).resolve().parent.parent / "unified_pipeline_clean" / "nif_downstream_code" / "pangenome_pipeline_consolidated2.py"
    spec = importlib.util.spec_from_file_location("pipeline_module", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_filter_keeps_only_gcf_complete():
    mod = load_pipeline_module()
    with tempfile.TemporaryDirectory() as tmpdir:
        tdir = Path(tmpdir)
        table = tdir / "assembly_quality.tsv"
        table.write_text(
            "assembly_accession\tassembly_level\n"
            "GCF_000001\tComplete Genome\n"
            "GCF_000002\tChromosome\n"
            "GCA_000003\tComplete Genome\n"
        )
        kept = mod.load_and_filter_assemblies(table, tdir)
        assert kept == ["GCF_000001"]
        filtered = tdir / "assembly_quality.filtered.tsv"
        assert filtered.exists()
        filtered_acc = tdir / "filtered_accessions.txt"
        assert filtered_acc.read_text().strip().split("\n") == ["GCF_000001"]
        summary = json.loads((tdir / "filter_summary.json").read_text())
        assert summary["kept"] == 1


def test_filter_empty_rejected():
    mod = load_pipeline_module()
    mod.config.require_refseq_gcf = True
    mod.config.require_complete_genome = True
    with tempfile.TemporaryDirectory() as tmpdir:
        tdir = Path(tmpdir)
        table = tdir / "assembly_quality.tsv"
        table.write_text(
            "assembly_accession\tassembly_level\n"
            "GCA_000003\tScaffold\n"
        )
        try:
            mod.load_and_filter_assemblies(table, tdir)
        except RuntimeError as exc:
            assert "No assemblies passed" in str(exc)
        else:
            raise AssertionError("Expected RuntimeError when nothing passes filter")
