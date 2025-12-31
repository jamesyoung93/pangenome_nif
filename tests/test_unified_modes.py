import importlib.util
from pathlib import Path


def load_pipeline_module():
    module_path = Path(__file__).resolve().parent.parent / "unified_pipeline_clean" / "nif_downstream_code" / "pangenome_pipeline_consolidated2.py"
    spec = importlib.util.spec_from_file_location("pipeline_module", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_full_mode_runs_directionality():
    mod = load_pipeline_module()
    steps = [num for num, _, _ in mod.pipeline_steps_for_mode("full")]
    assert 5 in steps


def test_experiment_mode_skips_directionality():
    mod = load_pipeline_module()
    steps = [num for num, _, _ in mod.pipeline_steps_for_mode("experiment")]
    assert 5 not in steps


def test_pipeline_alias_maps_to_full():
    mod = load_pipeline_module()
    assert mod.resolve_mode("pipeline") == "full"
