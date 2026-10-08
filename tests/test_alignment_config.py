"""Configuration boundaries, command workflows and reproducible diagnostics."""
from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from deepstructgenomics.alignment.config import AlignmentConfig, alignment_cell_count, load_alignment_config, load_sequence_inputs
from deepstructgenomics.alignment import align_sequences, AlignmentResult
from deepstructgenomics.alignment.benchmark import generate_test_pair


def test_priority_and_no_environment_mutation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    env = tmp_path / '.env'
    env.write_text('DSG_MAX_ALIGNMENT_CELLS=101\nDSG_ALIGNMENT_MATCH="3"\n', encoding='utf-8')
    system = {'DSG_MAX_ALIGNMENT_CELLS': '102'}
    assert load_alignment_config(environ={}).max_cells == 101
    assert load_alignment_config(environ=system).max_cells == 102
    assert load_alignment_config(environ=system, cli={'max_cells': 103}).max_cells == 103
    effective = load_alignment_config(environ=system, cli={'max_cells': 103}, gui={'max_cells': 104})
    assert effective.max_cells == 104 and effective.match_score == 3
    assert effective.max_length == AlignmentConfig().max_length
    assert system == {'DSG_MAX_ALIGNMENT_CELLS': '102'}
    env.unlink()
    assert load_alignment_config(environ={}) == AlignmentConfig()


@pytest.mark.parametrize('value', ['-1', '0', 'abc', '1.5'])
def test_invalid_limits_fail_fast(tmp_path, value):
    env = tmp_path / '.env'
    env.write_text('')
    with pytest.raises(ValueError, match='max_cells|DSG_MAX_ALIGNMENT_CELLS'):
        load_alignment_config(env, environ={'DSG_MAX_ALIGNMENT_CELLS': value})


@pytest.mark.parametrize('kwargs', [{'max_length': 0}, {'max_cells': True}, {'enabled': 'false'},
                                    {'match_score': float('nan')}, {'gap_open_score': float('inf')},
                                    {'overflow_policy': 'automatic'}])
def test_invalid_model(kwargs):
    with pytest.raises(ValueError):
        AlignmentConfig(**kwargs)


def test_limit_boundary_and_explicit_fallback(monkeypatch):
    config = AlignmentConfig(max_cells=100)
    assert align_sequences('ACGUACGUAC', 'ACGUACGUAC', config)
    monkeypatch.setattr('deepstructgenomics.alignment.pairwise.PairwiseAligner', lambda **kw: pytest.fail('must reject before allocation'))
    with pytest.raises(ValueError, match='110 cellules'):
        align_sequences('A' * 10, 'A' * 11, config)
    assert align_sequences('A' * 10, 'A' * 11, replace(config, overflow_policy='positional')) is None
    assert align_sequences('A' * 10, 'A' * 11, replace(config, enabled=False)) is None
    assert alignment_cell_count(10, 11) == 110


def test_source_layers_replace_sources_per_side(tmp_path):
    env = tmp_path / '.env'
    env.write_text('DSG_WT_FASTA=old.fa\nDSG_MUT_SEQUENCE=ACGU\n')
    assert load_sequence_inputs(env, environ={'DSG_WT_SEQUENCE': 'ACGA'}, cli={'sequence': 'GGGG'}) == {'sequence': 'GGGG', 'mutant_sequence': 'ACGU'}
    with pytest.raises(ValueError, match='seule source'):
        load_sequence_inputs(env, environ={}, cli={'sequence': 'A', 'fasta_path': 'a.fa'})


@pytest.mark.parametrize('pattern', ['deterministic-random', 'repetitive'])
def test_generator_reproducible_without_global_random_state(pattern):
    import random
    state = random.getstate()
    first = generate_test_pair(100, 10, 3, 2, 42, pattern)
    assert first == generate_test_pair(100, 10, 3, 2, 42, pattern)
    assert len(first[0]) == 100 and len(first[1]) == 101
    assert random.getstate() == state
    assert first != generate_test_pair(100, 10, 3, 2, 43, pattern)


def test_persist_effective_config_and_fallback(tmp_path):
    from deepstructgenomics.gui.services import run_rna, read_rna
    config = AlignmentConfig(max_cells=100, match_score=3)
    entry = run_rna(tmp_path, {'sequence': 'ACGU', 'mutant_sequence': 'ACGGU', 'alignment_config': config.to_dict()})
    data, _ = read_rna(entry['path'])
    assert data['alignment_run']['configuration'] == config.to_dict()
    assert data['alignment_run']['matrix']['cells'] == 20
    result = AlignmentResult.from_dict(data['alignment'], 'ACGU', 'ACGGU')
    assert result.run_metadata == data['alignment_run'] and result.ambiguous
    config = replace(config, max_cells=10, overflow_policy='positional')
    entry = run_rna(tmp_path, {'sequence': 'ACGU', 'mutant_sequence': 'ACGGU', 'alignment_config': config.to_dict()})
    data, _ = read_rna(entry['path'])
    assert data['alignment'] is None
    assert data['alignment_run']['status'] == 'overflow_positional'
    assert data['alignment_run']['configuration'] == config.to_dict()
    assert '20 cellules' in data['alignment_run']['reason']


def test_cli_fasta_env_override_and_short_output(tmp_path):
    root = Path(__file__).resolve().parents[1]
    wt, mut = tmp_path / 'wt.fa', tmp_path / 'mut.fa'
    wt.write_text('>wt\nACGU\n')
    mut.write_text('>mut\nACGGU\n')
    config = tmp_path / '.env'
    config.write_text('DSG_MAX_ALIGNMENT_CELLS=19\n')
    environment = {key: value for key, value in os.environ.items() if not key.startswith('DSG_')}
    environment['PYTHONUTF8'] = '1'
    command = [sys.executable, str(root / 'scripts/run_alignment.py'), '--wt-fasta', str(wt), '--mut-fasta', str(mut)]
    rejected = subprocess.run(command, cwd=tmp_path, env=environment, capture_output=True, text=True, encoding='utf-8')
    assert rejected.returncode == 2 and 'rejected' in rejected.stdout
    passed = subprocess.run(command + ['--max-cells', '20'], cwd=tmp_path, env=environment, capture_output=True, text=True, encoding='utf-8')
    assert passed.returncode == 0 and 'aligned' in passed.stdout
    assert 'ACGGU' not in passed.stdout
    environment['DSG_MAX_ALIGNMENT_CELLS'] = '22'
    shown = subprocess.run([sys.executable, str(root / 'scripts/run_alignment.py'), '--show-config'], cwd=tmp_path, env=environment, capture_output=True, text=True)
    assert json.loads(shown.stdout)['max_cells'] == 22


def test_benchmark_cli_is_synthetic_and_does_not_fold(tmp_path, monkeypatch, capsys):
    from deepstructgenomics.alignment.cli import main
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('deepstructgenomics.rna.secondary_structure.predict_secondary_structure', lambda *a: pytest.fail('no folding'))
    report = tmp_path / 'benchmark.json'
    assert main(['--length', '10', '--max-cells', '100', '--json-output', str(report)], benchmark=True) == 0
    payload = json.loads(report.read_text())
    assert payload['synthetic'] and payload['cells'] == 100 and payload['status'] == 'aligned'
    assert payload['time_seconds'] >= 0
    assert 'ACGUACGUAC' not in capsys.readouterr().out
