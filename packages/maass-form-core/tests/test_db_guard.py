"""Regression tests for the maass_form_core.database [db]-extras guard."""

import importlib
import sys

import pytest


def _force_reimport(monkeypatch, missing):
    monkeypatch.setitem(sys.modules, missing, None)
    monkeypatch.delitem(sys.modules, "maass_form_core.database", raising=False)
    monkeypatch.delitem(sys.modules, "maass_form_core.database.models", raising=False)
    monkeypatch.delitem(sys.modules, "maass_form_core.database.queryset", raising=False)


@pytest.mark.parametrize("missing", ["mongoengine", "comp_manager"])
def test_guard_fires_when_extras_missing(monkeypatch, missing):
    _force_reimport(monkeypatch, missing)
    with pytest.raises(ModuleNotFoundError, match=r"\[db\] extras"):
        importlib.import_module("maass_form_core.database")


@pytest.mark.parametrize("missing", ["mongoengine", "comp_manager"])
def test_guard_message_names_install_command(monkeypatch, missing):
    _force_reimport(monkeypatch, missing)
    with pytest.raises(ModuleNotFoundError, match=r"pip install 'maass_form_core\[db\]'"):
        importlib.import_module("maass_form_core.database")


@pytest.mark.parametrize("missing", ["mongoengine", "comp_manager"])
def test_guard_preserves_original_cause(monkeypatch, missing):
    _force_reimport(monkeypatch, missing)
    with pytest.raises(ModuleNotFoundError) as excinfo:
        importlib.import_module("maass_form_core.database")
    assert excinfo.value.__cause__ is not None
    assert isinstance(excinfo.value.__cause__, ImportError)
