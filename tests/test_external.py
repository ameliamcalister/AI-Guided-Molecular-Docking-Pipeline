"""ChEMBL and LLM calls, with the network and the Anthropic client replaced by fakes."""
import json
from types import SimpleNamespace

import pytest

import docking_pipeline as dp


class FakeResponse:
    def __init__(self, payload):
        self._body = json.dumps(payload).encode()

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def fake_urlopen(payload, seen=None):
    def _open(req, timeout=None):
        if seen is not None:
            seen.append(req)
        return FakeResponse(payload)
    return _open


def failing_urlopen(*a, **k):
    raise OSError("network down")


# ── ChEMBL ───────────────────────────────────────────────────────────────────
def test_fetch_chembl_inhibitors_parses_and_skips_incomplete(monkeypatch):
    seen = []
    payload = {"activities": [
        {"canonical_smiles": "CCO", "pchembl_value": "9.0", "standard_type": "IC50", "molecule_pref_name": "ETHANOL"},
        {"canonical_smiles": "", "pchembl_value": "8.0", "standard_type": "Ki"},
        {"canonical_smiles": "CCN", "pchembl_value": None, "standard_type": "Ki"},
        {"canonical_smiles": "CCC", "pchembl_value": "7.0", "standard_type": "Kd", "molecule_chembl_id": "CHEMBL1"},
    ]}
    monkeypatch.setattr(dp.urllib.request, "urlopen", fake_urlopen(payload, seen))
    out = dp.fetch_chembl_inhibitors("CHEMBL4523", max_compounds=5)
    assert [c["smiles"] for c in out] == ["CCO", "CCC"]
    assert out[0]["ic50_nm"] == 1.0 and out[1]["ic50_nm"] == 100.0
    assert out[1]["name"] == "CHEMBL1"
    assert "target_chembl_id=CHEMBL4523" in seen[0]
    assert "limit=5" in seen[0]


def test_fetch_chembl_inhibitors_returns_empty_on_network_error(monkeypatch):
    monkeypatch.setattr(dp.urllib.request, "urlopen", failing_urlopen)
    assert dp.fetch_chembl_inhibitors("CHEMBL4523") == []


def test_search_chembl_target_returns_first_hit(monkeypatch):
    payload = {"targets": [{"pref_name": "EGFR", "target_chembl_id": "CHEMBL203"}]}
    monkeypatch.setattr(dp.urllib.request, "urlopen", fake_urlopen(payload))
    assert dp.search_chembl_target("EGFR") == "CHEMBL203"


def test_search_chembl_target_no_hits_or_error_returns_none(monkeypatch):
    monkeypatch.setattr(dp.urllib.request, "urlopen", fake_urlopen({"targets": []}))
    assert dp.search_chembl_target("zzz") is None
    monkeypatch.setattr(dp.urllib.request, "urlopen", failing_urlopen)
    assert dp.search_chembl_target("EGFR") is None


def test_format_chembl_context():
    assert dp.format_chembl_context([]) == ""
    comps = [{"name": "A", "smiles": "CCO", "activity_type": "IC50", "pchembl": 9.0, "ic50_nm": 1.0},
             {"name": "", "smiles": "CCN", "activity_type": "Ki", "pchembl": 8.0, "ic50_nm": 10.0}]
    text = dp.format_chembl_context(comps, max_show=1)
    assert "A" in text and "CCO" in text and "IC50: 1.0 nM" in text
    assert "CCN" not in text


# ── LLM routing ──────────────────────────────────────────────────────────────
def test_call_llm_routes_backend_and_default_models(monkeypatch):
    calls = []
    monkeypatch.setattr(dp, "_call_claude", lambda p, s, m: calls.append(("claude", m)) or "c")
    monkeypatch.setattr(dp, "_call_ollama", lambda p, s, m: calls.append(("ollama", m)) or "o")
    assert dp.call_llm("p", "s") == "c"
    assert dp.call_llm("p", "s", backend="ollama") == "o"
    assert dp.call_llm("p", "s", backend="claude", model="custom-model") == "c"
    assert calls == [("claude", dp.DEFAULT_CLOUD), ("ollama", dp.DEFAULT_LOCAL), ("claude", "custom-model")]


def test_call_claude_sends_system_and_user_prompt(monkeypatch):
    captured = {}

    class FakeMessages:
        def create(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace(content=[SimpleNamespace(text="reply")])

    fake_module = SimpleNamespace(Anthropic=lambda: SimpleNamespace(messages=FakeMessages()))
    monkeypatch.setattr(dp, "ANTHROPIC_AVAILABLE", True)
    monkeypatch.setattr(dp, "anthropic", fake_module, raising=False)
    assert dp._call_claude("user text", "system text", "some-model") == "reply"
    assert captured["model"] == "some-model" and captured["system"] == "system text"
    assert captured["messages"] == [{"role": "user", "content": "user text"}]


def test_call_claude_without_sdk_raises_clear_error(monkeypatch):
    monkeypatch.setattr(dp, "ANTHROPIC_AVAILABLE", False)
    with pytest.raises(ImportError, match="pip install anthropic"):
        dp._call_claude("p", "s", "m")


def test_call_ollama_builds_payload_and_returns_content(monkeypatch):
    seen = []
    monkeypatch.setattr(dp.urllib.request, "urlopen",
                        fake_urlopen({"message": {"content": "local reply"}}, seen))
    assert dp._call_ollama("user", "system", "qwen2.5:1.5b") == "local reply"
    body = json.loads(seen[0].data)
    assert body["model"] == "qwen2.5:1.5b" and body["stream"] is False
    assert [m["role"] for m in body["messages"]] == ["system", "user"]
