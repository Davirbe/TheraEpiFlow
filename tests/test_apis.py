"""
API and tool availability tests for TheraEPIflow.

Tests each external dependency and prints the output format
so we can validate the pipeline logic before building modules.

Run with:
    python tests/test_apis.py
"""

import sys
import json
import requests

# ── Test sequences ────────────────────────────────────────────────────────────
TEST_PEPTIDES  = ["SLYNTVATLY", "GILGFVFTL", "NLVPMVATV"]
TEST_ALLELE    = "HLA-A*02:01"
TEST_ACCESSION = "NP_085468.1"   # ZIKV envelope protein


# ── Helpers ───────────────────────────────────────────────────────────────────

def section(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")

def ok(msg):    print(f"  ✅ {msg}")
def fail(msg):  print(f"  ✗  {msg}")
def info(msg):  print(f"  ℹ  {msg}")


# ── 1. Biopython Entrez ───────────────────────────────────────────────────────

def test_entrez():
    section("1. Biopython Entrez — GenBank fetch")
    try:
        from Bio import Entrez, SeqIO
        Entrez.email = "test@theraEPIflow.dev"

        handle = Entrez.efetch(
            db="protein", id=TEST_ACCESSION,
            rettype="fasta", retmode="text"
        )
        record = SeqIO.read(handle, "fasta")

        ok(f"Fetched: {record.id}")
        ok(f"Description: {record.description[:60]}...")
        ok(f"Length: {len(record.seq)} aa")
        info(f"First 30 aa: {str(record.seq[:30])}")
        return True

    except ImportError:
        fail("Biopython not installed — run: conda install biopython")
        return False
    except Exception as e:
        fail(f"Entrez error: {e}")
        return False


# ── 2. IEDB MHC-I Prediction API ─────────────────────────────────────────────

def test_iedb_mhci():
    section("2. IEDB MHC-I Prediction API")
    url = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
    payload = {
        "method":        "recommended",
        "sequence_text": TEST_PEPTIDES[0],
        "allele":        TEST_ALLELE,
        "length":        "9",
    }
    try:
        response = requests.post(url, data=payload, timeout=30)
        ok(f"Status: {response.status_code}")

        if response.status_code == 200:
            lines = response.text.strip().split("\n")
            ok(f"Lines returned: {len(lines)}")
            info(f"Header: {lines[0]}")
            if len(lines) > 1:
                info(f"First result: {lines[1]}")
        else:
            fail(f"Unexpected status: {response.status_code}")
            info(f"Response: {response.text[:200]}")
        return response.status_code == 200

    except Exception as e:
        fail(f"IEDB MHC-I API error: {e}")
        return False


# ── 3. IEDB Cluster API ───────────────────────────────────────────────────────

def test_iedb_cluster():
    section("3. IEDB Cluster API")
    url = "https://tools-cluster-interface.iedb.org/tools_api/cluster/"
    peptide_text = "\n".join(TEST_PEPTIDES)
    payload = {
        "peptides": peptide_text,
        "cutoff":   "0.9",
    }
    try:
        response = requests.post(url, data=payload, timeout=30)
        ok(f"Status: {response.status_code}")

        if response.status_code == 200:
            lines = response.text.strip().split("\n")
            ok(f"Lines returned: {len(lines)}")
            info(f"Header: {lines[0]}")
            if len(lines) > 1:
                info(f"First result: {lines[1]}")
        else:
            fail(f"Unexpected status: {response.status_code}")
            info(f"Response: {response.text[:300]}")
        return response.status_code == 200

    except Exception as e:
        fail(f"IEDB Cluster API error: {e}")
        return False


# ── 4. IEDB Conservancy API ───────────────────────────────────────────────────

def test_iedb_conservancy():
    section("4. IEDB Conservancy API")
    url = "https://tools-cluster-interface.iedb.org/tools_api/conservancy/"

    # Minimal test: epitope against a mock reference sequence
    test_fasta = ">seq1\nSLYNTVATLYGILGFVFTLNLVPMVATV\n"
    payload = {
        "query_sequence":      TEST_PEPTIDES[0],
        "reference_sequences": test_fasta,
        "similarity_threshold": 1.0,
    }
    try:
        response = requests.post(url, data=payload, timeout=30)
        ok(f"Status: {response.status_code}")
        info(f"Response preview: {response.text[:300]}")
        return response.status_code == 200

    except Exception as e:
        fail(f"IEDB Conservancy API error: {e}")
        return False


# ── 5. MHCFlurry ─────────────────────────────────────────────────────────────

def test_mhcflurry():
    section("5. MHCFlurry 2.0")
    try:
        from mhcflurry import Class1PresentationPredictor
        predictor = Class1PresentationPredictor.load()
        ok("Predictor loaded")

        result = predictor.predict(
            peptides=TEST_PEPTIDES,
            alleles=[TEST_ALLELE] * len(TEST_PEPTIDES),
        )
        ok(f"Predictions returned: {len(result)} rows")
        info(f"Columns: {list(result.columns)}")
        info(f"Sample:\n{result.head(2).to_string()}")
        return True

    except ImportError:
        fail("MHCFlurry not installed — run: pip install mhcflurry && mhcflurry-downloads fetch")
        return False
    except Exception as e:
        fail(f"MHCFlurry error: {e}")
        return False


# ── 6. NetMHCpan ─────────────────────────────────────────────────────────────

def test_netmhcpan():
    section("6. NetMHCpan 4.1")
    import subprocess
    try:
        result = subprocess.run(
            ["netMHCpan", "-h"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0 or "Usage" in result.stdout + result.stderr:
            ok("netMHCpan found and responsive")
            info(f"Version info: {(result.stdout + result.stderr)[:200]}")
            return True
        else:
            fail("netMHCpan returned unexpected output")
            return False

    except FileNotFoundError:
        fail("netMHCpan not found in PATH")
        info("Download from: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/")
        info("Requires academic registration and manual installation")
        return False
    except Exception as e:
        fail(f"netMHCpan error: {e}")
        return False


# ── 7. ToxinPred2 ─────────────────────────────────────────────────────────────

def test_toxinpred2():
    section("7. ToxinPred2")
    import subprocess

    # Try as CLI tool
    for cmd in ["toxinpred2", "toxinpred2.py", "python -m toxinpred2"]:
        try:
            result = subprocess.run(
                cmd.split() + ["--help"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0 or "usage" in (result.stdout + result.stderr).lower():
                ok(f"ToxinPred2 found via: {cmd}")
                return True
        except FileNotFoundError:
            continue

    # Try as Python module
    try:
        import toxinpred2
        ok("ToxinPred2 available as Python module")
        return True
    except ImportError:
        pass

    fail("ToxinPred2 not found")
    info("Options:")
    info("  pip install toxinpred2")
    info("  Or download from: https://webs.iiitd.edu.in/raghava/toxinpred2/")
    return False


# ── 8. Population Coverage (pickle) ──────────────────────────────────────────

def test_population_coverage():
    section("8. Population Coverage Pickle")
    try:
        import sys
        from pathlib import Path
        pickle_dir = Path("existing_scripts/population_coverage_pickle")
        sys.path.insert(0, str(pickle_dir.parent))

        from population_coverage_pickle import population_coverage, country_ethnicity
        ok("Pickle loaded successfully")
        ok(f"Populations available: {len(population_coverage)}")
        sample_keys = list(population_coverage.keys())[:3]
        info(f"Sample populations: {sample_keys}")
        return True

    except Exception as e:
        fail(f"Population coverage pickle error: {e}")
        info("Check existing_scripts/population_coverage_pickle/__init__.py")
        return False


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "="*60)
    print("  TheraEPIflow — API & Tool Availability Tests")
    print("="*60)

    results = {
        "Biopython Entrez":         test_entrez(),
        "IEDB MHC-I API":           test_iedb_mhci(),
        "IEDB Cluster API":         test_iedb_cluster(),
        "IEDB Conservancy API":     test_iedb_conservancy(),
        "MHCFlurry":                test_mhcflurry(),
        "NetMHCpan":                test_netmhcpan(),
        "ToxinPred2":               test_toxinpred2(),
        "Population Coverage":      test_population_coverage(),
    }

    section("SUMMARY")
    all_pass = True
    for name, passed in results.items():
        status = "✅ OK" if passed else "✗  FAIL"
        print(f"  {status}  {name}")
        if not passed:
            all_pass = False

    print()
    if all_pass:
        print("  All tests passed — pipeline ready to build.")
    else:
        failed = [n for n, p in results.items() if not p]
        print(f"  {len(failed)} dependency/dependencies need attention before proceeding.")
    print()


if __name__ == "__main__":
    main()
