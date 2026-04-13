"""
Global default configuration for TheraEPIflow pipeline.
These values can be overridden per project in project_config.json.
"""

# ── Peptide lengths ────────────────────────────────────────────────────────────
DEFAULT_PEPTIDE_LENGTHS = [9]
MURINE_PEPTIDE_LENGTHS  = [8, 9]

# ── NetMHCpan thresholds ───────────────────────────────────────────────────────
NETMHC_STRONG_BINDER_EL_RANK = 0.5    # % — strong binder cutoff
NETMHC_WEAK_BINDER_EL_RANK   = 2.0    # % — weak binder cutoff
NETMHC_IC50_CUTOFF            = 500    # nM

# ── MHCFlurry thresholds ───────────────────────────────────────────────────────
FLURRY_AFFINITY_CUTOFF        = 500    # nM
FLURRY_PRESENTATION_MIN       = 0.5   # presentation score minimum

# ── ToxinPred2 ─────────────────────────────────────────────────────────────────
TOXICITY_SCORE_THRESHOLD      = 0.6   # above = toxic

# ── IEDB Clustering ────────────────────────────────────────────────────────────
CLUSTER_IDENTITY_CUTOFF       = 0.9   # 90% sequence identity

# ── Conservation ──────────────────────────────────────────────────────────────
CONSERVATION_HIGH_THRESHOLD     = 0.90
CONSERVATION_MODERATE_THRESHOLD = 0.70
CONSERVATION_SIMILARITY_DEFAULT = 1.0  # pairwise identity cutoff (1.0 = 100%)

# ── Population coverage ────────────────────────────────────────────────────────
COVERAGE_MINIMUM_GLOBAL       = 0.80  # 80% minimum global coverage

# ── Murine prediction ─────────────────────────────────────────────────────────
MURINE_STRONG_BINDER_EL_RANK  = 2.0   # % — strong binder cutoff for H-2
MURINE_PROMISCUOUS_MIN_ALLELES = 2    # minimum H-2 alleles for promiscuity

MURINE_ALLELES = {
    "C57BL/6": ["H-2Db", "H-2Kb"],
    "BALB/c":  ["H-2Dd", "H-2Kd", "H-2Ld"],
    "complete": ["H-2Db", "H-2Kb", "H-2Dd", "H-2Kd", "H-2Ld"],
}

# ── Variant search ─────────────────────────────────────────────────────────────
VARIANTS_MAX_RESULTS          = 100

# ── NCBI Entrez ────────────────────────────────────────────────────────────────
ENTREZ_EMAIL                  = ""    # set by user on project creation
ENTREZ_DATABASE               = "protein"
