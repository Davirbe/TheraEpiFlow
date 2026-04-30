"""
Global default configuration for TheraEPIflow pipeline.
These values can be overridden per project in project_config.json.
"""

# ── Default HLA alleles ────────────────────────────────────────────────────────
# 27 MHC-I alleles covering global population diversity.
# Subset of IEDB Calis 2013-supported alleles; compatible with NetMHCpan EL and MHCFlurry 2.0.
# Each project can override this list in project_config['hla_alleles'].
DEFAULT_HLA_ALLELES = [
    # HLA-A (15 alleles)
    'HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06',
    'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02',
    'HLA-A*26:01', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02',
    'HLA-A*31:01', 'HLA-A*33:01', 'HLA-A*68:01',
    # HLA-B (12 alleles)
    'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*18:01',
    'HLA-B*27:05', 'HLA-B*35:01', 'HLA-B*39:01', 'HLA-B*40:01',
    'HLA-B*44:02', 'HLA-B*51:01', 'HLA-B*57:01', 'HLA-B*58:01',
]

# ── Peptide lengths ────────────────────────────────────────────────────────────
DEFAULT_PEPTIDE_LENGTHS = [9]
MURINE_PEPTIDE_LENGTHS  = [8, 9]

# ── NetMHCpan thresholds ───────────────────────────────────────────────────────
NETMHC_STRONG_BINDER_EL_RANK = 0.5    # % — strong binder cutoff
NETMHC_WEAK_BINDER_EL_RANK   = 2.0    # % — weak binder cutoff
NETMHC_IC50_CUTOFF            = 500    # nM

# ── MHCFlurry thresholds ───────────────────────────────────────────────────────
FLURRY_AFFINITY_CUTOFF        = 500    # nM
FLURRY_PRESENTATION_MIN       = 0.5    # presentation score minimum

# ── Step 04 — Consensus filter (two-stage) ─────────────────────────────────────
# Stage 1: parallel intersection of presentation tools (peptide must pass both).
CONSENSUS_NETMHCPAN_EL_RANK_MAX_PERCENT       = 2.0   # NetMHCpan EL %Rank ≤ 2%
CONSENSUS_MHCFLURRY_PRESENTATION_PERCENTILE_MAX = 2.0  # MHCflurry presentation %ile ≤ 2%

# Stage 2: sequential Calis immunogenicity (applied only to Stage 1 survivors).
CONSENSUS_CALIS_IMMUNOGENICITY_SCORE_MIN      = 0.0   # > 0 means immunogenic

# IEDB Calis et al. 2013 immunogenicity tool (HTTP POST endpoint)
IEDB_IMMUNOGENICITY_API_URL = 'http://tools-cluster-interface.iedb.org/tools_api/immunogenicity/'

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
