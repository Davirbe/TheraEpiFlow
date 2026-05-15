"""
Global default configuration for TheraEPIflow pipeline.
These values can be overridden per project in project_config.json.
"""

# ── Target host ────────────────────────────────────────────────────────────────
# Fixed: the pipeline targets human HLA-I, so the host is always Homo sapiens.
# Murine prediction uses a separate H-2 allele set (see MURINE_ALLELES in this file).
TARGET_HOST = 'Homo sapiens'

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

# ── ToxinPred3 ─────────────────────────────────────────────────────────────────
TOXICITY_SCORE_THRESHOLD      = 0.38  # above = toxic (ToxinPred3 default, calibrated for peptides)

# ── IEDB Clustering ────────────────────────────────────────────────────────────
CLUSTER_IDENTITY_CUTOFF       = 0.8   # 80% sequence identity

# ── Conservation ──────────────────────────────────────────────────────────────
CONSERVATION_HIGH_THRESHOLD     = 0.90
CONSERVATION_MODERATE_THRESHOLD = 0.80
CONSERVATION_SIMILARITY_DEFAULT = 1.0  # pairwise identity cutoff (1.0 = 100%)

# ── Population coverage ────────────────────────────────────────────────────────
COVERAGE_MINIMUM_GLOBAL       = 0.80  # 80% minimum global coverage

# ── Murine prediction ─────────────────────────────────────────────────────────
# Four-tier binder labels (parallel to the human pipeline percentile thresholds):
#   optimal     percentile ≤ MURINE_OPTIMAL_BINDER_RANK_MAX        ( ≤ 0.5 )
#   good        percentile ≤ MURINE_STRONG_BINDER_EL_RANK          ( ≤ 2.0 )
#   borderline  percentile ≤ MURINE_BORDERLINE_BINDER_RANK_MAX     ( ≤ 2.5 )
#   non_binder  percentile > MURINE_BORDERLINE_BINDER_RANK_MAX     (  > 2.5 )
MURINE_OPTIMAL_BINDER_RANK_MAX     = 0.5
MURINE_STRONG_BINDER_EL_RANK       = 2.0
MURINE_BORDERLINE_BINDER_RANK_MAX  = 2.5
MURINE_PROMISCUOUS_MIN_ALLELES     = 2    # kept for backward reference; not used by predict_murine

# Five-group H-2 class I strain map. Only alleles consolidated in both
# NetMHCpan-4.1 and MHCFlurry-2.0 are listed. SJL (H-2s) and DBA/1 (H-2q)
# are deliberately omitted. Format: "H-2-Db" (canonical for both predictors).
MURINE_ALLELES = {
    "default":      ["H-2-Db", "H-2-Kb",                       # C57BL/6
                     "H-2-Dd", "H-2-Kd", "H-2-Ld"],            # BALB/c   → 5 alleles
    "C57BL/6":      ["H-2-Db", "H-2-Kb"],                      # 2
    "BALB/c":       ["H-2-Dd", "H-2-Kd", "H-2-Ld"],            # 3
    "CBA/C3H/AKR":  ["H-2-Dk", "H-2-Kk"],                      # 2  (haplotype H-2k)
    "all":          ["H-2-Db", "H-2-Kb",
                     "H-2-Dd", "H-2-Kd", "H-2-Ld",
                     "H-2-Dk", "H-2-Kk"],                      # 7  (maximum consolidated set)
}

# ── Variant search ─────────────────────────────────────────────────────────────
VARIANTS_MAX_RESULTS          = 100

# ── NCBI Entrez ────────────────────────────────────────────────────────────────
ENTREZ_EMAIL                  = ""    # set by user on project creation
ENTREZ_DATABASE               = "protein"
