"""integrate_data — pure pandas core.

Stacks every track's CURATE_MURINE_{track_id}.csv into a project-wide table,
derives organism/protein columns from project_config, and projects the
configurable VIEW. No Rich, no openpyxl, no input() — testable in isolation.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from utils.naming import COLUMN_PEPTIDE, get_step_filename, parse_track_id


# ── Column catalog (single source of truth) ───────────────────────────────────
#
# Each entry: (internal_name, display_header, group, default_on, required)
#   group: 'identity' | 'binding' | 'conservation' | 'coverage' | 'murine'
#   required: True → always present, not user-toggleable
#   default_on: True → pre-checked in customization prompt
#
# `coverage_*` is dynamic (one per population) and is appended at runtime by
# `build_column_catalog`. The Coverage entries here are placeholders that get
# expanded once we know what populations actually exist in the data.

_STATIC_CATALOG: list[tuple[str, str, str, bool, bool]] = [
    # identity (always on)
    ('organism',                       'Organism',                        'identity',     True,  True),
    ('protein',                        'Protein',                         'identity',     True,  True),
    ('peptide',                        'Peptide',                         'identity',     True,  True),
    # binding (default on)
    ('best_combined_percentile',       'Best percentile (consensus)',     'binding',      True,  False),
    ('num_alleles_united',             'HLA alleles bound',               'binding',      True,  False),
    ('alleles_united',                 'HLA alleles',                     'binding',      True,  False),
    # binding (optional)
    ('netmhcpan_el_percentile',        'NetMHCpan EL %ile (best allele)', 'binding',      False, False),
    ('mhcflurry_presentation_percentile', 'MHCflurry %ile (best allele)', 'binding',      False, False),
    ('calis_score',                    'Calis immunogenicity score',      'binding',      False, False),
    # conservation (default on)
    ('pct_identity_100',               'Conservation 100%',               'conservation', True,  False),
    ('pct_identity_threshold',         'Conservation @ threshold',        'conservation', True,  False),
    ('conservation_label',             'Conservation label',              'conservation', True,  False),
    # conservation (optional)
    ('pct_identity_90',                'Conservation 90%',                'conservation', False, False),
    ('pct_identity_80',                'Conservation 80%',                'conservation', False, False),
    ('n_mutations_safe',               'Mutations: safe (excellent+tol)', 'conservation', False, False),
    ('n_mutations_risky',              'Mutations: likely lost',          'conservation', False, False),
    # murine (default on)
    ('murine_label',                   'Murine binding',                  'murine',       True,  False),
    ('murine_best_percentile',         'Murine best %ile',                'murine',       True,  False),
    ('num_murine_alleles_bound',       'Murine alleles bound',            'murine',       True,  False),
]


def build_column_catalog(coverage_populations: list[str]) -> list[dict]:
    """Returns the full catalog of selectable VIEW columns, with coverage
    columns expanded for the populations present in the data.

    Each entry is a dict ready for the prompt UI:
        {name, header, group, default_on, required}
    """
    catalog: list[dict] = []
    for name, header, group, default_on, required in _STATIC_CATALOG:
        catalog.append({
            'name': name, 'header': header, 'group': group,
            'default_on': default_on, 'required': required,
        })
        if name == 'num_alleles_united':
            # Coverage columns slot in right after the HLA breadth column so
            # related metrics stay visually adjacent.
            for population in coverage_populations:
                col_name = f'coverage_{population}'
                catalog.append({
                    'name': col_name,
                    'header': f'Coverage ({population})',
                    'group': 'coverage',
                    'default_on': True,
                    'required': False,
                })
    return catalog


def default_column_selection(catalog: list[dict]) -> list[str]:
    """Subset of catalog entries that are on by default."""
    return [entry['name'] for entry in catalog if entry['default_on'] or entry['required']]


# ── Per-track loading ─────────────────────────────────────────────────────────

def load_track_full_table(
    intermediate_dir: Path,
    track_id: str,
    project_config: dict,
) -> pd.DataFrame | None:
    """Reads CURATE_MURINE_{track_id}.csv and prepends organism/protein columns.

    Returns None when the file is missing (the caller decides whether to warn-and-skip
    or error). Organism/protein come from project_config['tracks'][track_id] so the
    naming is canonical even for organisms whose labels contain dashes.
    """
    curate_csv = intermediate_dir / track_id / 'murine' / get_step_filename('CURATE_MURINE', track_id)
    if not curate_csv.exists():
        return None

    df = pd.read_csv(curate_csv)
    organism, protein = parse_track_id(track_id, project_config)

    df.insert(0, 'track_id', track_id)
    df.insert(1, 'organism', organism)
    df.insert(2, 'protein',  protein)
    return df


def load_anchor_aggregates(intermediate_dir: Path, track_id: str) -> pd.DataFrame | None:
    """Aggregates CONSERVATION_MUTATIONS_{track_id}.xlsx per peptide into
    n_mutations_safe + n_mutations_risky. Returns None when the file is absent.

    Only computed when the user opts into the anchor columns — the read is
    skipped otherwise to keep the default path cheap.
    """
    mutations_xlsx = intermediate_dir / track_id / 'conservation' / get_step_filename(
        'CONSERVATION_MUTATIONS', track_id, ext='xlsx',
    )
    if not mutations_xlsx.exists():
        return None
    try:
        mutations_df = pd.read_excel(mutations_xlsx)
    except Exception:
        return None
    if mutations_df.empty or 'peptide' not in mutations_df.columns or 'mhc_verdict' not in mutations_df.columns:
        return mutations_df.assign(n_mutations_safe=0, n_mutations_risky=0).head(0)

    safe_verdicts  = {'excellent_match', 'tolerated'}
    risky_verdicts = {'likely_lost'}
    grouped = mutations_df.groupby('peptide')['mhc_verdict'].agg(
        n_mutations_safe  = lambda series: int(series.isin(safe_verdicts).sum()),
        n_mutations_risky = lambda series: int(series.isin(risky_verdicts).sum()),
    ).reset_index()
    return grouped


# ── Project-wide stacking ─────────────────────────────────────────────────────

def stack_tracks(
    intermediate_dir: Path,
    project_config: dict,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Concatenates every track's CURATE_MURINE table.

    Returns (full_df, present_tracks, skipped_tracks). Skipped tracks are
    those without a CURATE_MURINE file — the caller writes them to audit.
    """
    track_ids = list(project_config.get('tracks', {}).keys())
    frames: list[pd.DataFrame] = []
    present: list[str] = []
    skipped: list[str] = []

    for track_id in track_ids:
        df = load_track_full_table(intermediate_dir, track_id, project_config)
        if df is None:
            skipped.append(track_id)
            continue
        frames.append(df)
        present.append(track_id)

    if not frames:
        empty = pd.DataFrame(columns=['track_id', 'organism', 'protein', COLUMN_PEPTIDE])
        return empty, present, skipped

    full_df = pd.concat(frames, ignore_index=True, sort=False)
    return full_df, present, skipped


def coverage_populations_in(full_df: pd.DataFrame) -> list[str]:
    """Returns the population names found in `coverage_<population>` columns,
    preserving column order so downstream UI is deterministic."""
    return [
        col.removeprefix('coverage_')
        for col in full_df.columns
        if col.startswith('coverage_')
    ]


def attach_anchor_aggregates(
    full_df: pd.DataFrame,
    intermediate_dir: Path,
    present_tracks: list[str],
) -> pd.DataFrame:
    """Left-joins n_mutations_safe / n_mutations_risky into full_df.

    Each track's mutations file is aggregated independently then concatenated
    (peptide can repeat across tracks but with different verdicts).
    """
    per_track_aggregates: list[pd.DataFrame] = []
    for track_id in present_tracks:
        agg = load_anchor_aggregates(intermediate_dir, track_id)
        if agg is None or agg.empty:
            continue
        agg = agg.copy()
        agg.insert(0, 'track_id', track_id)
        per_track_aggregates.append(agg)

    if not per_track_aggregates:
        return full_df.assign(n_mutations_safe=0, n_mutations_risky=0)

    anchor_df = pd.concat(per_track_aggregates, ignore_index=True)
    merged = full_df.merge(anchor_df, on=['track_id', 'peptide'], how='left')
    merged['n_mutations_safe']  = merged['n_mutations_safe'].fillna(0).astype(int)
    merged['n_mutations_risky'] = merged['n_mutations_risky'].fillna(0).astype(int)
    return merged


# ── VIEW projection ───────────────────────────────────────────────────────────

def project_view(
    full_df: pd.DataFrame,
    selected_columns: list[str],
    catalog: list[dict],
) -> pd.DataFrame:
    """Returns the VIEW dataframe: only `selected_columns` that actually exist
    in `full_df`, in catalog order (not user-typed order — keeps display stable).

    `track_id` is explicitly dropped from VIEW since identity in the union is
    organism + protein + peptide (per plan).
    """
    catalog_order = [entry['name'] for entry in catalog]
    ordered_selection = [c for c in catalog_order if c in selected_columns and c in full_df.columns]
    ordered_selection = [c for c in ordered_selection if c != 'track_id']
    return full_df[ordered_selection].copy()


def view_headers(selected_columns: list[str], catalog: list[dict]) -> dict[str, str]:
    """Maps internal column names → user-facing headers for the columns in VIEW."""
    catalog_index = {entry['name']: entry['header'] for entry in catalog}
    return {col: catalog_index.get(col, col) for col in selected_columns}
