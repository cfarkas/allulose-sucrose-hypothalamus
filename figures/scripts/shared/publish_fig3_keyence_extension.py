#!/usr/bin/env python3
"""Publish validated Figure 3 outputs; preserve reviewed A-D and unrelated figures."""
from pathlib import Path
import argparse, hashlib, json, re, shutil
from datetime import datetime, timezone
import pandas as pd
from PIL import Image, ImageChops
Image.MAX_IMAGE_PIXELS = None
PAPER = Path(__file__).resolve().parents[2]
PREFIX = 'Figure_3_cFos_NPY'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=PAPER)
    parser.add_argument('--analysis-dir', type=Path, default=PAPER/'analyses/fig3_keyence_extension_20260920')
    parser.add_argument('--statistics-dir', type=Path, help='Endpoint-specific main-analysis tables; defaults to analysis directory')
    parser.add_argument('--original-cohort', action='store_true', help='Publish only the original nine-animal experiment')
    parser.add_argument('--extended-water-cohort', action='store_true', help='Publish the cohort extended with the 2026 Water controls')
    parser.add_argument('--figure4-reference-json', type=Path, help='Figure4 hashes captured at the start of this run')
    args = parser.parse_args()
    root = args.root.resolve(); analysis = args.analysis_dir.resolve()
    statistics = args.statistics_dir.resolve() if args.statistics_dir else analysis
    single_cohort = args.original_cohort or args.extended_water_cohort
    iqr_main = statistics != analysis or single_cohort
    provenance_prefix = ('extended_water_cohort_' if args.extended_water_cohort else
                         'original_cohort_' if args.original_cohort else 'keyence_extension_')
    if args.extended_water_cohort:
        # The composition is read from the bundle so that adding another animal needs no code change.
        recorded = read_json(analysis/'analysis_validation.json')['conditions']
        expected_conditions = {condition:int(recorded[condition]) for condition in ('Allulose','Sucrose','Water')}
    elif args.original_cohort:
        expected_conditions = {'Allulose':3,'Sucrose':3,'Water':3}
    else:
        expected_conditions = {'Allulose':5,'Sucrose':4,'Water':3}
    stage = analysis/'assembled/Fig3'; current = root/'Fig3'; archive = analysis/'archive/Fig3'
    validations = {name: read_json(analysis/name) for name in
                   ('analysis_validation.json','spatial_permanova_validation.json','panel_render_validation.json')}
    validations['spatial_permanova_validation.json'] = read_json(statistics/'spatial_permanova_validation.json')
    if iqr_main:
        assert (statistics/'eligibility.csv').exists() and (statistics/'outlier_feature_ledger.csv').exists()
    assert validations['analysis_validation.json']['status'] == 'PASS'
    assert validations['spatial_permanova_validation.json']['status'] == 'PASS'
    assert validations['panel_render_validation.json']['all_text_inside_canvas']
    for name, expected in validations['spatial_permanova_validation.json']['inputs'].items():
        assert digest(statistics/name) == expected, ('Stale statistics', name)
    rendered = read_json(analysis/'panel_render_receipt.json')
    assert validations['panel_render_validation.json']['render_receipt_sha256'] == digest(analysis/'panel_render_receipt.json'), 'Stale panel validation'
    for group in ('input_sha256','files'):
        for name, expected in rendered[group].items():
            assert digest(analysis/name) == expected, ('Stale renderer', name)
    assembly = read_json(stage/'provenance/complete_conditions_20260908.json')
    assert assembly['quantitative_analysis_changed'] and assembly['nonspatial_analysis_changed']
    assert set(assembly['panel_map']) == set('ABCDEFGHI')
    assert digest(stage/(PREFIX+'.pdf')) == next(r['sha256'] for r in assembly['figures'] if r['language']=='en')
    for figure in assembly['figures']:
        for panel in figure['panels']:
            assert digest(root/panel['source']) == panel['sha256'], ('Stale assembly', panel['source'])
    # Validate the established output contract before any canonical writes.
    import importlib.util
    spec = importlib.util.spec_from_file_location('figure_contract',root/'scripts/utilities/07_validate_figure_outputs.py')
    validator = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = validator; spec.loader.exec_module(validator)
    staged_validation = validator.validate('Fig3', stage)
    write_json(analysis/'logs/staged_validation.json', staged_validation)
    before = read_json(analysis/'before_change_hashes.json')
    figure4 = read_json(args.figure4_reference_json or analysis/'figure4_before_hashes.json')
    for name, expected in figure4.items():
        assert digest(root/name) == expected, ('Figure 4 changed', name)
    checks = []
    preserved = {}
    for letter in 'ABCD':
        for language in ('','_spanish'):
            for suffix in ('.pdf','.png'):
                relative = 'Fig3/panels/'+PREFIX+'_Panel_'+letter+language+suffix
                target = root/relative
                assert digest(target) == before[relative], ('Reviewed panel changed', relative)
                preserved[relative] = digest(target)
                if suffix == '.png':
                    with Image.open(stage/'panels'/target.name) as a, Image.open(target) as b:
                        identical = a.size == b.size and ImageChops.difference(a.convert('RGB'),b.convert('RGB')).getbbox() is None
                    assert identical, ('Staged reviewed panel differs', relative)
                    checks.append(dict(panel=relative,identical_pixels=True))
    # Preserve every overwritten source/provenance file on the first publication.
    already_published = (analysis/'publication_receipt.json').exists()
    def preserve(destination):
        saved = archive/destination.relative_to(current)
        if destination.is_file() and not saved.exists() and not already_published:
            saved.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(destination,saved)
    published = []
    def install(source, destination):
        preserve(destination); destination.parent.mkdir(parents=True,exist_ok=True)
        temporary = destination.with_name(destination.name+'.publishing')
        shutil.copy2(source, temporary); temporary.replace(destination)
        published.append(destination)
    def write_frame(frame, destination):
        preserve(destination); frame.to_csv(destination,index=False); published.append(destination)
    for source in current.joinpath('provenance').glob('*'):
        if source.is_file(): preserve(source)
    for source in stage.iterdir():
        if source.is_file(): install(source,current/source.name)
    for folder in ('panels','legends','provenance'):
        for source in (stage/folder).iterdir():
            if not source.is_file(): continue
            match = re.search(r'_Panel_([A-Z])(?:_|\.)',source.name)
            if folder=='panels' and match and match.group(1) in 'ABCD': continue
            install(source,current/folder/source.name)
    # Replace live spatial aliases as well as the lettered publication panels.
    spatial_source = analysis/'spatial_panels/Fig3'
    for language in ('en','es'):
        for source in (spatial_source/language).glob('*'):
            if source.is_file(): install(source,current/'spatial_panels'/language/source.name)
    for source in spatial_source.glob('Spatial_shape_definition_NPY*'):
        if source.is_file(): install(source,current/'shape_cartoon'/source.name)
    ring_directory = current/'ring_cartoon'
    if ring_directory.exists():
        for source in ring_directory.rglob('*'):
            if source.is_file():
                preserve(source)
                assert (archive/source.relative_to(current)).exists()
                source.unlink()
        for directory in sorted(ring_directory.rglob('*'),reverse=True):
            if directory.is_dir(): directory.rmdir()
        ring_directory.rmdir()
    sd = current/'source_data'
    counts = pd.read_csv(statistics/'animal_plot_values.csv')
    assert counts.groupby('condition').size().to_dict() == expected_conditions
    if single_cohort: assert not counts.animal.str.contains('keyence',case=False).any()
    install(statistics/'animal_plot_values.csv',sd/(PREFIX+'_animal_plot_values.csv'))
    old_counts = pd.read_csv(archive/'source_data'/f'{PREFIX}_spatial_animal_counts.csv').set_index('animal_id')
    full_counts = counts.copy()
    full_counts['experimental_unit'] = 'biological animal'
    full_counts['field_count'] = full_counts.sections
    full_counts['sample_id'] = full_counts.animal_id.map(old_counts.sample_id).fillna(full_counts.animal_id)
    full_counts['animal_id_status'] = full_counts.animal_id.map(old_counts.animal_id_status).fillna('cohort-specific Keyence identity verified')
    full_counts['cfos_npy_double_positive_nuclei'] = full_counts.double_positive_nuclei
    full_counts['cfos_npy_over_npy'] = full_counts.double_over_npy
    full_counts['cfos_npy_over_dapi'] = full_counts.double_positive_nuclei/full_counts.dapi_nuclei
    write_frame(full_counts,sd/(PREFIX+'_spatial_animal_counts.csv'))
    stats = pd.read_csv(statistics/'abundance_statistics.csv')
    stats['p_value_raw'] = stats.p_value
    stats['comparison'] = stats.group_a.fillna('Water|Sucrose|Allulose') + stats.group_b.fillna('').map(lambda s: '|'+s if s else '')
    stats['experimental_unit'] = 'biological animal'; stats['analysis_set'] = ('extended_water_cohort_after_IQR_check' if args.extended_water_cohort else 'original_nine_after_IQR_check' if args.original_cohort else ('within_condition_IQR_filtered' if iqr_main else 'all_twelve_animals'))
    if 'outliers_excluded' not in stats: stats['outliers_excluded'] = False if single_cohort else iqr_main
    write_frame(stats[stats.endpoint.ne('npy_over_dapi')],sd/(PREFIX+'_anova_mwu_statistics_used.csv'))
    for name in ('abundance_robust_statistics.csv','abundance_robust_plan.json'):
        source = statistics/name
        if source.exists(): install(source,sd/(PREFIX+'_'+name))
    supplemental = counts[['animal','condition','sample','platform','biological_cohort','npy_over_dapi']].copy()
    supplemental['npy_over_dapi_percent'] = supplemental.npy_over_dapi*100
    write_frame(supplemental,sd/(PREFIX+'_supplemental_NPY_DAPI_animal_values.csv'))
    write_frame(stats[stats.endpoint.eq('npy_over_dapi')],sd/(PREFIX+'_supplemental_NPY_DAPI_anova_mwu_statistics.csv'))
    supplemental_readme = sd/(PREFIX+'_supplemental_NPY_DAPI_README.txt')
    preserve(supplemental_readme)
    composition_text = '/'.join(str(expected_conditions[c]) for c in ('Water','Sucrose','Allulose'))
    supplemental_readme.write_text((f'Supplementary NPY/DAPI abundance for the {sum(expected_conditions.values())} cohort animals (Water/Sucrose/Allulose n={composition_text}).' if single_cohort else 'Supplementary NPY/DAPI abundance for all 12 animals (Water/Sucrose/Allulose n=3/4/5).') + ' Not a panel in the master. Animal points are count-pooled within animal, with no excluded outliers. ANOVA and exact two-sided Mann-Whitney p-values are nominal, unadjusted and exploratory across acquisition/cohort. See Figure_3_cFos_NPY_spatial_shape_README.md for methods and limitations.\n')
    published.append(supplemental_readme)
    for name in ('animal_metrics','unit_metrics','condition_summary','stratum_metrics','cluster_counts_by_condition','contour_manifest','bandwidth_sensitivity_condition_summary','treatment_contrasts'):
        source = statistics/(name+'.csv')
        if not source.exists(): source = analysis/(name+'.csv')
        if not source.exists(): continue
        target_name = PREFIX+'_spatial_shape_'+name+'.csv'
        if iqr_main and not single_cohort and name == 'bandwidth_sensitivity_condition_summary': target_name = PREFIX+'_unfiltered_spatial_shape_'+name+'.csv'
        install(source,sd/target_name)
    install(statistics/'cluster_membership.csv.gz',sd/(PREFIX+'_spatial_shape_cluster_membership.csv.gz'))
    for pattern in ('spatial_permanova*.csv*','spatial_sugar_platform*.csv*','spatial_dispersion*.csv','spatial_feature_eligibility.csv','animal_roster.csv'):
        for source in statistics.glob(pattern): install(source,sd/(PREFIX+'_'+source.name))
    if iqr_main:
        for pattern in ('eligibility.csv','outlier_feature_ledger.csv','spatial_eligibility.csv','*summary*.csv','*comparison*.csv'):
            for source in statistics.glob(pattern): install(source,sd/(PREFIX+'_iqr_'+source.name))
        if not single_cohort:
            for name in ('animal_plot_values.csv','abundance_statistics.csv','spatial_permanova.csv','spatial_dispersion.csv','spatial_sugar_platform_sensitivity.csv','animal_metrics.csv','unit_metrics.csv','condition_summary.csv','cluster_counts_by_condition.csv'):
                install(analysis/name,sd/(PREFIX+'_unfiltered_'+name))
        stale = sd/(PREFIX+'_spatial_shape_bandwidth_sensitivity_condition_summary.csv')
        if stale.exists() and not single_cohort:
            assert digest(stale) == digest(analysis/'bandwidth_sensitivity_condition_summary.csv')
            stale.unlink()
    install(analysis/'README.md',sd/(PREFIX+'_spatial_shape_README.md'))
    if not single_cohort:
        install(analysis/'inputs/keyence_association_audit.json',sd/(PREFIX+'_keyence_association_audit.json'))
    # Retire the superseded inference and explicitly label the nine-animal model audit.
    retired = [] if single_cohort else [sd/(PREFIX+'_spatial_shape_treatment_contrasts.csv')]
    if single_cohort:
        retired.extend(current.joinpath('provenance').glob('keyence_extension_*'))
        retired.extend(sd.glob('*keyence*'))
        retired.extend(sd.glob('*_unfiltered_*'))
    if args.extended_water_cohort:
        retired.extend(current.joinpath('provenance').glob('original_cohort_*'))
        retired.extend(current.joinpath('provenance').glob('water4_cohort_*'))
        # The joint two-measure test and its dispersion diagnostic are no longer produced.
        retired.extend(sd.glob('*spatial_dispersion*'))
        retired.extend(sd.glob('*spatial_permanova_per_feature*'))
    legacy_audit = sd/(PREFIX+'_spatial_segmentation_model_setting_audit.csv')
    if legacy_audit.exists():
        preserve(legacy_audit)
        install(legacy_audit,sd/(PREFIX+'_legacy_nine_segmentation_model_setting_audit.csv'))
        retired.append(legacy_audit)
    for name in ('Figure_3_cFos_NPY_provenance.json','Figure_3_cFos_NPY_manifest.csv','REPRODUCED_CONTENT_MANIFEST.json','README_REPRODUCTION_PASS.json'):
        retired.append(current/'provenance'/name)
    for path in retired:
        if path.exists():
            # A file first published by a later run was never captured by the
            # first-publication snapshot, so archive it now rather than lose it.
            saved = archive/path.relative_to(current)
            if not saved.exists():
                saved.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(path,saved)
            assert saved.exists(), ('Missing archive',path)
            path.unlink()
    for name in ('analysis_metadata.json','analysis_validation.json','spatial_permanova_plan.json','spatial_permanova_validation.json','spatial_platform_diagnostics.json','panel_render_receipt.json','panel_render_validation.json'):
        source = statistics/name if name.startswith('spatial_') else analysis/name
        install(source,current/'provenance'/(provenance_prefix+name))
    if iqr_main:
        for source in statistics.glob('*.json'):
            install(source,current/'provenance'/(provenance_prefix+'iqr_'+source.name))
        if (analysis/'allulose_example_provenance.json').exists():
            install(analysis/'allulose_example_provenance.json',current/'provenance/keyence_extension_allulose_example.json')
    if (analysis/'examples/example_allulose_npy_provenance.json').exists():
        install(analysis/'examples/example_allulose_npy_provenance.json',current/'provenance'/(provenance_prefix+'allulose_example.json'))
        install(analysis/'examples/example_allulose_npy.npz',sd/(PREFIX+'_allulose_example.npz'))
    # Existing generic shape metadata now refers to this current cohort.
    install(analysis/'analysis_metadata.json',current/'provenance/spatial_shape_analysis_metadata.json')
    install(analysis/'analysis_validation.json',current/'provenance/spatial_shape_validation.json')
    if (analysis/'source_records/fixed_keyence_inputs.json').exists():
        install(analysis/'source_records/fixed_keyence_inputs.json',current/'provenance/keyence_extension_fixed_inputs.json')
    readme = current/'README.txt'; preserve(readme)
    if args.extended_water_cohort:
        water_n = expected_conditions['Water']
        added = read_json(analysis/'analysis_metadata.json').get('added_animals', [])
        spatial_counts = pd.read_csv(statistics/'inclusion_counts.csv')
        npy_water = spatial_counts[(spatial_counts.domain=='spatial') & (spatial_counts.endpoint=='marker_cfos') & (spatial_counts.condition=='Water')]
        estimable = int(npy_water.n_included.iloc[0]) if len(npy_water) else water_n
        readme.write_text(f"""FIGURE 3 — NPY-ASSOCIATED c-FOS LABELING AND SPATIAL CLUSTERING

Figure 3 compares {sum(expected_conditions.values())} animals: {water_n} Water, {expected_conditions['Sucrose']} Sucrose and {expected_conditions['Allulose']} Allulose.

Panels
  A    Representative NPY-GFP and c-FOS microscopy, merges and enlarged fields.
  B–D  DAPI nuclear maps with NPY/c-FOS assignments for the three conditions.
  E    c-FOS-positive nuclei as a fraction of DAPI nuclei.
  F    c-FOS/NPY double-positive nuclei as a fraction of NPY-positive nuclei.
  G    Positive-cell pair-counting illustration in Allulose field E8_FR6-4.
  H    Spatial clustering of total c-FOS-positive nuclei.
  I    Spatial clustering of NPY-associated c-FOS-positive nuclei.

Abundance fractions pool counts within each animal and weight animals equally.
The analysis includes an exactly enumerated Welch omnibus statistic, a
Brown–Forsythe variance test, exact studentized contrasts, Hodges–Lehmann shifts
and Cliff's delta. Panel E/F values and comparison definitions are provided in
source_data/Figure_3_cFos_NPY_abundance_robust_statistics.csv. IQR flags do not
exclude animals from inference.

Panels H/I count positive-cell pairs within acquisitions and tissue sides at
distances of 20–150 µm. Conditional random labeling fixes cell positions,
strata and positive counts. Each curve expresses observed clustering in
standard deviations from its own null. A studentized maximum-deviation global
envelope compares conditions across the distance range; Holm adjustment covers
the two endpoints. Distances measure proximity, without establishing anatomical
homology or connectivity. The independent unit for these comparisons is the
animal; cage membership is undocumented.

Total c-FOS analysis uses all {sum(expected_conditions.values())} animals. NPY-associated curves require at least
two double-positive nuclei and use {estimable}/{expected_conditions['Sucrose']}/{expected_conditions['Allulose']} Water/Sucrose/Allulose animals.
WATER_NPY3 has zero and WATER_NPY4 one such nucleus; both contribute abundance
measurements. WATER_NPY5 has two and contributes to the spatial comparison.

WATER_NPY4–5 and the sugar animals share a Zeiss LSM 780 and 25×/0.8 objective.
The former images use a 2×2 mean reduction to match the sugar images' analysis
pitch of 0.664213 µm. The water group also includes Leica acquisitions.
Physical calibration does not remove staining, detection or sampling differences.

Files
  Figure_3_cFos_NPY.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/       Individual A–I panels in English and Spanish, PNG and PDF.
  legends/      Figure and panel captions in both languages.
  source_data/  Plotted values, tests and portable cell-coordinate inputs.
  provenance/   Cohort, acquisition and figure validation records.

Recompute neuronal statistics from the GitHub repository without microscopy:
  python figures/scripts/shared/recheck_current_statistics.py

Rebuild and install the figure from a complete local Paper/ analysis tree:
  ./scripts/run_fig3_new_water_cohort.sh
The runner requires analyses/fig3_water_extended_20260924 and the accepted
microscopy, mask and panel inputs. --from-raw also extracts CZI channels and
runs GPU Cellpose segmentation. The downloadable Zenodo data do not include
the raw WATER_NPY4–5 acquisitions; GitHub supplies their reviewed cell-coordinate
inputs for numerical reproduction. Fig3/07_rebuild_complete_figure3.sh routes
to the same cohort workflow.
""")
    else:
        readme.write_text("""FIGURE 3 — ORIGINAL NPY EXPERIMENT, CLOSED DENSITY CONTOURS AND SPATIAL CLUSTERS

Figure_3_cFos_NPY.pdf and its 600-dpi PNG are the complete English master.
A–I panels and legends are available individually in English and Spanish.

Only the original nine animals are included: 3 Water, 3 Sucrose, 3 Allulose.
Existing reviewed DAPI/NPY/c-FOS ROIs are fixed. No original animal is flagged
by the within-condition 1.5-IQR rule. Water3 has zero double-positive cells;
its abundance stays included and its conditional spatial metrics are undefined.
A–D retain their reviewed files. G uses Allulose E8_FR6-4, the accepted tissue
and third-ventricle contours, and the median eminence at the bottom.
H–I show nominal_p and Holm symbol labels on each comparison bracket; numeric
p-values and PERMANOVA details are in bilingual legends and source tables.
NS means p>=0.05 for the indicated nominal or adjusted value; NE means
not estimable. Cluster comparisons use observed-minus-expected excess.

Rebuild and publish from the Paper root:
  ./scripts/run_fig3_original_cohort.sh
The usual Fig3/07_rebuild_complete_figure3.sh entry routes to this analysis.

Methods, original-only inputs, IQR ledger and results:
  ../analyses/fig3_original_cohort_20260920/README.md
Current receipts: provenance/original_cohort_*.json
Previous canonical files are archived in that analysis bundle's archive/Fig3.
Dated older provenance files describe historical revisions.
""")
    published.append(readme)
    for relative, expected in preserved.items(): assert digest(root/relative)==expected
    for relative, expected in figure4.items(): assert digest(root/relative)==expected
    # The assembler checked a snapshot; publication deliberately updates tables.
    published_assembly_path = current/'provenance/complete_conditions_20260908.json'
    published_assembly = read_json(published_assembly_path)
    if 'source_data_sha256_unchanged' in published_assembly:
        published_assembly['source_data_sha256_at_assembly_start'] = published_assembly.pop('source_data_sha256_unchanged')
    published_assembly['assembly_only_source_data_unchanged'] = True
    published_assembly['source_data_sha256_after_publication'] = {str(p.relative_to(root)):digest(p) for p in sorted(sd.iterdir()) if p.is_file()}
    published_assembly['publication_analysis_directory'] = str(analysis.relative_to(root))
    write_json(published_assembly_path,published_assembly)
    canonical_validation = validator.validate('Fig3', current)
    write_json(analysis/'logs/canonical_validation.json',canonical_validation)
    write_json(analysis/'nonspatial_panel_comparison.json',checks)
    receipt = dict(schema=('fig3-extended-water-cohort-publication-v1' if args.extended_water_cohort else 'fig3-original-cohort-publication-v1' if args.original_cohort else 'fig3-keyence-publication-v1'),status='PASS',
                   published_utc=datetime.now(timezone.utc).isoformat(),
                   conditions={condition:expected_conditions[condition] for condition in ('Water','Sucrose','Allulose')},
                   figure4_unchanged=True,figure4_hashes=figure4,
                   reviewed_A_D_files_preserved=preserved,
                   staged_reviewed_panel_pixels_identical=True,
                   analysis_directory=str(analysis.relative_to(root)),
                   source_masks_unchanged=True,
                   main_statistics_directory=str(statistics.relative_to(root)),
                   outlier_policy='within-condition 1.5-IQR exclusions by endpoint' if iqr_main else 'none',
                   published_files={str(p.relative_to(root)):digest(p) for p in sorted(set(published))},
                   output_contract=canonical_validation)
    write_json(analysis/'publication_receipt.json',receipt)
    shutil.copy2(analysis/'publication_receipt.json',current/'provenance'/(provenance_prefix+'publication_receipt.json'))
    print(json.dumps(dict(status='PASS',published_files=len(receipt['published_files']),reviewed_files_preserved=len(preserved),figure4_unchanged=True,master=str(current/(PREFIX+'.pdf'))),indent=2))

if __name__=='__main__':
    main()
