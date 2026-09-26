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
        readme.write_text(f"""FIGURE 3 — NPY EXPERIMENT WITH EXTENDED WATER CONTROLS, DENSITY CONTOURS AND CLUSTERS

Figure_3_cFos_NPY.pdf and its 600-dpi PNG are the complete English master.
A–I panels and legends are available individually in English and Spanish.

{sum(expected_conditions.values())} animals are included: {water_n} Water, {expected_conditions['Sucrose']} Sucrose, {expected_conditions['Allulose']} Allulose.
The added Water animals ({', '.join(added)}; specimens N132-5 and N132-3,
24 September 2026) were acquired on the same Zeiss LSM 780 and the same
LD LCI Plan-Apochromat 25x/0.8 objective as the sugar animals. Their native
0.332106 um sampling was reduced by an exact 2x2 block mean to the sugar
cohort's own 0.664213 um pitch before segmentation, so every animal is
segmented at one physical scale. The original nine animals' counts and spatial
values are byte-identical to the previous revision.

The Water group is far more dispersed than either sugar group, and its spread
coincides with acquisition batch. The abundance endpoints therefore report
Welch's robust F with an exactly enumerated permutation null as the primary
omnibus, a Brown–Forsythe test of equal variance, and exact studentized
permutation tests per contrast, with Hodges–Lehmann shifts and Cliff's delta as
effect sizes. The Welch F approximation, the classical equal-variance ANOVA and
the exact Mann–Whitney tests are retained beside them. Each exact test also
records the smallest p it can return at these group sizes. The E and F panels
print the enumerated permutation p and the exact Mann–Whitney p-values; the
rest of that family is in the source table.
See source_data/Figure_3_cFos_NPY_abundance_robust_statistics.csv.

Existing reviewed DAPI/NPY/c-FOS ROIs are fixed. No animal is flagged by the
within-condition 1.5-IQR rule. Water animals with fewer than two NPY/c-FOS
double-positive nuclei stay in abundance while their conditional NPY spatial
metrics are undefined; that endpoint has {estimable} of {water_n} estimable Water animals.
A–D retain their reviewed files. G uses Allulose E8_FR6-4, the accepted tissue
and third-ventricle contours, and the median eminence at the bottom.
H and I are scale resolved. Each shows, for every animal, the clustering of
activated cells at distances from 20 to 150 um, expressed in standard
deviations from that animal's own random-labelling null, together with a
studentized maximum-deviation global envelope test comparing conditions across
the whole range. One p-value covers every distance, so no radius is singled
out, and the plot shows at which distances the conditions separate.
See scale_resolved_plan.json in the analysis bundle.

Core enrichment and clustering excess are retained and tested separately in
source_data/Figure_3_cFos_NPY_spatial_permanova.csv. There is no joint
two-measure test, and no multivariate dispersion diagnostic: the single-measure
analogue is the Brown-Forsythe equal-variance test in
source_data/Figure_3_cFos_NPY_abundance_robust_statistics.csv.

E and F print the exact permutation p and the exact Mann-Whitney p-values.
Every displayed value is a number, not a significance symbol. NE means not
estimable.

Water acquisition remains confounded with condition and batch; no test removes
that. The within-batch Allulose–Sucrose contrast is the only unconfounded one.

Rebuild and publish from the Paper root:
  ./scripts/run_fig3_new_water_cohort.sh
Add --from-raw to re-derive the 2026 Water animals from their CZI files.
The usual Fig3/07_rebuild_complete_figure3.sh entry routes to this analysis.

Methods, cohort inputs, IQR ledger and results:
  ../analyses/fig3_water_extended_20260924/README.md
Current receipts: provenance/extended_water_cohort_*.json
The previous nine-animal revision remains in
  ../analyses/fig3_original_cohort_20260920/
Previous canonical files are archived in this analysis bundle's archive/Fig3.
Dated older provenance files describe historical revisions.
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
