#!/usr/bin/env python3
"""Replace Figure 4 panels I and J with the ARC scale-resolved envelope test.

Panels I and J previously showed two scalar spatial summaries, core enrichment
inside the 80% POMC density contour and the clustered fraction, in ARC and in
ME. Those summaries answer "are activated cells closer together than chance?"
only at the one contour level and the one neighbour cap they were built from,
and the ME rows were largely non-estimable. Both panels now report, for every
cage, how strongly activated cells cluster at each distance from 20 to 150
micrometres against that animal's own exact random-labelling null, and compare
conditions with a studentized maximum-deviation global envelope test whose
single p-value already covers the whole distance range. This is the same
analysis Figure 3 panels H and I carry; only the population, the region and the
experimental unit differ. ME is not analysed and no ME row is drawn.

Panels A-H, their standalone files and every source-data table are preserved
byte for byte, and the master above the I/J row is checked to be pixel
identical to the version this revision archived. The scalar contour and cluster
tables remain in Fig4/source_data as before.

Run from the Paper root:

    python -B Fig4/09_relayout_arc_scale_resolved.py            # stage only
    python -B Fig4/09_relayout_arc_scale_resolved.py --publish  # stage, verify, install
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import sys

os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import fitz
import numpy as np
from PIL import Image

# Locate the Paper root from this file's own position, so the script works
# from Fig4/, from its scripts/Fig4/ mirror, and from any unpacked copy of the
# published tree whatever it is called.
HERE = Path(__file__).resolve()
PAPER = next((path for path in HERE.parents
              if path.name == 'Paper' or ((path / 'scripts/setup').is_dir() and (path / 'README.txt').is_file())), None)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f'Could not locate the Paper root above {HERE}')
DEFAULT_REVISION = PAPER / 'analyses/fig4_scale_resolved_20260925'
SHARED = PAPER / 'scripts/shared'

# The two panels this revision touches, with the standalone file suffix the
# published tree already uses and the native stem the renderer writes.
SPATIAL_PANELS = (('I', 'spatial_cfos', 'Figure4_Spatial_cFOS_occurrence'),
                  ('J', 'spatial_cfos_pomc', 'Figure4_Spatial_cFOS_POMC_occurrence'))
# Directories copied wholesale when the pre-revision state is archived.
ARCHIVED_FOLDERS = ('panels', 'legends', 'provenance', 'source_data', 'spatial_panels')
ARCHIVED_FILES = ('Figure_4.pdf', 'Figure_4.png', 'README.txt')
# Geometry of the plot box on a standalone I/J panel page, in points. The
# letter sits in the left gutter above and left of this box; only the box
# height changes, because the new plot has its own aspect ratio.
PLOT_LEFT_PT = 14.28
PLOT_TOP_PT = 21.45
PLOT_RIGHT_MARGIN_PT = 2.80
PLOT_BOTTOM_MARGIN_PT = 2.88
# Legend files that carry a "(I) ..." or "(J) ..." line.
LEGEND_FILES = ('Figure_4_LEGEND.txt', 'Figure_4_LEGEND_spanish.txt', 'Figure_4_caption.txt',
                'Figure_ARC_ME_multipanel_LEGEND.txt', 'Figure_ARC_ME_multipanel_caption.txt',
                'Figure_ARC_ME_panel_I_spatial_cfos_LEGEND.txt',
                'Figure_ARC_ME_panel_I_spatial_cfos_LEGEND_spanish.txt',
                'Figure_ARC_ME_panel_J_spatial_cfos_pomc_LEGEND.txt',
                'Figure_ARC_ME_panel_J_spatial_cfos_pomc_LEGEND_spanish.txt')
# The comparison-bracket note that the 20 September display revision appended.
# The new panels carry no brackets, so the note is removed from every legend.
BRACKET_NOTE_OPENERS = ('Comparison brackets join Allulose to Water and Sucrose',
                        'Los corchetes comparan Alulosa con Agua y Sacarosa')

README = """FIGURE 4 — POMC/c-FOS IN ARC AND ME, WITH SCALE-RESOLVED SPATIAL PANELS (25 SEPTEMBER 2026)

Figure_4.pdf and its 600-dpi PNG are the complete English master.
Individual A–J panels and their legends are available in English and Spanish.

Panels A–H retain the original microscopy, masks, quantitative summaries and
standalone panel files, and still report ARC and ME separately.

Panels I and J are scale resolved. For every animal the number of pairs of
activated cells closer than r is counted within each acquisition and tissue
side, for r from 20 to 150 micrometres, and expressed in standard deviations
from that animal's own conditional random-labelling null. The null holds every
cell position, stratum and activated-cell count fixed and permutes only which
cells are activated, so its mean and variance are exact and the field geometry
cancels between the observed and the expected count; no edge correction enters.
The independent unit is the biological cage, so animal curves are averaged with
equal weight within a cage. Conditions are compared with a studentized
maximum-deviation global envelope test (Myllymäki et al. 2017, J. R. Stat. Soc.
B 79:381) over an exhaustive enumeration of cage-label allocations. Taking the
maximum over distances controls the error across the whole curve, so the single
p-value printed on each panel needs no further correction across distances;
Holm then covers the two endpoints. This is the analysis Figure 3 panels H and
I carry; only the population, the region and the experimental unit differ.

ME is not analysed in I and J. It retains 74 marker-associated cells with 28
activated across 11 animals, and only 6 of those animals have two or more
activated cells, which is too sparse to define a curve. Panels G and H still
report both regions.

The earlier scalar spatial summaries, core enrichment inside the 80% POMC
density contour and the clustered fraction on the six-neighbour graph, are no
longer drawn. Their tables are unchanged in source_data/, and the panels they
produced are archived with their comparison receipts under
../analyses/fig4_scale_resolved_20260925/archive/Fig4/.

Rebuild and publish Figure 4 from the Paper root:

  ./scripts/run_fig4_scale_resolved.sh

That script re-derives the cage curves and the envelope test from the archived
POMC cell table, re-renders the bilingual I/J panels, re-assembles the master
and installs it here. Without --publish on the final step, outputs stay in
analyses/fig4_scale_resolved_20260925/assembled/Fig4 and nothing here changes.

Scale-resolved methods, per-cage curves, envelope bands and test results:
  ../analyses/fig4_scale_resolved_20260925/README.md
Contour and cluster methods, parameters, condition summaries and tests:
  ../analyses/spatial_shape_20260920/README.md
Historical shell source tables and earlier README:
  ../analyses/spatial_shape_20260920/archive/Fig4/
The original spatial_reanalysis_20260919 package remains the input reference.
"""


SPATIAL_PANELS_NOTE = """Native Figure 4 spatial panels.

Figure4_Spatial_cFOS_occurrence(.pdf/.png) and
Figure4_Spatial_cFOS_POMC_occurrence(.pdf/.png), with their _spanish variants
and legends.json, are the scale-resolved panels drawn as I and J. They are
produced by scripts/shared/fig4_scale_resolved_panels.py from
analyses/fig4_scale_resolved_20260925.

The four *_comparisons.json files belong to the superseded contour-and-cluster
panels that carried significance brackets. The current panels carry no
brackets, so those receipts describe nothing in the published figure and are
retained only as provenance for the panels archived under
../../analyses/fig4_scale_resolved_20260925/archive/Fig4/panels/.
"""


def paper_relative(path) -> str:
    """Record paths relative to the Paper root, so receipts do not carry the
    directory this tree happened to occupy when they were written."""
    path = Path(path).resolve()
    try:
        return str(path.relative_to(PAPER))
    except ValueError:
        return str(path)


def sha256(path) -> str:
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def unique_letter(page, letter) -> fitz.Rect:
    words = [word for word in page.get_text('words') if word[4] == letter]
    if len(words) != 1:
        raise ValueError(f'Expected one panel letter {letter}, got {len(words)}')
    return fitz.Rect(words[0][:4])


def save_pdf_png(document, pdf, raster_size, dpi):
    document.save(pdf, garbage=4, deflate=True)
    pixmap = document[0].get_pixmap(dpi=dpi, alpha=False)
    raster = Image.frombytes('RGB', (pixmap.width, pixmap.height), pixmap.samples)
    # PDF rasterization rounds page edges upward; reconcile only white edge pixels.
    target_width, target_height = raster_size
    for crop in ((target_width, 0, raster.width, raster.height),
                 (0, target_height, raster.width, raster.height)):
        if crop[0] < crop[2] and crop[1] < crop[3]:
            edge = np.asarray(raster.crop(crop))
            if edge.size and edge.min() < 250:
                raise ValueError('Cannot trim nonwhite page-rounding edge')
    canvas = Image.new('RGB', raster_size, 'white')
    canvas.paste(raster.crop((0, 0, min(raster.width, target_width), min(raster.height, target_height))), (0, 0))
    canvas.save(Path(pdf).with_suffix('.png'), dpi=(dpi, dpi))
    return canvas


def archive_canonical(canonical: Path, archive: Path, manifest_path: Path) -> dict:
    """Copy the published Figure 4 aside once, then verify it on every later run.

    The archive is what the rebuilt master is diffed against, so it has to be
    the state that existed before this revision first ran. Creating it is a
    copy; nothing in Fig4 is moved or removed.
    """
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for item in manifest['files']:
            if sha256(archive / item['path']) != item['sha256']:
                raise ValueError('Archived pre-revision input changed: ' + item['path'])
        return manifest
    archive.mkdir(parents=True, exist_ok=True)
    for folder in ARCHIVED_FOLDERS:
        shutil.copytree(canonical / folder, archive / folder, dirs_exist_ok=True)
    for name in ARCHIVED_FILES:
        shutil.copy2(canonical / name, archive / name)
    files = [dict(path=str(path.relative_to(archive)), sha256=sha256(path))
             for path in sorted(archive.rglob('*')) if path.is_file()]
    manifest = dict(schema='figure4_arc_scale_resolved_prerevision_v1',
                    source=paper_relative(canonical), files=files)
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def render_native_panels(analysis_dir: Path, native: Path, dpi: int) -> dict:
    sys.path.insert(0, str(SHARED))
    spec = importlib.util.spec_from_file_location('fig4_scale_resolved_panels',
                                                  SHARED / 'fig4_scale_resolved_panels.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.generate(analysis_dir, native, dpi)


def rewrite_legends(folder: Path, legends: dict) -> list[str]:
    """Swap the (I) and (J) lines and drop the retired comparison-bracket note."""
    changed = []
    for name in LEGEND_FILES:
        path = folder / name
        lines = path.read_text(encoding='utf-8').splitlines()
        language = 'es' if 'spanish' in name else 'en'
        replaced = 0
        for index, line in enumerate(lines):
            for letter in ('I', 'J'):
                if line.startswith(f'({letter}) '):
                    lines[index] = f'({letter}) ' + legends[language][letter]
                    replaced += 1
        if not replaced:
            raise ValueError(f'No (I)/(J) legend line in {path}')
        for index, line in enumerate(lines):
            if line.startswith(BRACKET_NOTE_OPENERS):
                lines = lines[:index]
                break
        text = '\n'.join(lines).rstrip('\n') + '\n'
        path.write_text(text, encoding='utf-8')
        changed.append(name)
    return changed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--revision-dir', type=Path, default=DEFAULT_REVISION)
    parser.add_argument('--dpi', type=int, default=600)
    parser.add_argument('--skip-render', action='store_true',
                        help='Reuse the native panels already in <revision>/panels/Fig4.')
    parser.add_argument('--publish', action='store_true',
                        help='After validation, install the revised I/J assets, master, legends and provenance into Fig4.')
    args = parser.parse_args()

    revision = args.revision_dir.resolve()
    canonical = PAPER / 'Fig4'
    archive = revision / 'archive/Fig4'
    stage = revision / 'assembled/Fig4'
    native = revision / 'panels/Fig4'
    manifest_path = revision / 'archive_manifest.json'
    archive_canonical(canonical, archive, manifest_path)

    if args.skip_render:
        panel_receipt = json.loads((native / 'panel_receipt.json').read_text())
    else:
        panel_receipt = render_native_panels(revision, native, args.dpi)
    legends = json.loads((native / 'legends.json').read_text(encoding='utf-8'))

    stage.mkdir(parents=True, exist_ok=True)
    for folder in ARCHIVED_FOLDERS:
        shutil.copytree(archive / folder, stage / folder, dirs_exist_ok=True)
    for name in ARCHIVED_FILES:
        shutil.copy2(archive / name, stage / name)
    for path in native.iterdir():
        if path.is_file():
            shutil.copy2(path, stage / 'spatial_panels' / path.name)
    (stage / 'spatial_panels/README.txt').write_text(SPATIAL_PANELS_NOTE, encoding='utf-8')
    (stage / 'README.txt').write_text(README, encoding='utf-8')
    rewrite_legends(stage / 'legends', legends)

    master_source = fitz.open(archive / 'Figure_4.pdf')
    master_page = master_source[0]
    offsets, panel_receipts, rebuilt = {}, [], {}
    old_panel_height = new_panel_height = None
    for letter, suffix, stem in SPATIAL_PANELS:
        for language in ('en', 'es'):
            locale = '_spanish' if language == 'es' else ''
            name = f'Figure_ARC_ME_panel_{letter}_{suffix}{locale}.pdf'
            old = fitz.open(archive / 'panels' / name)
            old_page = old[0]
            old_letter = unique_letter(old_page, letter)
            if old_letter.x1 >= PLOT_LEFT_PT or old_letter.y0 >= PLOT_TOP_PT:
                raise ValueError(f'Panel {letter} letter does not sit in the expected gutter')
            plot_x0, plot_y0 = PLOT_LEFT_PT, PLOT_TOP_PT
            plot_x1 = old_page.rect.width - PLOT_RIGHT_MARGIN_PT
            source = fitz.open(native / (stem + locale + '.pdf'))
            plot_height = (plot_x1 - plot_x0) * source[0].rect.height / source[0].rect.width
            page_height = plot_y0 + plot_height + PLOT_BOTTOM_MARGIN_PT
            if old_panel_height is None:
                old_panel_height, new_panel_height = old_page.rect.height, page_height
            new = fitz.open()
            page = new.new_page(width=old_page.rect.width, height=page_height)
            # Preserve the existing letter vector and its physical size.
            letter_clip = fitz.Rect(0, 0, plot_x0 - .1, 26)
            page.show_pdf_page(letter_clip, old, 0, clip=letter_clip)
            plot_rect = fitz.Rect(plot_x0, plot_y0, plot_x1, plot_y0 + plot_height)
            page.show_pdf_page(plot_rect, source, 0, keep_proportion=True)
            with Image.open(archive / 'panels' / name.replace('.pdf', '.png')) as original_png:
                size = (original_png.width, round(page_height / 72 * args.dpi))
            canvas = save_pdf_png(new, stage / 'panels' / name, size, args.dpi)
            panel_receipts.append(dict(panel=letter, language=language,
                                       old_height_pt=old_page.rect.height, new_height_pt=page_height,
                                       width_pt=old_page.rect.width, png_size=list(size),
                                       plot_rect_pt=[plot_rect.x0, plot_rect.y0, plot_rect.x1, plot_rect.y1],
                                       isotropic_plot_scale=plot_rect.width / source[0].rect.width,
                                       native_source=paper_relative(native / (stem + locale + '.pdf')),
                                       native_sha256=sha256(native / (stem + locale + '.pdf'))))
            if language == 'en':
                master_letter = unique_letter(master_page, letter)
                offsets[letter] = (master_letter.x0 - old_letter.x0, master_letter.y0 - old_letter.y0)
                rebuilt[letter] = (new, canvas)

    cutoff = min(offset[1] for offset in offsets.values())
    if max(offset[1] for offset in offsets.values()) - cutoff > 1e-6:
        raise ValueError('Panels I and J are not on a common baseline in the master')
    new_master_height = master_page.rect.height - old_panel_height + new_panel_height
    master = fitz.open()
    page = master.new_page(width=master_page.rect.width, height=new_master_height)
    preserved = fitz.Rect(0, 0, master_page.rect.width, cutoff)
    page.show_pdf_page(preserved, master_source, 0, clip=preserved)
    for letter, (document, _) in rebuilt.items():
        x, y = offsets[letter]
        page.show_pdf_page(fitz.Rect(x, y, x + document[0].rect.width, y + document[0].rect.height), document, 0)
    master.save(stage / 'Figure_4.pdf', garbage=4, deflate=True)

    # Keep the A-H raster pixels unchanged as well as their original geometry.
    original_raster = Image.open(archive / 'Figure_4.png').convert('RGB')
    shrink_pixels = round((old_panel_height - new_panel_height) / 72 * args.dpi)
    master_raster = Image.new('RGB', (original_raster.width, original_raster.height - shrink_pixels), 'white')
    cut_pixel = round(cutoff / 72 * args.dpi)
    master_raster.paste(original_raster.crop((0, 0, original_raster.width, cut_pixel)), (0, 0))
    for letter, (_, canvas) in rebuilt.items():
        x, y = offsets[letter]
        master_raster.paste(canvas, (round(x / 72 * args.dpi), round(y / 72 * args.dpi)))
    master_raster.save(stage / 'Figure_4.png', dpi=(args.dpi, args.dpi))

    unchanged = []
    for path in sorted((archive / 'panels').glob('*')):
        if re.search(r'panel_[A-H]_', path.name):
            assert sha256(path) == sha256(stage / 'panels' / path.name), path.name
            unchanged.append(str((stage / 'panels' / path.name).relative_to(stage)))
    for path in sorted((archive / 'source_data').glob('*')):
        if path.is_file():
            assert sha256(path) == sha256(stage / 'source_data' / path.name), path.name
    # Compare the preserved upper master at identical rasterization settings.
    clip = fitz.Rect(0, 0, master_page.rect.width, cutoff - 1)
    old_pixels = master_page.get_pixmap(matrix=fitz.Matrix(2, 2), clip=clip, alpha=False)
    new_pixels = page.get_pixmap(matrix=fitz.Matrix(2, 2), clip=clip, alpha=False)
    assert old_pixels.width == new_pixels.width and old_pixels.height == new_pixels.height
    difference = int(np.max(np.abs(np.frombuffer(old_pixels.samples, np.uint8).astype(np.int16)
                                   - np.frombuffer(new_pixels.samples, np.uint8).astype(np.int16))))
    if difference > 1:
        raise AssertionError('Upper master vector appearance changed: ' + str(difference))
    for language in ('en', 'es'):
        for letter in ('I', 'J'):
            if 'ARC' not in legends[language][letter]:
                raise ValueError(f'Panel {letter} ({language}) legend must name the analysed region')
            if re.search(r'\bME\b', legends[language][letter]):
                raise ValueError(f'Panel {letter} ({language}) legend still reports ME')

    spatial_receipt = dict(schema='figure4_arc_scale_resolved_spatial_v1',
                           analysis_dir=paper_relative(revision), region='ARC', me_analysed=False,
                           unit='biological cage', panel_renderer='scripts/shared/fig4_scale_resolved_panels.py',
                           panel_receipt=panel_receipt,
                           shape_spatial_source=paper_relative(stage / 'spatial_panels'),
                           legend_source=paper_relative(stage / 'spatial_panels/legends.json'),
                           legend_sha256=sha256(stage / 'spatial_panels/legends.json'),
                           panels={f'{letter}_{language}': dict(
                               source=paper_relative(stage / 'spatial_panels' / f'{stem}{"_spanish" if language == "es" else ""}.png'),
                               sha256=sha256(stage / 'spatial_panels' / f'{stem}{"_spanish" if language == "es" else ""}.png'))
                               for letter, _, stem in SPATIAL_PANELS for language in ('en', 'es')})
    (stage / 'provenance/shape_spatial_source_manifest.json').write_text(json.dumps(spatial_receipt, indent=2) + '\n')

    receipt = dict(status='PASS',
                   scope='Only Figure 4 panels I and J, the master and the matching legends changed',
                   spatial_analysis='scale-resolved global envelope test in ARC; ME not analysed',
                   A_H_standalone_files_byte_identical=len(unchanged),
                   A_H_master_png_pixels_unchanged_above_y=cut_pixel,
                   A_H_master_pdf_max_pixel_difference_at_144dpi=difference,
                   A_H_physical_coordinates_and_scale_unchanged=True,
                   source_tables_byte_identical=True,
                   nonspatial_analysis_changed=False,
                   old_panel_height_pt=old_panel_height, new_panel_height_pt=new_panel_height,
                   master_size_pt=[page.rect.width, page.rect.height],
                   master_png_size=list(master_raster.size),
                   panels=panel_receipts,
                   scale_resolved_test_sha256=sha256(revision / 'scale_resolved_test.csv'),
                   scale_resolved_curves_sha256=sha256(revision / 'scale_resolved_curves.csv'),
                   scale_resolved_envelope_sha256=sha256(revision / 'scale_resolved_envelope.csv'),
                   render_script_sha256=sha256(SHARED / 'fig4_scale_resolved_panels.py'),
                   analysis_script_sha256=sha256(SHARED / 'fig3_scale_resolved.py'),
                   assembly_script_sha256=sha256(Path(__file__)))
    (stage / 'provenance/fig4_arc_scale_resolved_revision.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (revision / 'validation.json').write_text(json.dumps(receipt, indent=2) + '\n')

    if args.publish:
        for relative in unchanged:
            assert sha256(canonical / relative) == sha256(archive / relative), relative
        for path in (archive / 'source_data').glob('*'):
            if path.is_file():
                assert sha256(path) == sha256(canonical / 'source_data' / path.name), path.name
        publish_files = [stage / 'Figure_4.pdf', stage / 'Figure_4.png', stage / 'README.txt']
        publish_files.extend(path for path in (stage / 'panels').glob('*') if re.search(r'panel_[IJ]_', path.name))
        publish_files.extend(stage / 'legends' / name for name in LEGEND_FILES)
        publish_files.extend(stage / 'provenance' / name for name in
                             ('fig4_arc_scale_resolved_revision.json', 'shape_spatial_source_manifest.json'))
        publish_files.extend(path for path in (stage / 'spatial_panels').glob('*') if path.is_file())
        published = {}
        for source in publish_files:
            destination = canonical / source.relative_to(stage)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            assert sha256(source) == sha256(destination)
            published[str(destination.relative_to(PAPER))] = sha256(destination)
        for relative in unchanged:
            assert sha256(canonical / relative) == sha256(archive / relative), relative
        for path in (archive / 'source_data').glob('*'):
            if path.is_file():
                assert sha256(path) == sha256(canonical / 'source_data' / path.name), path.name
        publication = dict(status='PASS', published=True,
                           scope='Only Figure 4 I/J, the master and matching documentation',
                           preserved_A_H_files=len(unchanged), source_data_preserved_byte_identical=True,
                           files=published, revision_validation_sha256=sha256(revision / 'validation.json'),
                           script_sha256=sha256(Path(__file__)))
        (revision / 'publication_receipt.json').write_text(json.dumps(publication, indent=2) + '\n')
        shutil.copy2(revision / 'publication_receipt.json',
                     canonical / 'provenance/fig4_arc_scale_resolved_publication.json')
        print('Published ARC scale-resolved Figure 4:', len(published), 'files; A-H and source data unchanged')
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
