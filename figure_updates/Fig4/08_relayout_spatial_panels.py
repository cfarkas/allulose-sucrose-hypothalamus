#!/usr/bin/env python3
"""Render taller Figure 4 I/J only, retaining the reviewed A-H publication assets."""
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
import pandas as pd
from PIL import Image

PAPER = Path(__file__).resolve().parents[1]
DEFAULT = PAPER / 'analyses/fig4_height_20260920'


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def unique_letter(page, letter):
    words = [word for word in page.get_text('words') if word[4] == letter]
    if len(words) != 1:
        raise ValueError(f'Expected one panel letter {letter}, got {len(words)}')
    return fitz.Rect(words[0][:4])


def visible_spatial_image(page):
    images = [item for item in page.get_image_info() if page.rect.contains(fitz.Rect(item['bbox']))]
    if len(images) != 1:
        raise ValueError('Expected one fully visible canonical spatial raster')
    return fitz.Rect(images[0]['bbox'])


def save_pdf_png(document, pdf, raster_size, dpi):
    document.save(pdf, garbage=4, deflate=True)
    pixmap = document[0].get_pixmap(dpi=dpi, alpha=False)
    raster = Image.frombytes('RGB', (pixmap.width, pixmap.height), pixmap.samples)
    # PDF rasterization rounds page edges upward; reconcile only white edge pixels.
    target_width, target_height = raster_size
    if raster.width > target_width:
        edge = np.asarray(raster.crop((target_width, 0, raster.width, raster.height)))
        if edge.size and edge.min() < 250:
            raise ValueError('Cannot trim nonwhite horizontal page-rounding edge')
    if raster.height > target_height:
        edge = np.asarray(raster.crop((0, target_height, raster.width, raster.height)))
        if edge.size and edge.min() < 250:
            raise ValueError('Cannot trim nonwhite vertical page-rounding edge')
    canvas = Image.new('RGB', raster_size, 'white')
    canvas.paste(raster.crop((0, 0, min(raster.width, target_width), min(raster.height, target_height))), (0, 0))
    canvas.save(pdf.with_suffix('.png'), dpi=(dpi, dpi))
    return canvas


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--revision-dir', type=Path, default=DEFAULT)
    parser.add_argument('--analysis-dir', type=Path, default=PAPER/'analyses/spatial_shape_20260920')
    parser.add_argument('--height-factor', type=float, default=1.33)
    parser.add_argument('--dpi', type=int, default=600)
    parser.add_argument('--publish', action='store_true', help='After validation, publish only revised I/J assets, the master, legends and provenance; preserve A-H and source data.')
    args = parser.parse_args()
    revision = args.revision_dir.resolve(); archive = revision/'archive/Fig4'
    stage = revision/'assembled/Fig4'; native = revision/'spatial_panels'
    if args.height_factor <= 1 or not np.isfinite(args.height_factor):
        raise ValueError('The spatial panel height factor must be finite and greater than one')
    if not (archive/'Figure_4.pdf').is_file():
        raise FileNotFoundError('A verified pre-revision archive is required: '+str(archive))
    before = json.loads((revision/'archive_manifest.json').read_text())
    for item in before['files']:
        if sha256(item['archive']) != item['sha256']:
            raise ValueError('Archived canonical input changed: '+item['archive'])
    source_paths = [args.analysis_dir/'unit_metrics.csv', args.analysis_dir/'treatment_contrasts.csv']
    source_hashes = {str(path): sha256(path) for path in source_paths}
    stage.mkdir(parents=True, exist_ok=True)
    for folder in ('panels', 'legends', 'provenance', 'source_data'):
        shutil.copytree(archive/folder, stage/folder, dirs_exist_ok=True)
    for path in archive.iterdir():
        if path.is_file() and path.suffix == '.txt':
            shutil.copy2(path, stage/path.name)

    renderer_path = PAPER/'scripts/shared/shape_spatial_figures.py'
    spec = importlib.util.spec_from_file_location('fig4_spatial_render', renderer_path)
    renderer = importlib.util.module_from_spec(spec); spec.loader.exec_module(renderer)
    renderer.set_style()
    units, contrasts = (pd.read_csv(path) for path in source_paths)
    for language in ('en', 'es'):
        for endpoint in ('cfos', 'marker_cfos'):
            renderer.plot_pomc(units, contrasts, endpoint, language, native, args.dpi, height_factor=args.height_factor)
    legends = json.loads((args.analysis_dir/'panels/Fig4/legends.json').read_text())
    note_en = ('Comparison brackets join Allulose to Water and Sucrose and report nominal_p from the existing unadjusted p-value and Holm from the existing Holm-adjusted p-value. Each uses the same symbol key: NS, p ≥ 0.05; *, p < 0.05; **, p < 0.01; ***, p < 0.001; ****, p < 0.0001. NE means no estimable comparison and is never treated as NS. The clustering brackets test observed-minus-expected cluster fraction (cluster_excess); plotted points retain the observed and random-label-expected fractions. Exact p-values and estimability reasons remain in the source-data table.')
    note_es = ('Los corchetes comparan Alulosa con Agua y Sacarosa y muestran nominal_p a partir del valor p existente sin ajuste y Holm a partir del valor p existente ajustado con Holm. Ambos usan los mismos símbolos: NS, p ≥ 0,05; *, p < 0,05; **, p < 0,01; ***, p < 0,001; ****, p < 0,0001. NE indica una comparación no estimable y nunca se interpreta como NS. Los corchetes de agrupación evalúan la fracción observada menos la esperada (cluster_excess); los puntos conservan las fracciones observada y esperada por etiquetado aleatorio. Los valores p exactos y los motivos de no estimabilidad se conservan en los datos fuente.')
    for language, note in [('en', note_en), ('es', note_es)]:
        for letter in ('I', 'J'):
            legends[language][letter] += ' '+note
    (native/'Fig4/legends.json').write_text(json.dumps(legends, indent=2, ensure_ascii=False)+'\n')
    for path in (stage/'legends').glob('*.txt'):
        if re.search(r'panel_[IJ]_', path.name) or path.name.startswith(('Figure_4_', 'Figure_ARC_ME_multipanel_')):
            path.write_text(path.read_text()+'\n\n'+(note_es if 'spanish' in path.name else note_en)+'\n')

    shutil.copytree(native/'Fig4', stage/'spatial_panels', dirs_exist_ok=True)
    spatial_receipt_path = stage/'provenance/shape_spatial_source_manifest.json'
    spatial_receipt = {'schema':'figure4_bottom_height_display_revision_v1','analysis_changed':False,'nonspatial_analysis_changed':False,'shape_spatial_source':str(stage/'spatial_panels'),'legend_source':str(stage/'spatial_panels/legends.json'),'legend_sha256':sha256(stage/'spatial_panels/legends.json'),'height_factor_only_I_J':args.height_factor,'panels':{}}
    for letter,stem in [('I','Figure4_Spatial_cFOS_occurrence'),('J','Figure4_Spatial_cFOS_POMC_occurrence')]:
        for language in ('en','es'):
            source = stage/'spatial_panels'/(stem+('_spanish' if language=='es' else '')+'.png')
            spatial_receipt['panels'][letter+'_'+language]={'source':str(source),'sha256':sha256(source)}
    spatial_receipt_path.write_text(json.dumps(spatial_receipt,indent=2)+'\n')
    multipanel_path = stage/'provenance/multipanel_source_manifest.csv'
    multipanel = pd.read_csv(multipanel_path,keep_default_na=False)
    multipanel['shape_spatial_source']=str(stage/'spatial_panels')
    multipanel['shape_spatial_receipt']=str(spatial_receipt_path)
    multipanel['shape_spatial_receipt_sha256']=sha256(spatial_receipt_path)
    multipanel['bottom_panel_height_factor']=args.height_factor
    multipanel['bottom_panel_display_only_revision']=True
    multipanel.to_csv(multipanel_path,index=False)
    translation_path = stage/'provenance/multipanel_spanish_translation_receipt.json'
    translation=json.loads(translation_path.read_text())
    translation['spanish_spatial_sources']=[spatial_receipt['panels'][letter+'_es']['source'] for letter in ('I','J')]
    translation['bottom_panel_height_factor']=args.height_factor
    translation_path.write_text(json.dumps(translation,indent=2,ensure_ascii=False)+'\n')
    readme=stage/'README.txt'
    text=readme.read_text()
    start=text.index('Rebuild and assemble both figures')
    end=text.index('Detailed methods',start)
    text=text[:start]+('Display revision: only panels I/J are 33% taller. A–H standalone files and\nupper-master geometry remain unchanged. Each comparison bracket shows\nnominal_p and Holm symbols from the existing raw and adjusted p-values;\nNE remains distinct when a test is not estimable.\nExact p-values and full methods remain in source data and bilingual legends.\n\nReproduce this layout from its verified archived inputs (no reanalysis):\n  /home/server/anaconda3/envs/lazyslide311/bin/python -B Fig4/08_relayout_spatial_panels.py --publish\nRun from Paper. Without --publish, outputs are staged in\nanalyses/fig4_height_20260920/assembled/Fig4 without changing canonical files.\n\n')+text[end:]
    readme.write_text(text)

    master_source = fitz.open(archive/'Figure_4.pdf')
    master_page = master_source[0]
    positions = {}; panel_receipts = []; new_panels = {}
    for letter, suffix, native_stem in [('I', 'spatial_cfos', 'Figure4_Spatial_cFOS_occurrence'), ('J', 'spatial_cfos_pomc', 'Figure4_Spatial_cFOS_POMC_occurrence')]:
        for language in ('en', 'es'):
            locale_suffix = '_spanish' if language == 'es' else ''
            name = f'Figure_ARC_ME_panel_{letter}_{suffix}{locale_suffix}.pdf'
            old = fitz.open(archive/'panels'/name); old_page = old[0]
            image_rect = visible_spatial_image(old_page)
            width, height = old_page.rect.width, old_page.rect.height * args.height_factor
            new = fitz.open(); page = new.new_page(width=width, height=height)
            # Preserve the existing letter vector and its physical size.
            letter_clip = fitz.Rect(0, 0, image_rect.x0 - .1, 26)
            page.show_pdf_page(letter_clip, old, 0, clip=letter_clip)
            source = fitz.open(native/'Fig4'/(native_stem+locale_suffix+'.pdf'))
            new_image_rect = fitz.Rect(image_rect.x0, image_rect.y0, image_rect.x1, image_rect.y0+image_rect.height*args.height_factor)
            page.show_pdf_page(new_image_rect, source, 0, keep_proportion=True)
            with Image.open(archive/'panels'/name.replace('.pdf','.png')) as original_png:
                size = (original_png.width, round(original_png.height*args.height_factor))
            canvas = save_pdf_png(new, stage/'panels'/name, size, args.dpi)
            panel_receipts.append({'panel':letter,'language':language,'old_width_pt':old_page.rect.width,'new_width_pt':width,'old_height_pt':old_page.rect.height,'new_height_pt':height,'height_ratio':height/old_page.rect.height,'png_size':list(size),'isotropic_plot_scale':new_image_rect.width/source[0].rect.width,'statistics_recomputed':False})
            if language == 'en':
                old_letter = unique_letter(old_page, letter); master_letter = unique_letter(master_page, letter)
                positions[letter] = (master_letter.x0-old_letter.x0, master_letter.y0-old_letter.y0)
                new_panels[letter] = (new, canvas, old_page.rect.height)
    old_height = max(item[2] for item in new_panels.values())
    increment = old_height*(args.height_factor-1)
    master = fitz.open(); page = master.new_page(width=master_page.rect.width, height=master_page.rect.height+increment)
    cutoff = min(position[1] for position in positions.values())
    preserved = fitz.Rect(0, 0, master_page.rect.width, cutoff)
    page.show_pdf_page(preserved, master_source, 0, clip=preserved)
    for letter, (document, canvas, _) in new_panels.items():
        x,y=positions[letter]
        page.show_pdf_page(fitz.Rect(x,y,x+document[0].rect.width,y+document[0].rect.height),document,0)
    master.save(stage/'Figure_4.pdf',garbage=4,deflate=True)
    # Keep A-H raster pixels unchanged as well as their original PDF geometry.
    original_raster = Image.open(archive/'Figure_4.png').convert('RGB')
    master_raster = Image.new('RGB',(original_raster.width,original_raster.height+round(increment/72*args.dpi)),'white')
    cut_pixel = round(cutoff/72*args.dpi)
    master_raster.paste(original_raster.crop((0,0,original_raster.width,cut_pixel)),(0,0))
    for letter, (_, canvas, _) in new_panels.items():
        x,y=positions[letter];master_raster.paste(canvas,(round(x/72*args.dpi),round(y/72*args.dpi)))
    master_raster.save(stage/'Figure_4.png',dpi=(args.dpi,args.dpi))
    unchanged=[]
    for path in sorted((archive/'panels').glob('*')):
        if re.search(r'panel_[A-H]_',path.name):
            target=stage/'panels'/path.name
            assert sha256(path)==sha256(target),str(target)
            unchanged.append(str(target.relative_to(stage)))
    for path in (archive/'source_data').glob('*'):
        if path.is_file():assert sha256(path)==sha256(stage/'source_data'/path.name),str(path)
    assert all(sha256(path)==value for path,value in source_hashes.items())
    # Compare the preserved upper master at identical PDF rasterization settings.
    clip=fitz.Rect(0,0,master_page.rect.width,cutoff-1)
    old_pixels=master_page.get_pixmap(matrix=fitz.Matrix(2,2),clip=clip,alpha=False)
    new_pixels=page.get_pixmap(matrix=fitz.Matrix(2,2),clip=clip,alpha=False)
    assert old_pixels.width==new_pixels.width and old_pixels.height==new_pixels.height
    old_array=np.frombuffer(old_pixels.samples,np.uint8);new_array=np.frombuffer(new_pixels.samples,np.uint8)
    maximum_difference=int(np.max(np.abs(old_array.astype(np.int16)-new_array.astype(np.int16))))
    receipt={'status':'PASS','scope':'Only Figure 4 I/J heights and comparison annotations changed','height_factor':args.height_factor,'A_H_standalone_files_byte_identical':len(unchanged),'A_H_master_png_pixels_unchanged_above_y':cut_pixel,'A_H_master_pdf_max_pixel_difference_at_144dpi':maximum_difference,'A_H_physical_coordinates_and_scale_unchanged':True,'source_tables_byte_identical':True,'no_segmentation_or_statistical_recalculation':True,'source_table_hashes':source_hashes,'panels':panel_receipts,'master_height_increase_pt':increment,'master_size_pt':[page.rect.width,page.rect.height],'master_png_size':list(master_raster.size),'comparison_source':str(source_paths[1]),'render_script_sha256':sha256(renderer_path),'assembly_script_sha256':sha256(Path(__file__)),'render_only_revision':True}
    if maximum_difference>1:raise AssertionError('Upper master vector appearance changed: '+str(maximum_difference))
    (stage/'provenance/fig4_bottom_height_revision.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (revision/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    comparison_count = {name:{'NS':0,'NE':0,'significant':0} for name in ('nominal_p','Holm')}
    for path in sorted((stage/'spatial_panels').glob('*_comparisons.json')):
        for annotation in json.loads(path.read_text())['annotations']:
            rows=contrasts[(contrasts.cohort=='POMC') & (contrasts.endpoint==annotation['endpoint']) & (contrasts.region==annotation['region']) & (contrasts.metric==annotation['tested_metric']) & (contrasts.group_a==annotation['group_a']) & (contrasts.group_b==annotation['group_b'])]
            assert len(rows)==1
            for display_name,value_key,symbol_key in [('nominal_p','p_raw','nominal_symbol'),('Holm','p_holm','holm_symbol')]:
                value=float(rows.iloc[0][value_key]);symbol=annotation[symbol_key]
                expected='NE' if not np.isfinite(value) else ('****' if value<.0001 else '***' if value<.001 else '**' if value<.01 else '*' if value<.05 else 'NS')
                assert symbol==expected
                if np.isfinite(value):assert value==annotation[value_key]
                else:assert annotation[value_key] is None
                comparison_count[display_name][symbol if symbol in ('NS','NE') else 'significant']+=1
            assert annotation['label']==f"nominal_p: {annotation['nominal_symbol']}; Holm: {annotation['holm_symbol']}"
    receipt['bilingual_comparison_labels_verified_against_existing_table']=comparison_count
    (revision/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (stage/'provenance/fig4_bottom_height_revision.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if args.publish:
        canonical=PAPER/'Fig4'
        for relative in unchanged:
            assert sha256(canonical/relative)==sha256(archive/relative),relative
        for path in (archive/'source_data').glob('*'):
            if path.is_file():assert sha256(path)==sha256(canonical/'source_data'/path.name),str(path)
        publish_files=[stage/'Figure_4.pdf',stage/'Figure_4.png',stage/'README.txt']
        publish_files.extend(path for path in (stage/'panels').glob('*') if re.search(r'panel_[IJ]_',path.name))
        publish_files.extend(path for path in (stage/'legends').glob('*.txt') if re.search(r'panel_[IJ]_',path.name) or path.name.startswith(('Figure_4_','Figure_ARC_ME_multipanel_')))
        publish_files.extend(stage/'provenance'/name for name in ('fig4_bottom_height_revision.json','shape_spatial_source_manifest.json','multipanel_source_manifest.csv','multipanel_spanish_translation_receipt.json'))
        publish_files.extend(path for path in (stage/'spatial_panels').glob('*') if path.is_file())
        published={}
        for source in publish_files:
            destination=canonical/source.relative_to(stage);destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source,destination);assert sha256(source)==sha256(destination)
            published[str(destination.relative_to(PAPER))]=sha256(destination)
        for relative in unchanged:assert sha256(canonical/relative)==sha256(archive/relative),relative
        for path in (archive/'source_data').glob('*'):
            if path.is_file():assert sha256(path)==sha256(canonical/'source_data'/path.name),str(path)
        publication={'status':'PASS','published':True,'scope':'Only Figure 4 I/J, master and matching display documentation','preserved_A_H_files':len(unchanged),'source_data_preserved_byte_identical':True,'files':published,'revision_validation_sha256':sha256(revision/'validation.json'),'script_sha256':sha256(Path(__file__))}
        (revision/'publication_receipt.json').write_text(json.dumps(publication,indent=2)+'\n')
        shutil.copy2(revision/'publication_receipt.json',canonical/'provenance/fig4_bottom_height_publication.json')
        print('Published validated bottom-row revision:',len(published),'files; A-H and source data unchanged')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
