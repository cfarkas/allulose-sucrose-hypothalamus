#!/usr/bin/env python3
"""Assemble one complete English Figure 3 with bilingual individual panels.

Use after 03_make_figure_3_cfos_npy.py. The accepted sucrose HIL polygon is
mandatory. Quantitative inputs retain their original D/E identities internally;
the final paper labels are A–I, including an explanatory ring cartoon. No analysis or source measurement is changed.
"""
from pathlib import Path
import argparse,hashlib,json,re,shutil,sys,importlib.util
SHARED = next(p for p in Path(__file__).resolve().parents if (p/"scripts/shared").is_dir())/"scripts/shared"
sys.path.insert(0,str(SHARED))
from spatial_ring_cartoon import save_cartoon, caption as ring_caption
import fitz
from PIL import Image

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def clean(source):
    doc=fitz.open(source);page=doc[0];oldtext=page.get_text()
    streams=sorted(hashlib.sha256(doc.xref_stream_raw(x[0])).hexdigest() for x in page.get_images(full=True))
    removed=[]
    # Historical certified sources retain an imprecise density label. Correct
    # only that display text in the current presentation, preserving all data.
    axis_replacements=[]
    for old,new in [("Equal-DAPI-density radial shell","DAPI-count quantile ring"),
                    ("Anillo radial de igual densidad DAPI","Anillo por cuantiles DAPI")]:
        for rect in page.search_for(old):
            axis_replacements.append((rect,new))
            page.add_redact_annot(rect,fill=False)
    for b in page.get_text('dict',flags=fitz.TEXTFLAGS_DICT & ~fitz.TEXT_PRESERVE_IMAGES)['blocks']:
        for line in b.get('lines',[]):
            for s in line['spans']:
                if re.fullmatch('[A-J]',s['text'].strip()) and s['size']>=14:
                    rect=fitz.Rect(s['bbox']);x=(rect.x0+rect.x1)/2;y=(rect.y0+rect.y1)/2
                    page.add_redact_annot(fitz.Rect(x-.1,y-.1,x+.1,y+.1),fill=False);removed.append(s['text'])
    if removed or axis_replacements:page.apply_redactions(images=0,graphics=0,text=0)
    for rect,new in axis_replacements:
        size=7.8
        x=(rect.x0+rect.x1-fitz.get_text_length(new,fontname="helv",fontsize=size))/2
        page.insert_text((x,rect.y1-2),new,fontname="helv",fontsize=size)

    assert streams==sorted(hashlib.sha256(doc.xref_stream_raw(x[0])).hexdigest() for x in page.get_images(full=True))
    number=lambda s:sorted(re.findall(r'(?<![A-Za-z])\d+(?:[.,]\d+)?',s))
    assert number(oldtext)==number(page.get_text())
    return doc,removed

LEGENDS = {'en': {'title': 'NPY-associated c-FOS labeling and spatial profiles in the arcuate region.', 'body': 'A, Representative NPY-GFP, c-FOS and merged images from water, sucrose and allulose animals. B–D, DAPI nuclear masks showing NPY and c-FOS assignments in representative water, sucrose and allulose fields, respectively. Black lines delineate ventricular boundaries reviewed manually against each field’s DAPI image; gray outlines show the DAPI-supported tissue envelope used for display. E, c-FOS-positive nuclei as a percentage of all DAPI nuclei in the complete field. F, c-FOS-positive NPY nuclei as a percentage of NPY-positive nuclei. E–F include three animals per condition; points represent animals, bars mean ± SD, and marked outliers were retained. One-way ANOVA and nominal exact two-sided Mann–Whitney p values are shown. G–H, c-FOS/DAPI and double-positive/DAPI occurrence across six quantile shells with approximately equal numbers of DAPI nuclei; thin traces represent animals, thick traces group means, and heat maps show shell means. The c-FOS profile did not separate globally (PERMANOVA p=0.121; q=0.121). The double-positive profile differed (p=0.0107; q=0.0214), but multivariate dispersion also differed (p=0.0143; q=0.0286). Shells represent normalized nuclear positions rather than anatomically registered layers, and the spatial double-positive denominator is DAPI rather than NPY. The finding therefore combines abundance and position and remains sensitive to dispersion and acquisition differences. Scale bars are retained on the corresponding microscopy panels.', 'panels': {'A': 'Representative NPY-GFP (green), c-FOS (red), merged and magnified merged images from water, sucrose and allulose animals. Dashed boxes locate the magnified fields; black outlines identify NPY ROIs assigned to double-positive DAPI nuclei. Scale bars, 100 µm in full fields and 50 µm in magnified fields. Display intensities are adjusted independently and do not support comparisons of brightness across conditions.', 'B': 'Water: DAPI nuclear masks displaying NPY-positive nuclei, c-FOS-positive nuclei and combined assignments. Visible NPY/c-FOS/double-positive counts are 205/113/57. Black ventricular boundaries were reviewed manually against the original registered DAPI. Gray external outlines are DAPI-supported tissue envelopes used only for presentation. Neither outline is an analytical ROI. Scale bars, 100 µm.', 'C': 'Sucrose, E7_FR7-5: DAPI nuclear masks displaying NPY-positive nuclei, c-FOS-positive nuclei and combined assignments. Visible NPY/c-FOS/double-positive counts are 292/66/19; full native-field counts are 302/81/19, with 3,241 DAPI nuclei. The display crop does not define the analytical denominator. Black ventricular boundaries were reviewed manually against the original registered DAPI. Gray external outlines are DAPI-supported tissue envelopes used only for presentation. Neither outline is an analytical ROI. Scale bars, 100 µm.', 'D': 'Allulose, E8_FR6-4: DAPI nuclear masks displaying NPY-positive nuclei, c-FOS-positive nuclei and combined assignments. Visible NPY/c-FOS/double-positive counts are 493/271/131. Black ventricular boundaries were reviewed manually against the original registered DAPI. Gray external outlines are DAPI-supported tissue envelopes used only for presentation. Neither outline is an analytical ROI. Scale bars, 100 µm.', 'E': 'c-FOS-positive nuclei as a percentage of all DAPI nuclei in each complete field. One-way ANOVA p=0.0951; exact two-sided Mann–Whitney p values: water–sucrose 0.2, water–allulose 1.0, sucrose–allulose 0.1. Three animals per condition. Points represent animals, bars mean ± SD; flagged outliers were retained. Pairwise p values are nominal, without multiplicity correction.', 'F': 'c-FOS-positive NPY nuclei as a percentage of NPY-positive nuclei in each complete field. One-way ANOVA p=0.0969; exact two-sided Mann–Whitney p values: water–sucrose 0.7, water–allulose 0.7, sucrose–allulose 0.1. Three animals per condition. Points represent animals, bars mean ± SD; flagged outliers were retained. Pairwise p values are nominal, without multiplicity correction.', 'G': 'c-FOS/DAPI occurrence across six shells defined from the DAPI nuclear distribution. Exact animal-level PERMANOVA p=0.121; BH q=0.121. The corresponding dispersion test gave p=0.0643; q=0.0643. Thin traces represent individual animals (three per condition), thick traces group means and heatmaps condition means. Shells are normalized nuclear-distribution coordinates, not anatomically registered layers. Cell counts and shells are not independent biological replicates. These occurrence profiles combine abundance and position; acquisition differences also constrain interpretation.', 'H': 'c-FOS/NPY double-positive occurrence relative to DAPI nuclei across six shells defined from the DAPI nuclear distribution. Exact animal-level PERMANOVA p=0.0107; BH q=0.0214. Dispersion also differed (p=0.0143; q=0.0286), so the result may reflect differences in centroid, dispersion or both. Thin traces represent individual animals (three per condition), thick traces group means and heatmaps condition means. Shells are normalized nuclear-distribution coordinates, not anatomically registered layers. Cell counts and shells are not independent biological replicates. These occurrence profiles combine abundance and position; acquisition differences also constrain interpretation.'}}, 'es': {'title': 'Marcación c-FOS asociada a NPY y perfiles espaciales en la región del núcleo arqueado.', 'body': 'A, Imágenes representativas de NPY-GFP, c-FOS y su combinación en animales de agua, sacarosa y alulosa. B–D, Máscaras nucleares DAPI con asignaciones NPY y c-FOS en campos de agua, sacarosa y alulosa, respectivamente. Las líneas negras delimitan los ventrículos revisados manualmente sobre DAPI; los bordes grises muestran el soporte de tejido estimado a partir de DAPI para la presentación. E, Porcentaje de núcleos c-FOS positivos respecto de DAPI en el campo completo. F, Porcentaje de núcleos NPY positivos con c-FOS. En E–F se incluyen tres animales por condición; los puntos corresponden a animales, las barras a media ± DE y los valores marcados como atípicos se conservaron. Se muestran ANOVA de una vía y valores p nominales de Mann–Whitney exacta bilateral. G–H, Ocurrencia c-FOS/DAPI y doble positiva/DAPI en seis capas definidas a partir de la distribución nuclear DAPI de cada campo; las líneas finas representan animales, las gruesas medias grupales y los mapas de calor medias por capa. El perfil c-FOS no difirió globalmente (PERMANOVA p=0,121; q=0,121). El perfil doble positivo difirió (p=0,0107; q=0,0214), junto con la dispersión multivariada (p=0,0143; q=0,0286). Las capas representan posiciones nucleares normalizadas, sin registro anatómico común; el denominador espacial doble positivo es DAPI. El resultado combina abundancia y posición y debe interpretarse considerando la dispersión y las diferencias de adquisición. Las barras de escala se indican en los paneles correspondientes.', 'panels': {'A': 'Imágenes representativas de NPY-GFP (verde), c-FOS (rojo), combinación y ampliación en agua, sacarosa y alulosa. Los recuadros discontinuos localizan las ampliaciones; los contornos negros señalan regiones NPY asignadas a núcleos DAPI dobles positivos. Barras de escala: 100 µm en campos completos y 50 µm en ampliaciones. La intensidad se ajusta por imagen para su visualización y no permite comparar brillo entre condiciones.', 'B': 'Agua: máscaras nucleares DAPI con asignaciones NPY, c-FOS y combinadas. Recuentos visibles NPY/c-FOS/dobles positivos: 205/113/57. Los contornos ventriculares negros se revisaron manualmente sobre el DAPI original registrado. Los bordes externos grises proceden del soporte de tejido estimado mediante DAPI y se utilizan solo para presentación; ninguno de estos contornos es una región de análisis. Barras de escala: 100 µm.', 'C': 'Sacarosa, E7_FR7-5: máscaras nucleares DAPI con asignaciones NPY, c-FOS y combinadas. Recuentos visibles NPY/c-FOS/dobles positivos: 292/66/19; campo nativo completo: 302/81/19, con 3.241 núcleos DAPI. El recorte de presentación no define el denominador analítico. Los contornos ventriculares negros se revisaron manualmente sobre el DAPI original registrado. Los bordes externos grises proceden del soporte de tejido estimado mediante DAPI y se utilizan solo para presentación; ninguno de estos contornos es una región de análisis. Barras de escala: 100 µm.', 'D': 'Alulosa, E8_FR6-4: máscaras nucleares DAPI con asignaciones NPY, c-FOS y combinadas. Recuentos visibles NPY/c-FOS/dobles positivos: 493/271/131. Los contornos ventriculares negros se revisaron manualmente sobre el DAPI original registrado. Los bordes externos grises proceden del soporte de tejido estimado mediante DAPI y se utilizan solo para presentación; ninguno de estos contornos es una región de análisis. Barras de escala: 100 µm.', 'E': 'Porcentaje de núcleos c-FOS positivos respecto de DAPI en el campo completo. ANOVA de una vía p=0,0951; Mann–Whitney exacta bilateral: agua–sacarosa p=0,2, agua–alulosa p=1,0 y sacarosa–alulosa p=0,1. Tres animales por condición. Los puntos representan animales y las barras media ± DE; se conservaron los valores atípicos señalados. Los valores p entre pares son nominales, sin corrección por multiplicidad.', 'F': 'Porcentaje de núcleos NPY positivos con c-FOS en el campo completo. ANOVA de una vía p=0,0969; Mann–Whitney exacta bilateral: agua–sacarosa p=0,7, agua–alulosa p=0,7 y sacarosa–alulosa p=0,1. Tres animales por condición. Los puntos representan animales y las barras media ± DE; se conservaron los valores atípicos señalados. Los valores p entre pares son nominales, sin corrección por multiplicidad.', 'G': 'Ocurrencia c-FOS/DAPI en seis capas definidas por la distribución de núcleos DAPI. PERMANOVA exacta por animal p=0,121; q de BH=0,121. La prueba de dispersión dio p=0,0643; q=0,0643. Las líneas finas representan animales (tres por condición), las gruesas medias grupales y los mapas de calor medias por condición. Las capas son coordenadas normalizadas de la distribución nuclear, sin registro anatómico común. Los núcleos y las capas no son réplicas biológicas independientes. Estos perfiles combinan abundancia y posición; las diferencias de adquisición también limitan su interpretación.', 'H': 'Ocurrencia doble positiva c-FOS/NPY respecto de DAPI en seis capas definidas por la distribución nuclear. PERMANOVA exacta por animal p=0,0107; q de BH=0,0214. La dispersión también difirió (p=0,0143; q=0,0286); el resultado puede reflejar diferencias de centroide, dispersión o ambas. Las líneas finas representan animales (tres por condición), las gruesas medias grupales y los mapas de calor medias por condición. Las capas son coordenadas normalizadas de la distribución nuclear, sin registro anatómico común. Los núcleos y las capas no son réplicas biológicas independientes. Estos perfiles combinan abundancia y posición; las diferencias de adquisición también limitan su interpretación.'}}}

for lang in ("en", "es"):
    panels=LEGENDS[lang]["panels"]
    panels["I"]=panels.pop("H")
    panels["H"]=panels.pop("G")
    panels["G"]=ring_caption("NPY",lang)
    LEGENDS[lang]["panels"]={letter:panels[letter] for letter in "ABCDEFGHI"}
    LEGENDS[lang]["body"]=LEGENDS[lang]["body"].replace("G–H,", "G, "+ring_caption("NPY",lang)+" H–I,")

def write_legends(out):
    folder=out/'legends';folder.mkdir(exist_ok=True)
    for lang,data in LEGENDS.items():
        suffix='' if lang=='en' else '_spanish';heading='Figure' if lang=='en' else 'Figura'
        body=heading+' 3. '+data['title']+'\n\n'+data['body']+'\n'
        (folder/f'Figure_3_cFos_NPY_LEGEND{suffix}.txt').write_text(body)
        if lang=='en':(out/'Figure_3_cFos_NPY_caption.txt').write_text(body)
        for letter,paragraph in data['panels'].items():
            (folder/f'Figure_3_cFos_NPY_Panel_{letter}_LEGEND{suffix}.txt').write_text(heading+' 3, panel '+letter+'. '+paragraph+'\n')

def refresh_spatial_panels(root, out, dpi):
    """Render current axis labels/layout from the certified numerical tables."""
    result={}
    for lang,suffix,folder in [('en','','spatial_english'),('es','_spanish','spatial_spanish')]:
        path=root/'Fig3'/f'05_analyze_spatial_distributions{suffix}.py'
        spec=importlib.util.spec_from_file_location('figure3_spatial_plot_'+lang,path)
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        source=root/'analyses/Fig3/results'/folder
        features=module.pd.read_csv(source/'animal_radial_occurrence_features.csv')
        summary=module.pd.read_csv(source/'exact_spatial_permanova.csv')
        target=out/'spatial_panels'/lang;target.mkdir(parents=True,exist_ok=True)
        for endpoint in module.ENDPOINTS:
            module.save_panel(features,summary,endpoint,target,dpi)
        result[lang]=target
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path.cwd());p.add_argument('--output-dir',type=Path);p.add_argument('--base-panels-dir',type=Path);p.add_argument('--dpi',type=int,default=600);args=p.parse_args()
    if args.dpi<300:raise ValueError('Publication exports require at least 300 dpi')
    root=args.root.resolve();out=args.output_dir or root/'Fig3';analysis=root/'analyses/Fig3/results';stage=analysis/'complete_conditions_20260908';quant=args.base_panels_dir or stage/'quantitative_inputs'
    if args.base_panels_dir:
        cache=stage/'quantitative_inputs';cache.mkdir(parents=True,exist_ok=True)
        for letter in 'DE':
            for suffix in ('','_spanish'):
                src=quant/f'Figure_3_cFos_NPY_Panel_{letter}{suffix}.pdf'
                if src.resolve() != (cache/src.name).resolve():shutil.copy2(src,cache/src.name)
        quant=cache
    hil=root/'Fig3/hil_review/sucrose_ventricle_20260908_v1';rec=json.loads((hil/'sucrose_ventricle_receipt.json').read_text());manifest=json.loads((hil/'review_inputs/REVIEW_INPUT_MANIFEST.json').read_text())
    assert rec['status']=='accepted' and rec['reviewer'] and rec['session']
    assert rec['mask_file_sha256']==sha(hil/'sucrose_ventricle_mask.png')
    assert rec['review_manifest_sha256']==sha(hil/'review_inputs/REVIEW_INPUT_MANIFEST.json')
    for path,digest in manifest['source_sha256'].items():assert sha(root/path)==digest,path
    cartoon_receipt=json.loads((stage/'sucrose_panels/receipt.json').read_text())
    assert cartoon_receipt['ventricle_hil_receipt_sha256']==sha(hil/'sucrose_ventricle_receipt.json')
    source_tables=list((root/'Fig3/source_data').glob('*.csv'))
    before={str(f.relative_to(root)):sha(f) for f in source_tables}
    out.mkdir(parents=True,exist_ok=True);(out/'panels').mkdir(exist_ok=True);(out/'provenance').mkdir(exist_ok=True)
    save_cartoon(out/'ring_cartoon', marker='NPY', dpi=args.dpi)
    spatial_render=refresh_spatial_panels(root,out,args.dpi)
    records=[]
    for lang,suffix,folder,spatial in [('en','','final_run','spatial_english'),('es','_spanish','spanish_analysis_work','spatial_spanish')]:
        sources={
          'A':analysis/folder/'panels/Figure3_Panel_A_raw_microscopy.pdf',
          'B':analysis/folder/'panels/Figure3_Panel_B_Water_marker_positive_nuclei.pdf',
          'C':stage/f'sucrose_panels/Sucrose_NPY_cartoon_{lang}.pdf',
          'D':analysis/folder/'panels/Figure3_Panel_C_Allulose_marker_positive_nuclei.pdf',
          'E':quant/f'Figure_3_cFos_NPY_Panel_D{suffix}.pdf',
          'F':quant/f'Figure_3_cFos_NPY_Panel_E{suffix}.pdf',
          'G':out/'ring_cartoon'/f'Spatial_ring_definition_NPY{suffix}.pdf',
          'H':spatial_render[lang]/'Figure3_Spatial_cFOS_occurrence.pdf',
          'I':spatial_render[lang]/'Figure3_Spatial_cFOS_NPY_occurrence.pdf'}
        cleaned={};source_records=[]
        for letter,source in sources.items():
            doc,removed=clean(source);cleaned[letter]=doc
            source_records.append({'panel':letter,'source':str(source.relative_to(root)),'sha256':sha(source),'removed_previous_letter':removed})
        width=1152;mar=18;gap=20;labelgap=25;y=mar;layout={}
        for row,ratios in [('A',[1]),('BCD',[1,1,1]),('EFG',[.26,.26,.48]),('HI',[1,1])]:
            available=width-2*mar-gap*(len(row)-1)
            cells=[available*r/sum(ratios) for r in ratios]
            heights=[cell*cleaned[letter][0].rect.height/cleaned[letter][0].rect.width for letter,cell in zip(row,cells)]
            x=mar
            for letter,cell,h in zip(row,cells,heights):
                layout[letter]=fitz.Rect(x,y+labelgap,x+cell,y+labelgap+h)
                x+=cell+gap
            y+=labelgap+max(heights)+gap
        master=fitz.open();page=master.new_page(width=width,height=y-gap+mar)
        for letter,box in layout.items():
            page.show_pdf_page(box,cleaned[letter],0);page.insert_text((box.x0,box.y0-7),letter,fontsize=20,fontname='hebo')
            # Standalone final panels retain the same complete source content.
            item=fitz.open();src=cleaned[letter][0];ip=item.new_page(width=src.rect.width,height=src.rect.height+labelgap)
            ip.show_pdf_page(fitz.Rect(0,labelgap,src.rect.width,src.rect.height+labelgap),cleaned[letter],0);ip.insert_text((4,18),letter,fontsize=18,fontname='hebo')
            target=out/'panels'/f'Figure_3_cFos_NPY_Panel_{letter}{suffix}.pdf';item.save(target,garbage=4,deflate=True)
            ip.get_pixmap(dpi=args.dpi,alpha=False).save(target.with_suffix('.png'))
        record={'language':lang,'panels':source_records,'panel_boxes_pt':{k:list(v) for k,v in layout.items()},'size_pt':[width,page.rect.height],'dpi':args.dpi,'panel_sequence':list('ABCDEFGHI')}
        if lang == 'en':
            target=out/'Figure_3_cFos_NPY.pdf';master.save(target,garbage=4,deflate=True)
            page.get_pixmap(dpi=args.dpi,alpha=False).save(target.with_suffix('.png'))
            page.get_pixmap(matrix=fitz.Matrix(1500/width,1500/width),alpha=False).save(stage/'Figure3_complete_en_preview.png')
            record.update(pdf=str(target.relative_to(root)) if target.is_relative_to(root) else str(target),sha256=sha(target))
        records.append(record)
        print('Completed English master and panels' if lang=='en' else 'Completed Spanish individual panels',flush=True)
    assert before=={str(f.relative_to(root)):sha(f) for f in source_tables}
    receipt={'schema':'figure3_complete_three_conditions_v1','reviewer':rec['reviewer'],'sucrose_hil_receipt':str((hil/'sucrose_ventricle_receipt.json').relative_to(root)),'sucrose_hil_sha256':sha(hil/'sucrose_ventricle_receipt.json'),'native_sucrose_counts':manifest['native_counts'],'panel_map':{'A':'microscopy: water/sucrose/allulose','B':'water nuclear assignments','C':'sucrose nuclear assignments','D':'allulose nuclear assignments','E':'whole-field c-FOS/DAPI','F':'c-FOS-positive/NPY-positive','G':'DAPI / third-ventricle quantile ring schematic','H':'c-FOS/DAPI spatial profile','I':'double-positive/DAPI spatial profile'},'publication_language':'en', 'bilingual_outputs':'individual panels and legends only', 'sucrose_cartoon_receipt':str((stage/'sucrose_panels/receipt.json').relative_to(root)), 'sucrose_cartoon_sha256':sha(stage/'sucrose_panels/receipt.json'), 'ring_cartoon':'DAPI covariance-normalized quantile shells; illustration only','analysis_changed':False,'source_data_sha256_unchanged':before,'figures':records}
    (out/'provenance/complete_conditions_20260908.json').write_text(json.dumps(receipt,indent=2)+'\n')
    write_legends(out)
    print('Verified unchanged source measurements and complete A–I coverage; English master with bilingual panels and legends.')
if __name__=='__main__':main()
