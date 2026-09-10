"""Render S7/S8 from verified source crops, saved losses and pilot power curves."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.transforms import Bbox
from pomc_size_qc import filter_pomc_by_local_dapi

LEGENDS={
'S7':{
'en':'Figure S7. Human-supervised Cellpose training. A, Initial segmentation, human correction and model retraining. B, Source 512×512-pixel DAPI, c-FOS and POMC training crops (upper row) and the same crops with reviewed mask boundaries in green (lower row). POMC overlays include a post-segmentation size filter: objects smaller than the mean DAPI nuclear area in the same field and anatomical region were excluded. This revision did not retrain the historical models represented by the loss curves. Display contrast uses the 1st and 99.7th percentiles, identically within each signal/overlay pair. C, Saved training losses for five checkpoints, with one value per epoch and a logarithmic vertical axis. These curves describe fit to training annotations, not independent test performance. The source crop coordinates and original TIFF/mask hashes are preserved in provenance/training_crop_sources.json; the complete training inventory is in Supplementary Table S6.',
'es':'Figura S7. Entrenamiento de Cellpose con supervisión humana. A, Segmentación inicial, corrección humana y reentrenamiento del modelo. B, Recortes de 512×512 píxeles de las imágenes de entrenamiento DAPI, c-FOS y POMC (fila superior), y los mismos recortes con bordes de máscaras revisadas en verde (fila inferior). Las máscaras POMC incluyen un filtro posterior a la segmentación: se excluyeron objetos menores que el área media de los núcleos DAPI del mismo campo y región anatómica. Esta revisión no reentrenó los modelos históricos representados por las curvas de pérdida. El contraste de presentación utiliza los percentiles 1 y 99,7 de manera idéntica en cada par de señal y superposición. C, Pérdidas de entrenamiento guardadas para cinco modelos, con un valor por época y eje vertical logarítmico. Estas curvas describen el ajuste a las anotaciones de entrenamiento y no el rendimiento en un conjunto independiente. Las coordenadas y huellas digitales de los TIFF y máscaras originales se conservan en provenance/training_crop_sources.json; el inventario de entrenamiento se entrega en la Tabla suplementaria S6.'},
'S8':{
'en':'Figure S8. Pilot-based Monte Carlo sample-size estimates. A, Day-6 cumulative bottle-volume removal. B, Day-6 percentage body-weight change. Solid curves retain the pilot group-specific standard deviations and use Welch ANOVA; dashed curves use a pooled common variance and ordinary ANOVA. The horizontal line marks 80% power. Points identify the first balanced group size reaching that threshold. Each size used 500,000 normally distributed simulations, alpha=0.05 and PCG64 seed 20260908. The pilot comprised 12 independent cages (water/sucrose/allulose, 3/5/4). Group-specific minima were five cages per group for volume and fourteen for weight; common-variance minima were four and nine. Curves are displayed through 24 cages per group; the full table preserves all simulated sizes. These are conditional planning estimates for the pilot outcomes, not achieved power for neuronal or glial comparisons. Source values, assumptions and curves are in Supplementary Table S6.',
'es':'Figura S8. Estimación del tamaño muestral por Monte Carlo basada en el piloto. A, Volumen acumulado retirado de la botella al día 6. B, Cambio porcentual de peso corporal al día 6. Las curvas continuas conservan las desviaciones estándar propias de cada grupo del piloto y usan ANOVA de Welch; las discontinuas utilizan una varianza común combinada y ANOVA convencional. La línea horizontal señala 80% de potencia. Los puntos identifican el primer tamaño de grupo equilibrado que alcanza ese umbral. Cada tamaño utilizó 500.000 simulaciones normales, alfa=0,05 y semilla PCG64 20260908. El piloto incluyó 12 jaulas independientes (agua/sacarosa/alulosa, 3/5/4). Los mínimos con varianzas propias fueron cinco jaulas por grupo para volumen y catorce para peso; con varianza común fueron cuatro y nueve. Se muestran las curvas hasta 24 jaulas por grupo; la tabla completa conserva todos los tamaños simulados. Son estimaciones de planificación condicionadas al piloto y no potencia alcanzada para las comparaciones neuronales o gliales. Los valores, supuestos y curvas se entregan en la Tabla suplementaria S6.'}}

def verify_inputs(root):
 receipt=json.loads((root/'provenance/input_manifest.json').read_text())
 for row in receipt['files']:
  path=root/row['path']
  if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
   raise ValueError('Source input hash differs: '+str(path))
 return receipt

def training_size_masks(root):
 cleaned={};rows=[]
 for path in sorted((root/'raw_data/pomc_training_size_reference').glob('*.npz')):
  with np.load(path,allow_pickle=False) as z:
   mask,audit=filter_pomc_by_local_dapi(z['mask'],z['dapi'],z['regions'])
  cleaned[path.stem]=mask
  rows.extend({'field':path.stem,**obj} for obj in audit['objects'])
 if len(cleaned)!=10:raise ValueError('Ten POMC training fields with local DAPI references are required')
 return cleaned,rows

def training(root,lang):
 es=lang=='es';examples=[];size_masks,_=training_size_masks(root)
 for channel in ['DAPI','c-FOS','POMC']:
  with np.load(root/'raw_data'/(channel+'.npz'),allow_pickle=False) as data:crop=data['image'].astype(float);m=data['mask']
  if channel=='POMC':m=size_masks['FR6-1_S01_right'][208:720,0:512]
  assert crop.shape==m.shape==(512,512)
  lo,hi=np.percentile(crop,[1,99.7]);crop=np.clip((crop-lo)/max(hi-lo,1e-8),0,1)
  border=(m>0)&((m!=np.roll(m,1,0))|(m!=np.roll(m,-1,0))|(m!=np.roll(m,1,1))|(m!=np.roll(m,-1,1)))
  border=border|np.roll(border,1,0)|np.roll(border,1,1)
  rgb=np.repeat(crop[:,:,None],3,axis=2);rgb[border]=[.15,1,.28];examples.append((channel,crop,rgb))
 fig=plt.figure(figsize=(8,8.15));gs=fig.add_gridspec(4,3,height_ratios=[.5,1,1,1.18],hspace=.34,wspace=.1)
 ax=fig.add_subplot(gs[0,:]);ax.axis('off')
 labels=['TIFF original','Segmentación inicial','Corrección humana','Reentrenamiento'] if es else ['Source TIFF','Initial segmentation','Human correction','Model retraining']
 for i,label in enumerate(labels):
  x=.115+i*.255;ax.text(x,.6,label,ha='center',va='center',fontsize=9,bbox=dict(boxstyle='round,pad=.55',fc='#eaf1f5',ec='#7495a7'))
  if i<3:ax.add_patch(FancyArrowPatch((x+.091,.6),(x+.159,.6),arrowstyle='-|>',mutation_scale=12,color='#466579'))
 for j,(channel,crop,overlay) in enumerate(examples):
  ax=fig.add_subplot(gs[1,j]);ax.imshow(crop,cmap='gray',vmin=0,vmax=1);ax.set_title(channel,fontsize=11);ax.axis('off')
  ax=fig.add_subplot(gs[2,j]);ax.imshow(overlay);ax.axis('off')
 fig.text(.023,.681,'Imagen' if es else 'Image',rotation=90,va='center',fontsize=10)
 fig.text(.023,.459,'Máscara tras control' if es else 'Mask after QC',rotation=90,va='center',fontsize=10)
 ax=fig.add_subplot(gs[3,:]);losses=pd.read_csv(root/'source_data/training_losses.csv')
 for name,color in zip(['dapi_channel','dapi_Keyence4','cfos_channel3','cfos_Keyence1','POMC_2'],['#2166ac','#69a9d2','#b23a97','#e493c4','#c58915']):
  row=losses[losses.model==name];ax.plot(row.epoch,row.training_loss,lw=.9,color=color,label=name,alpha=.9)
 ax.set_yscale('log');ax.set_xlabel('Época de entrenamiento' if es else 'Training epoch',fontsize=10);ax.set_ylabel('Pérdida de entrenamiento\n(escala logarítmica)' if es else 'Training loss\n(log scale)',fontsize=9)
 ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12);ax.legend(ncol=3,fontsize=8,frameon=False,loc='upper right')
 fig.subplots_adjust(left=.09,right=.985,top=.975,bottom=.065)
 for letter,y in [('A',.985),('B',.855),('C',.329)]:fig.text(.009,y,letter,fontsize=13,fontweight='bold',va='top')
 boxes={'A':Bbox.from_extents(0,7.15,8,8.15),'B':Bbox.from_extents(0,2.80,8,7.12),'C':Bbox.from_extents(0,0,8,2.78)}
 return fig,boxes

def power(root,lang):
 data=pd.read_csv(root/'source_data/power_curves.csv');es=lang=='es';fig,axs=plt.subplots(1,2,figsize=(8,3.45))
 titles=['Volumen retirado al día 6','Cambio de peso al día 6'] if es else ['Volume removed at day 6','Body-weight change at day 6']
 for letter,ax,endpoint,title in zip('AB',axs,['bottle_volume_removed_day6','body_weight_change_day6'],titles):
  for scenario,color,style,label in [('common_pooled_variance','#9b5918','--','ANOVA: varianza común' if es else 'ANOVA: common variance'),('group_specific_variance','#185f8d','-','Welch: varianzas por grupo' if es else 'Welch: group-specific variances')]:
   a=data[(data.endpoint==endpoint)&(data.scenario==scenario)].sort_values('n_cages_per_group');ax.plot(a.n_cages_per_group,a.power*100,style,color=color,lw=1.7,label=label)
   first=a[a.power>=.8].iloc[0];ax.scatter(first.n_cages_per_group,first.power*100,color=color,s=25,zorder=5)
   ax.annotate(f'n = {int(first.n_cages_per_group)}',(first.n_cages_per_group,first.power*100),xytext=(4,10 if scenario.startswith('common') else -16),textcoords='offset points',fontsize=8,color=color)
  ax.axhline(80,color='#555555',ls=':',lw=1);ax.set(title=title,ylim=(0,104),xlabel='Jaulas independientes por grupo' if es else 'Independent cages per group');ax.set_ylabel('Potencia estimada (%)' if es else 'Estimated power (%)');ax.set_xlim(2,24);ax.set_xticks([2,5,10,15,20,24]);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.legend(loc='lower right',fontsize=6.8,frameon=False)
  ax.text(-.19,1.14,letter,transform=ax.transAxes,fontweight='bold',fontsize=13,va='top')
 fig.tight_layout(pad=1.2);fig.canvas.draw();renderer=fig.canvas.get_renderer()
 boxes={letter:ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted()).expanded(1.025,1.04) for letter,ax in zip('AB',axs)}
 return fig,boxes

def main(number,default_root):
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input-root',type=Path,default=default_root);parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--dpi',type=int,default=600);args=parser.parse_args()
 inputs=verify_inputs(args.input_root);out=args.output_dir
 if out.exists():raise FileExistsError('Use a fresh output directory: '+str(out))
 for folder in ['panels','legends','provenance','source_data']:(out/folder).mkdir(parents=True,exist_ok=True)
 if number=='S7':
  _,qc_rows=training_size_masks(args.input_root)
  pd.DataFrame(qc_rows).to_csv(out/'source_data/pomc_training_size_qc.csv',index=False)
 for lang in ['en','es']:
  fig,boxes=(training if number=='S7' else power)(args.input_root,lang);suffix='' if lang=='en' else '_spanish'
  if lang=='en':
   for ext in ['.png','.pdf']:fig.savefig(out/('Figure_'+number+ext),dpi=args.dpi)
  for letter,box in boxes.items():
   for ext in ['.png','.pdf']:fig.savefig(out/'panels'/f'Figure_{number}_Panel_{letter}{suffix}{ext}',dpi=args.dpi,bbox_inches=box)
   (out/'legends'/f'Figure_{number}_Panel_{letter}_LEGEND{suffix}.txt').write_text(LEGENDS[number][lang]+'\n'+('Isolated panel ' if lang=='en' else 'Panel aislado ')+letter+'.\n')
  (out/'legends'/f'Figure_{number}_LEGEND{suffix}.txt').write_text(LEGENDS[number][lang]+'\n');plt.close(fig)
 receipt={'figure':number,'status':'complete','inputs':inputs,'dpi':args.dpi,'master_language':'English','panel_languages':['English','Spanish'],'raw_images_or_quantification_modified':False,'outputs':[]}
 for p in sorted(out.rglob('*')):
  if p.is_file():receipt['outputs'].append({'path':str(p.relative_to(out)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 (out/'provenance/render_receipt.json').write_text(json.dumps(receipt,indent=2));print('Rendered',number,'from verified source inputs:',out)
