#!/usr/bin/env python3
"""Human review of actual microglial cell/soma masks and morphology classes."""
from __future__ import annotations
import argparse
import csv
import fcntl
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import sys
import threading
from datetime import datetime,timezone
from urllib.parse import urlparse
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi
from skimage.morphology import skeletonize

PAPER=next(p for p in Path(__file__).resolve().parents if (p/'Fig5/05_review_microglia_classes.py').is_file())
spec=importlib.util.spec_from_file_location('fig5_class_reviewer',PAPER/'Fig5/05_review_microglia_classes.py')
app=importlib.util.module_from_spec(spec);spec.loader.exec_module(app)
base=app.base
base.OUTPUT_FIELDS=tuple(base.OUTPUT_FIELDS)+('mask_reviewed','mask_revision','soma_mask_sha256','cell_mask_sha256')


def png_bytes(array):
    buffer=io.BytesIO();Image.fromarray(array).save(buffer,format='PNG');return buffer.getvalue()


def atomic_bytes(path,data):
    temporary=path.with_suffix('.tmp');temporary.write_bytes(data);os.replace(temporary,path)


def prepare(analysis,root,previous,count=300,seed=20260907):
    receipt=root/'review_batch.json'
    if receipt.is_file():
        meta=app.prepare(analysis,root,count,seed)
        for row in meta['selected_cells']:
            for kind in ['soma','cell']:
                if base.sha256_file(root/'masks'/f'{row["cell_uid"]}_{kind}.png')!=row[f'{kind}_mask_sha256']:
                    raise ValueError(f'{kind} mask changed outside the reviewer: {row["cell_uid"]}')
        return meta
    if root.exists() and any(root.iterdir()):raise ValueError('Use a fresh review directory.')
    cells=pd.read_csv(analysis/'per_microglia_cell_measurements.csv').set_index('cell_uid')
    with (analysis/base.TEMPLATE_NAME).open() as f:rows=list(csv.DictReader(f))
    by_uid={r['cell_uid']:r for r in rows}
    ordering=base.Review.__new__(base.Review);ordering.rows=rows
    order=ordering._review_order(seed)
    previous_rows=[];matches={};mapping=[]
    old_file=previous/base.DEFAULT_OUTPUT_NAME
    if old_file.is_file():
        with old_file.open() as f:previous_rows=list(csv.DictReader(f))
        old_meta=json.loads((previous/'review_batch.json').read_text())
        old_cells=pd.read_csv(Path(old_meta['source_analysis'])/'per_microglia_cell_measurements.csv').set_index('cell_uid')
        used=set()
        for row in sorted(previous_rows,key=lambda r:r['reviewed_at']):
            old=old_cells.loc[row['cell_uid']]
            candidates=cells[(cells['sample']==old['sample'])&(cells['section']==old['section'])]
            distances=np.hypot(candidates.soma_centroid_y-old.centroid_y,candidates.soma_centroid_x-old.centroid_x)
            uid=str(distances.idxmin()) if len(distances) else ''
            distance=float(distances.min()) if len(distances) else None
            mapped=bool(uid and distance<=10 and uid not in used)
            mapping.append(dict(previous_cell_uid=row['cell_uid'],previous_state=row['state'],new_cell_uid=uid if mapped else '',distance_px=distance,
                                status='requires_human_reconfirmation' if mapped else 'no_unambiguous_nearby_match'))
            if mapped:matches[uid]=row;used.add(uid)
    selected=[];tallies={s:0 for s in base.STATES};quota=count//4
    for uid in list(matches)+order:
        if uid in selected:continue
        state=by_uid[uid]['pseudo_state']
        if tallies[state]>=quota:continue
        selected.append(uid);tallies[state]+=1
        if len(selected)==count:break
    if len(selected)<count:
        selected.extend(u for u in order if u not in selected)
        selected=selected[:count]
    if len(selected)!=count:raise ValueError('Not enough cells for the requested review batch.')
    root.mkdir(parents=True)
    for directory in ['patches','masks','mask_history','mask_edits']:(root/directory).mkdir()
    audit=[];template=[]
    for uid in selected:
        row=dict(by_uid[uid]);measurement=cells.loc[uid]
        original=app.source_patch(analysis,row);target=root/'patches'/f'{uid}.png';shutil.copy2(original,target)
        patch=np.array(Image.open(target).convert('RGB'))
        for kind in ['soma','nucleus']:
            shutil.copy2(Path(measurement[f'{kind}_patch_path']),root/'masks'/f'{uid}_{kind}.png')
        cell_path=root/'masks'/f'{uid}_cell.png';cell_path.write_bytes(png_bytes((patch[...,1]>127).astype(np.uint8)*255))
        row['patch_path']=str(target);row['state']='';template.append(row)
        audit.append(dict(cell_uid=uid,patch_sha256=base.sha256_file(target),original_patch_sha256=base.sha256_file(original),source_patch=str(original),
                          soma_mask_sha256=base.sha256_file(root/'masks'/f'{uid}_soma.png'),cell_mask_sha256=base.sha256_file(cell_path),
                          patch_origin_y=int(measurement.patch_origin_y),patch_origin_x=int(measurement.patch_origin_x),
                          source_height=int(measurement.source_height),source_width=int(measurement.source_width),mask_revision=0,
                          previous_state=matches.get(uid,{}).get('state',''),previous_cell_uid=matches.get(uid,{}).get('cell_uid','')))
    with (root/base.TEMPLATE_NAME).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(template[0]));writer.writeheader();writer.writerows(template)
    if old_file.is_file():shutil.copy2(old_file,root/'previous_review_labels_preserved.csv')
    app.atomic_json(root/'previous_label_mapping.json',{'previous_labels':len(previous_rows),'mapping':mapping})
    meta={'source_analysis':str(analysis),'source_template_sha256':base.sha256_file(analysis/base.TEMPLATE_NAME),
          'source_measurements_sha256':base.sha256_file(analysis/'per_microglia_cell_measurements.csv'),
          'cells':count,'seed':seed,'sampling':'Balanced proposed classes, with spatially matched previous review targets prioritised for reconfirmation.',
          'proposed_class_counts':dict(app.collections.Counter(r['pseudo_state'] for r in template)),
          'animal_counts':dict(app.collections.Counter(r['animal_id'] for r in template)),
          'region_counts':dict(app.collections.Counter(r['region'] for r in template)),
          'previous_labels_preserved':len(previous_rows),'previous_targets_in_batch':sum(bool(r['previous_state']) for r in audit),
          'mask_review':True,'selected_cells':audit}
    app.atomic_json(receipt,meta)
    return meta


class Review(base.Review):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.meta=json.loads((self.analysis_dir/'review_batch.json').read_text())
        self.audit={r['cell_uid']:r for r in self.meta['selected_cells']}
        prior=[r['cell_uid'] for r in self.meta['selected_cells'] if r.get('previous_state')]
        self.order=prior+[u for u in self.order if u not in prior]

    def counts(self):
        with self.lock:
            return super().counts()

    def mask(self,uid,kind):
        return np.array(Image.open(self.analysis_dir/'masks'/f'{uid}_{kind}.png'))>127

    def cell(self,uid):
        with self.lock:
            result=super().cell(uid);row=self.audit[uid]
            result.update(cell_pixels=np.flatnonzero(self.mask(uid,'cell')).tolist(),soma_pixels=np.flatnonzero(self.mask(uid,'soma')).tolist(),
                          nucleus_pixels=np.flatnonzero(self.mask(uid,'nucleus')).tolist(),patch_sha256=row['patch_sha256'],
                          mask_revision=row.get('mask_revision',0),previous_state=row.get('previous_state',''))
            return result

    def apply_mask(self,payload):
        uid=str(payload.get('cell_uid',''));reviewer=str(payload.get('reviewer','')).strip();session=str(payload.get('session','')).strip()
        if uid not in self.audit or not reviewer or not session or len(reviewer)>200 or len(session)>200:raise ValueError('A known cell and reviewer name/session are required.')
        row=self.audit[uid]
        if payload.get('expected_patch_sha256')!=row['patch_sha256'] or payload.get('expected_mask_revision',0)!=row.get('mask_revision',0):raise ValueError('This mask changed; reload the cell before saving.')
        if 'cell_pixels' not in payload:return False
        patch_path=Path(self.by_uid[uid]['patch_path']);patch=np.array(Image.open(patch_path).convert('RGB'))
        shape=patch.shape[:2];size=shape[0]*shape[1]
        masks=[]
        for kind in ['cell','soma']:
            indices=payload.get(f'{kind}_pixels')
            if not isinstance(indices,list) or not indices or len(indices)>size:raise ValueError(f'{kind} mask must be nonempty.')
            if any(not isinstance(i,int) or isinstance(i,bool) or i<0 or i>=size for i in indices):raise ValueError('Invalid mask pixels.')
            mask=np.zeros(size,bool);mask[indices]=True;masks.append(mask.reshape(shape))
        cell,soma=masks
        if np.any(soma&~cell):raise ValueError('The soma must be inside the cell mask.')
        yy,xx=np.indices(shape);inside=(yy+row['patch_origin_y']>=0)&(xx+row['patch_origin_x']>=0)&(yy+row['patch_origin_y']<row['source_height'])&(xx+row['patch_origin_x']<row['source_width'])
        if np.any(cell&~inside):raise ValueError('Mask extends outside the source image.')
        if np.array_equal(cell,self.mask(uid,'cell')) and np.array_equal(soma,self.mask(uid,'soma')):return False
        distance=ndi.distance_transform_edt(cell);distance/=max(float(distance.max()),1e-6)
        patch[...,1]=cell.astype(np.uint8)*255
        patch[...,2]=np.clip((.65*skeletonize(cell)+.35*distance)*255,0,255).astype(np.uint8)
        revision=row.get('mask_revision',0)+1
        history=self.analysis_dir/'mask_history'/f'{uid}_r{revision:03d}';history.mkdir()
        destinations={'patch':patch_path,'cell':self.analysis_dir/'masks'/f'{uid}_cell.png','soma':self.analysis_dir/'masks'/f'{uid}_soma.png'}
        old_bytes={k:p.read_bytes() for k,p in destinations.items()}
        for kind,data in old_bytes.items():(history/f'before_{kind}.png').write_bytes(data)
        prior=dict(row);prior_label=self.labels.get(uid)
        try:
            atomic_bytes(patch_path,png_bytes(patch));atomic_bytes(destinations['cell'],png_bytes(cell.astype(np.uint8)*255));atomic_bytes(destinations['soma'],png_bytes(soma.astype(np.uint8)*255))
            row.update(patch_sha256=base.sha256_file(patch_path),soma_mask_sha256=base.sha256_file(destinations['soma']),cell_mask_sha256=base.sha256_file(destinations['cell']),mask_revision=revision)
            self.labels.pop(uid,None);self._write()
            app.atomic_json(self.analysis_dir/'review_batch.json',self.meta)
            np.savez_compressed(self.analysis_dir/'mask_edits'/f'{uid}.npz',cell_mask=cell,soma_mask=soma,
                                patch_origin_y=row['patch_origin_y'],patch_origin_x=row['patch_origin_x'],sample=self.by_uid[uid]['sample'],section=self.by_uid[uid]['section'])
            app.atomic_json(history/'edit.json',{'cell_uid':uid,'reviewer':reviewer,'session':session,'edited_at':datetime.now(timezone.utc).isoformat(),'before':prior,'after':dict(row),'prior_classification':prior_label})
        except Exception:
            for kind,path in destinations.items():atomic_bytes(path,old_bytes[kind])
            row.clear();row.update(prior)
            if prior_label:self.labels[uid]=prior_label
            self._write();app.atomic_json(self.analysis_dir/'review_batch.json',self.meta)
            raise
        return True

    def save(self,payload):
        if str(payload.get('state','')) not in base.STATES:raise ValueError('Choose a valid morphology class.')
        with self.lock:self.apply_mask(payload)
        result=super().save(payload)
        with self.lock:
            uid=payload['cell_uid'];row=self.audit[uid]
            self.labels[uid].update(mask_reviewed=True,mask_revision=row.get('mask_revision',0),soma_mask_sha256=row['soma_mask_sha256'],cell_mask_sha256=row['cell_mask_sha256'])
            self._write()
        return result


class Handler(app.Handler):
    def end_headers(self):
        self.send_header('Cache-Control','no-store');super().end_headers()

    def do_GET(self):
        if urlparse(self.path).path=='/':self.send_bytes(HTML.encode(),'text/html; charset=utf-8')
        else:super().do_GET()

    def _post(self):
        if urlparse(self.path).path!='/api/mask':return super()._post()
        if self.server.training.state['status']!='waiting_for_labels':return self.error_json(409,'Review is finished; masks are frozen for training.')
        try:
            if self.headers.get_content_type()!='application/json':raise ValueError('Content-Type must be application/json')
            size=int(self.headers.get('Content-Length','0'))
            if not 0<size<=base.MAX_BODY:raise ValueError('Invalid mask request size.')
            payload=json.loads(self.rfile.read(size));uid=payload.get('cell_uid','')
            with self.review.lock:self.review.apply_mask(payload)
            self.send_json({'ok':True,'cell':self.review.cell(uid),**self.review.counts()})
        except (ValueError,KeyError,TypeError) as exc:self.error_json(422,str(exc))


MASK_RENDER = r'''function maskPayload(){const p={expected_patch_sha256:cur.patch_sha256,expected_mask_revision:cur.mask_revision||0};if(dirty){p.cell_pixels=[];p.soma_pixels=[];for(let i=0;i<cellMask.length;i++){if(cellMask[i])p.cell_pixels.push(i);if(somaMask[i])p.soma_pixels.push(i)}}return p}
function edge(mask,x,y,w,h){const i=y*w+x;return mask[i]&&(x===0||y===0||x===w-1||y===h-1||!mask[i-1]||!mask[i+1]||!mask[i-w]||!mask[i+w])}
function drawMasks(){if(!imageReady)return;const w=$('overlay').width,h=$('overlay').height;for(const name of ['signal','overlay']){const ctx=$(name).getContext('2d'),out=ctx.createImageData(w,h);for(let y=0;y<h;y++)for(let x=0;x<w;x++){const i=y*w+x,j=i*4;let v=signalPixels[i],rgb=[v,v,v];if(name==='overlay'){if(!cellMask[i])rgb=rgb.map(v=>v*.65);if($('showprocess').checked&&cellMask[i])rgb=rgb.map((v,k)=>v*.85+[45,220,240][k]*.15);if($('showprocess').checked&&edge(cellMask,x,y,w,h))rgb=[45,225,240];if($('showbody').checked&&somaMask[i])rgb=rgb.map((v,k)=>v*.6+[255,185,40][k]*.4);if($('showbody').checked&&edge(somaMask,x,y,w,h))rgb=[255,185,40];if($('shownucleus').checked&&edge(nucleusMask,x,y,w,h))rgb=[220,100,255]}out.data[j]=rgb[0];out.data[j+1]=rgb[1];out.data[j+2]=rgb[2];out.data[j+3]=255}ctx.putImageData(out,0,0)}$('editstatus').textContent=dirty?'Mask changes pending — save the mask or choose a class.':`Mask revision ${cur.mask_revision||0}`;controls()}
function render(c){cur=c;dirty=false;imageReady=false;painting=false;lastPoint=null;controls();const im=new Image();im.onload=()=>{const tmp=document.createElement('canvas');tmp.width=im.width;tmp.height=im.height;const ctx=tmp.getContext('2d');ctx.drawImage(im,0,0);const data=ctx.getImageData(0,0,im.width,im.height).data;signalPixels=new Uint8Array(im.width*im.height);cellMask=new Uint8Array(signalPixels.length);somaMask=new Uint8Array(signalPixels.length);nucleusMask=new Uint8Array(signalPixels.length);for(let i=0;i<signalPixels.length;i++)signalPixels[i]=data[i*4];for(const i of c.cell_pixels||[])cellMask[i]=1;for(const i of c.soma_pixels||[])somaMask[i]=1;for(const i of c.nucleus_pixels||[])nucleusMask[i]=1;for(const n of ['signal','overlay']){$(n).width=im.width;$(n).height=im.height}imageReady=true;drawMasks();controls()};im.onerror=()=>{status('Image could not load. Reload before reviewing.');imageReady=false;controls()};im.src='/api/patch?cell_uid='+encodeURIComponent(c.cell_uid)+'&revision='+encodeURIComponent(c.patch_sha256);$('meta').textContent=`${c.cell_uid} · ${c.region}`+(c.previous_state?` · Previous label: ${c.previous_state}; reassess this revised mask.`:'');document.querySelectorAll('.state').forEach(b=>b.classList.toggle('selected',b.dataset.state===c.state));proposal();status(c.state?`Saved as ${c.state}. Editing the mask requires confirming the class again.`:'')}
function brushPoint(x,y){const w=$('overlay').width,h=$('overlay').height,r=Number($('brushsize').value),mode=$('masktool').value;for(let yy=Math.max(0,y-r);yy<=Math.min(h-1,y+r);yy++)for(let xx=Math.max(0,x-r);xx<=Math.min(w-1,x+r);xx++){if((xx-x)**2+(yy-y)**2>r*r)continue;const i=yy*w+xx;if(mode==='soma-add'){somaMask[i]=1;cellMask[i]=1}else if(mode==='soma-erase')somaMask[i]=0;else if(mode==='cell-add')cellMask[i]=1;else if(mode==='cell-erase'){cellMask[i]=0;somaMask[i]=0}}dirty=true}
function pointerPosition(e){const r=$('overlay').getBoundingClientRect();return [Math.max(0,Math.min($('overlay').width-1,Math.floor((e.clientX-r.left)/r.width*$('overlay').width))),Math.max(0,Math.min($('overlay').height-1,Math.floor((e.clientY-r.top)/r.height*$('overlay').height)))]}
$('overlay').addEventListener('pointerdown',e=>{if(busy||trainingActive||!imageReady||$('masktool').value==='inspect')return;e.preventDefault();painting=true;lastPoint=pointerPosition(e);$('overlay').setPointerCapture(e.pointerId);brushPoint(...lastPoint);drawMasks()});
$('overlay').addEventListener('pointermove',e=>{if(!painting)return;const p=pointerPosition(e),steps=Math.max(1,Math.ceil(Math.hypot(p[0]-lastPoint[0],p[1]-lastPoint[1])));for(let j=1;j<=steps;j++)brushPoint(Math.round(lastPoint[0]+(p[0]-lastPoint[0])*j/steps),Math.round(lastPoint[1]+(p[1]-lastPoint[1])*j/steps));lastPoint=p;drawMasks()});
for(const name of ['pointerup','pointercancel'])$('overlay').addEventListener(name,()=>painting=false);
for(const name of ['showbody','showprocess','shownucleus'])$(name).onchange=drawMasks;
$('resetmask').onclick=()=>{if(cur&&!busy)render(cur)};
$('savemask').onclick=async()=>{if(busy||trainingActive||!cur||!dirty)return;const reviewer=$('reviewer').value.trim();if(!reviewer){status('Enter your name before saving a mask.');$('reviewer').focus();return}busy=true;controls();try{const d=await api('/api/mask',{cell_uid:cur.cell_uid,reviewer,session,...maskPayload()});stats(d);render(d.cell);status('Corrected cell and soma masks saved. Choose the morphology class when ready.')}catch(e){status(e.message)}finally{busy=false;controls()}};
'''

HTML=app.HTML.replace('Figure 5 · Microglia review','Figure 5 · Microglia masks & classes')
HTML=HTML.replace('Same signal · target cell outlined in cyan','Cell/process mask: cyan · soma mask: amber')
HTML=HTML.replace('<div class="help">Classify the outlined cell. Each choice saves immediately and opens the next cell.</div>',
'''<div class="help">Check the target mask, correct it if needed, then choose a class. Your choice saves both the mask and classification.</div>
<div class="masktools"><label>Mask tool <select id="masktool"><option value="inspect">Inspect</option><option value="soma-add">Paint soma</option><option value="soma-erase">Erase soma</option><option value="cell-add">Paint cell/processes</option><option value="cell-erase">Erase cell/processes</option></select></label><label>Brush <input id="brushsize" type="range" min="1" max="5" value="1"></label><button id="savemask">Save mask</button><button id="resetmask">Undo unsaved mask edits</button></div>
<div class="maskviews"><label><input type="checkbox" id="showbody" checked> Soma</label><label><input type="checkbox" id="showprocess" checked> Cell/processes</label><label><input type="checkbox" id="shownucleus"> DAPI nucleus</label><span id="editstatus" class="muted"></span></div>''')
HTML=HTML.replace('</style>','''#overlay{touch-action:none;cursor:crosshair}.masktools,.maskviews{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin:12px 0}.masktools select{padding:6px;background:#22333f;color:var(--ink);border:1px solid var(--line);border-radius:4px}.masktools button{font-size:12px;padding:6px 9px}.masktools input[type=range]{width:75px}.maskviews{font-size:12px}.maskviews label{white-space:nowrap}select{font:inherit}</style>''')
HTML=HTML.replace("imageReady=false,trainingActive=false;","imageReady=false,trainingActive=false,trainingReady=false,dirty=false,signalPixels=null,cellMask=null,somaMask=null,nucleusMask=null,painting=false,lastPoint=null;")
start=HTML.index('function render(c){');end=HTML.index('\nfunction stats(d){',start)
HTML=HTML[:start]+MASK_RENDER+HTML[end:]
HTML=HTML.replace('cell_uid:uid,state:s,reviewer,session});','cell_uid:uid,state:s,reviewer,session,...maskPayload()});')
assert HTML.count('...maskPayload()') == 2
HTML=HTML.replace("if(busy||!cur)return;busy=true;controls();try{await next(cur.cell_uid)}", "if(busy||!cur)return;if(dirty){status('Save or undo mask edits before skipping.');return}busy=true;controls();try{await next(cur.cell_uid)}")
HTML=HTML.replace("if(busy||!history.length)return;busy=true;", "if(busy||!history.length)return;if(dirty){status('Save or undo mask edits before going back.');return}busy=true;")
HTML=HTML.replace('/input|textarea|button/i','/input|textarea|button|select/i')
HTML=HTML.replace("$('train').disabled=busy||!t.ready", "$('train').disabled=busy||dirty||!t.ready")
HTML=HTML.replace("const t=d.training;trainingActive=", "const t=d.training;trainingReady=t.ready;trainingActive=")
HTML=HTML.replace("$('back').disabled=busy||!history.length}", "$('back').disabled=busy||!history.length;$('savemask').disabled=busy||trainingActive||!cur||!imageReady||!dirty;$('resetmask').disabled=busy||trainingActive||!cur||!dirty;$('masktool').disabled=busy||trainingActive||!cur||!imageReady;$('brushsize').disabled=busy||trainingActive||!cur||!imageReady;$('train').disabled=busy||dirty||trainingActive||!trainingReady}")
HTML=HTML.replace("$('train').onclick=async()=>{busy=true;", "$('train').onclick=async()=>{if(busy||trainingActive||!trainingReady)return;if(dirty){status('Save or undo mask edits before finishing review.');return}busy=true;")
HTML=HTML.replace('Every choice saves immediately.', 'Every choice saves immediately.')
HTML=HTML.replace('Labels are saved locally after every choice. You can close this page and resume later.',
                  'Masks and labels are saved locally. Earlier classifications are preserved; matched cells are shown first for reassessment with their revised masks.')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--analysis-dir',type=Path,default=PAPER/'analyses/Fig5/results/soma_refined_20260907_v1')
    ap.add_argument('--review-dir',type=Path,default=PAPER/'Fig5/hil_review/microglia_masks_20260907_v2')
    ap.add_argument('--previous-review',type=Path,default=PAPER/'Fig5/hil_review/microglia_classes_20260906_v1')
    ap.add_argument('--port',type=int,default=8765)
    ap.add_argument('--count',type=int,default=300)
    ap.add_argument('--prepare-only',action='store_true')
    ap.add_argument('--training-python',default=sys.executable)
    ap.add_argument('--exit-after-training',action='store_true',help='Return to the reproduction workflow after training completes or fails.')
    args=ap.parse_args();root=args.review_dir.resolve()
    meta=prepare(args.analysis_dir.resolve(),root,args.previous_review.resolve(),args.count)
    if args.prepare_only:
        print(json.dumps({k:v for k,v in meta.items() if k!='selected_cells'},indent=2));return
    lock=(root/'server.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    review=Review(root,root/base.DEFAULT_OUTPUT_NAME,20260907)
    server=base.ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    server.daemon_threads=True;server.review=review;server.mutation_lock=threading.Lock()
    server.training=app.Training(review,args.training_python,200)
    app.atomic_json(root/'server.json',{'pid':os.getpid(),'url':f'http://127.0.0.1:{server.server_port}','port':server.server_port,'review_dir':str(root),'cells':len(review.rows)})
    print(f'Mask and class reviewer: http://127.0.0.1:{server.server_port}; {len(review.rows)} cells; {meta["previous_labels_preserved"]} earlier labels preserved.',flush=True)
    if review.counts()['labelled']==len(review.rows) and server.training.info()['ready'] and server.training.state['status']=='waiting_for_labels':server.training.start()
    print(f'MICROGLIAL_REVIEW_URL: http://127.0.0.1:{server.server_port}',flush=True)
    finished=threading.Event()
    def await_training():
        while not finished.wait(.5):
            if server.training.state['status'] in {'complete','failed'}:
                server.shutdown()
                return
    if args.exit_after_training:
        threading.Thread(target=await_training,daemon=True).start()
    try:server.serve_forever(.25)
    except KeyboardInterrupt:
        if args.exit_after_training:raise
    finally:finished.set();server.server_close();lock.close()
    if args.exit_after_training and server.training.state['status']!='complete':
        raise RuntimeError('Microglial retraining did not complete; choices are saved in the review directory.')


if __name__=='__main__':main()
