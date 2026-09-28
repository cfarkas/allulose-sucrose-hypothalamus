#!/usr/bin/env python3
"""Review the sucrose ventricular boundary on the registered native DAPI image.

The reviewer draws and saves the polygon in the browser. It is a display-only
annotation: no nucleus classification, spatial analysis, or statistic changes.
"""
from pathlib import Path
import argparse, ast, hashlib, importlib.util, io, json, shutil, sys, threading
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs
from uuid import uuid4
import cv2
import numpy as np
from PIL import Image

DEFAULT_REVIEW = Path('Fig3/hil_review/sucrose_ventricle_20260908_v1')
SCHEMA = 'fig3_sucrose_native_dapi_ventricle_hil_v1'

def module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return m

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def png(a):
    b=io.BytesIO();Image.fromarray(a).save(b,format='PNG');return b.getvalue()
def atomic(path,payload):
    p=path.with_name(path.name+'.'+uuid4().hex+'.tmp');p.write_bytes(payload);p.replace(path)

def prepare(root, output):
    m=module('sucrose_native_figure3',root/'Fig3/03_make_figure_3_cfos_npy.py')
    spec=next(s for s in m.REPRESENTATIVE_FIELDS if s['condition']=='Sucrose')
    base=root/'Fig3/raw/legacy_experiment_2025_07_28'
    labels={k:m._load_cellpose_labels(base/v,v)[0] for k,v in spec['mask_paths'].items()}
    *_,counts=m._analysis_positive_sets(labels['dapi'],labels['cfos'],labels['npy'])
    for key in ('dapi_nuclei','cfos_positive_nuclei','npy_positive_nuclei','dual_positive_nuclei'):
        if counts[key]!=spec['expected_analysis_counts'][key]:raise ValueError((key,counts))
    displayed=m._warp_labels_to_display(labels['dapi'],spec)
    raw_path=base/'E7-FR7-5_dapi.tif'
    raw=np.array(Image.open(raw_path))
    gray=raw.max(axis=2) if raw.ndim==3 else raw
    height,width=displayed.shape
    registered=cv2.warpPerspective(gray.astype(np.float32),np.asarray(spec['display_to_raw'],np.float64),(width,height),flags=cv2.INTER_LINEAR|cv2.WARP_INVERSE_MAP,borderMode=cv2.BORDER_CONSTANT,borderValue=0)
    positive=registered[registered>0]
    low,high=np.percentile(positive,[.5,99.7])
    norm=np.clip((registered-low)/max(high-low,1e-8),0,1)**.55
    background=np.stack((norm*.72,norm*.82,norm),axis=-1)
    background=np.rint(background*255).astype(np.uint8)
    overlay=np.zeros((height,width,4),np.uint8)
    edges=cv2.morphologyEx((displayed>0).astype(np.uint8),cv2.MORPH_GRADIENT,np.ones((3,3),np.uint8))>0
    overlay[edges]=[0,255,180,150]
    output.mkdir(parents=True,exist_ok=True)
    inputs=output/'review_inputs';inputs.mkdir(exist_ok=True)
    files={str(raw_path.relative_to(root)):sha(raw_path)}
    files.update({str((base/v).relative_to(root)):sha(base/v) for v in spec['mask_paths'].values()})
    bg=png(background);ov=png(overlay)
    manifest={'schema':SCHEMA,'condition':'Sucrose','animal_id':spec['sample'],'coordinate_frame_yx':[height,width],'display_to_native_xy':spec['display_to_raw'],'pixel_size_um':spec['um_per_px'],'source_sha256':files,'native_counts':counts,'registered_dapi_labels_sha256':hashlib.sha256(displayed.astype('<i4').tobytes()).hexdigest(),'background_sha256':hashlib.sha256(bg).hexdigest(),'background_display':{'percentiles':[.5,99.7],'limits':[float(low),float(high)],'gamma':.55},'automatic_ventricle_proposal':False}
    mp=inputs/'REVIEW_INPUT_MANIFEST.json'
    if mp.exists() and json.loads(mp.read_text())!=json.loads(json.dumps(manifest)):
        raise ValueError('Source/geometry changed: use a new review directory.')
    mp.write_text(json.dumps(manifest,indent=2)+'\n')
    (inputs/'registered_native_DAPI.png').write_bytes(bg)
    (inputs/'segmented_nuclei_overlay.png').write_bytes(ov)
    return manifest,bg,ov

def load_accepted(root, output, expected_shape=None):
    """Validate the saved human annotation, source identities and frame."""
    path=output/'sucrose_ventricle_receipt.json';mask_path=output/'sucrose_ventricle_mask.png'
    if not path.exists():return None
    rec=json.loads(path.read_text());manifest=json.loads((output/'review_inputs/REVIEW_INPUT_MANIFEST.json').read_text())
    if rec.get('schema')!=SCHEMA or rec.get('status')!='accepted':raise ValueError('Invalid HIL receipt')
    if rec['review_manifest_sha256']!=sha(output/'review_inputs/REVIEW_INPUT_MANIFEST.json'):raise ValueError('Review inputs changed')
    if rec['mask_file_sha256']!=sha(mask_path):raise ValueError('Accepted mask changed')
    for name,digest in manifest['source_sha256'].items():
        if sha(root/name)!=digest:raise ValueError('Native input changed: '+name)
    mask=np.array(Image.open(mask_path))
    if list(mask.shape)!=manifest['coordinate_frame_yx']:raise ValueError('Wrong coordinate frame')
    if expected_shape is not None and tuple(mask.shape)!=tuple(expected_shape):raise ValueError('Cartoon frame differs from HIL')
    if not rec.get('reviewer') or not rec.get('session'):raise ValueError('Missing human review identity')
    vertices=np.array(rec['polygon_vertices_xy'],float)
    check=np.zeros(mask.shape,np.uint8);cv2.fillPoly(check,[np.rint(vertices).astype(np.int32)],255)
    if not np.array_equal(mask,check):raise ValueError('Saved polygon and mask differ')
    return mask>0,rec

class State:
    def __init__(self,root,output):
        self.root,self.output=root,output
        self.manifest,self.background,self.overlay=prepare(root,output)
        self.lock=threading.Lock()
    def annotation(self):
        accepted=load_accepted(self.root,self.output)
        return {'found':False,'condition':'Sucrose'} if accepted is None else {'found':True,**accepted[1]}
    def config(self):
        found=self.annotation()['found'];h,w=self.manifest['coordinate_frame_yx']
        return {'sections':[{'condition':'Sucrose','animal_id':self.manifest['animal_id'],'width':w,'height':h,'accepted':found}],'accepted':int(found),'total':1,'output_dir':str(self.output),'figure_button_enabled':False}
    def accept(self,p):
        with self.lock:
            if p.get('condition')!='Sucrose':raise ValueError('Only the Sucrose section is reviewed here')
            reviewer=str(p.get('reviewer','')).strip();session=str(p.get('session','')).strip();notes=str(p.get('notes','')).strip()
            if not reviewer or not session or max(map(len,(reviewer,session,notes)))>1000:raise ValueError('Reviewer and session are required (maximum 1,000 characters)')
            v=np.asarray(p.get('polygon_vertices_xy',[]),float)
            if v.ndim!=2 or v.shape[1]!=2 or not 3<=len(v)<=2000 or not np.isfinite(v).all():raise ValueError('Draw one polygon with 3–2,000 finite vertices')
            h,w=self.manifest['coordinate_frame_yx'];v=np.round(v,3);rv=np.rint(v).astype(np.int32)
            if np.any(rv<0) or np.any(rv[:,0]>=w) or np.any(rv[:,1]>=h):raise ValueError('Keep all vertices inside the image')
            mask=np.zeros((h,w),np.uint8);cv2.fillPoly(mask,[rv],255)
            contours,hierarchy=cv2.findContours(mask,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
            if np.count_nonzero(mask)<50 or len(contours)!=1:raise ValueError('The polygon must enclose one continuous lumen, without holes')
            payload=png(mask);mp=self.output/'sucrose_ventricle_mask.png';rp=self.output/'sucrose_ventricle_receipt.json'
            if mp.exists() or rp.exists():
                archive=self.output/'history'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid4().hex[:8]);archive.mkdir(parents=True)
                for path in (mp,rp):
                    if path.exists():shutil.copy2(path,archive/path.name)
            rec={'schema':SCHEMA,'status':'accepted','condition':'Sucrose','animal_id':self.manifest['animal_id'],'reviewer':reviewer,'session':session,'notes':notes,'accepted_at_utc':datetime.now(timezone.utc).isoformat(),'coordinate_frame_yx':[h,w],'polygon_vertices_xy':v.tolist(),'mask_file_sha256':hashlib.sha256(payload).hexdigest(),'review_manifest_sha256':sha(self.output/'review_inputs/REVIEW_INPUT_MANIFEST.json'),'semantics':{'display_only':True,'analysis_mask':False,'changes_cell_calls_or_counts':False,'render':'black boundary of the manually delineated ventricle'}}
            atomic(mp,payload);atomic(rp,(json.dumps(rec,indent=2)+'\n').encode())
            load_accepted(self.root,self.output)
            print('HUMAN_POLYGON_ACCEPTED '+str(rp),flush=True)
            return {'status':'accepted','mask_sha256':rec['mask_file_sha256']}

class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args):pass
    def send(self,status,payload,kind='application/json'):
        if isinstance(payload,dict):payload=json.dumps(payload).encode()
        self.send_response(status);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(payload)));self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(payload)
    def do_GET(self):
        try:
            route=urlparse(self.path).path
            if route=='/':return self.send(200,self.server.html.encode(),'text/html; charset=utf-8')
            if route=='/api/config':return self.send(200,self.server.state.config())
            if route=='/api/annotation':return self.send(200,self.server.state.annotation())
            if route=='/api/image':return self.send(200,self.server.state.background,'image/png')
            if route=='/api/source-overlay':return self.send(200,self.server.state.overlay,'image/png')
            if route=='/api/figure-status':return self.send(200,{'status':'idle','running':False,'step':'El contorno guardado se incorporará en la tesis y en el artículo.'})
            self.send(404,{'error':'Unknown route'})
        except Exception as e:self.send(400,{'error':str(e)})
    def do_POST(self):
        try:
            if self.path!='/api/accept':return self.send(404,{'error':'Unknown route'})
            size=int(self.headers.get('Content-Length','0'))
            if size<=0 or size>200000:raise ValueError('Invalid request size')
            self.send(200,self.server.state.accept(json.loads(self.rfile.read(size))))
        except Exception as e:self.send(400,{'error':str(e)})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path.cwd());p.add_argument('--output-dir',type=Path);p.add_argument('--host',default='127.0.0.1');p.add_argument('--port',type=int,default=34248);p.add_argument('--prepare-only',action='store_true');args=p.parse_args()
    root=args.root.resolve();output=args.output_dir or root/DEFAULT_REVIEW
    state=State(root,output)
    if args.prepare_only:print(json.dumps(state.config(),indent=2));return
    ui_tree=ast.parse((root/'Fig3/02a_review_cfos_npy_ventricles.py').read_text())
    ui=next(ast.literal_eval(n.value) for n in ui_tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='HTML' for t in n.targets))
    ui=ui.replace('Figure 3 ventricle HIL','Sacarosa · delineación del ventrículo').replace('Draw one exact closed lumen polygon for Water and Allulose','DAPI original de E7-FR7-5 · dibuja y guarda el contorno ventricular').replace('0 / 2 accepted','0 / 1 accepted').replace("condition:'Water'","condition:'Sucrose'").replace("keep='Water'","keep='Sucrose'")
    ui=ui.replace('Show raw-DAPI proposal','Mostrar núcleos segmentados').replace('type="checkbox" checked','type="checkbox"').replace('The magenta proposal is guidance only. The accepted polygon fully controls the white lumen and black wall in all three cartoon views.','La imagen corresponde al canal DAPI original registrado. El contorno guardado se usará en las tres vistas de sacarosa.').replace('rgba(255,255,255,.66)','rgba(255,255,255,.12)')
    ui=ui.replace('Accept this HIL polygon','Guardar contorno ventricular').replace('Close polygon','Cerrar polígono').replace('Undo point','Deshacer punto').replace('Reopen</button>','Reabrir</button>').replace('Clear canvas','Limpiar dibujo').replace('st.canBuild=j.accepted===j.total','st.canBuild=false')
    ui=ui.replace('<div class="group"><div class="title">Final render</div>','<div class="group" style="display:none"><div class="title">Final render</div>')
    server=ThreadingHTTPServer((args.host,args.port),Handler);server.state=state;server.html=ui
    print('HIL_URL http://'+args.host+':'+str(args.port)+'/',flush=True);print('OUTPUT '+str(output),flush=True)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()
if __name__=='__main__':main()
