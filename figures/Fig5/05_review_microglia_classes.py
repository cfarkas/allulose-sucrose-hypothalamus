#!/usr/bin/env python3
"""Prepare a reproducible cell-review batch and serve it with supervised retraining."""
from __future__ import annotations
import argparse
import collections
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
from datetime import datetime, timezone
from urllib.parse import urlparse

PAPER = next(p for p in Path(__file__).resolve().parents if (p / 'scripts/shared/05_annotate_microglia.py').is_file())
spec = importlib.util.spec_from_file_location('microglia_review_base', PAPER / 'scripts/shared/05_annotate_microglia.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
base.STATE_HELP['Ramified'] = 'Small soma, long thin branched processes'
base.STATE_HELP['Amoeboid'] = 'Round soma, few or no visible processes'


def atomic_json(path, obj):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(obj, indent=2) + '\n')
    os.replace(temp, path)


def source_patch(analysis, row):
    # Promoted tables can retain deleted temporary build paths. Use their
    # stable sample/section/cell identity to resolve the shipped patch first.
    p = analysis / 'microglia_patches' / row['sample'] / row['section'] / (row['cell_uid'] + '.png')
    if not p.is_file():
        raise ValueError(f'Missing canonical patch: {p}')
    return p


def prepare(analysis, review_dir, count, seed):
    receipt = review_dir / 'review_batch.json'
    if receipt.is_file():
        meta = json.loads(receipt.read_text())
        if meta['source_analysis'] != str(analysis) or meta['cells'] != count or meta['seed'] != seed:
            raise ValueError('Existing review batch differs; use a new review directory.')
        for row in meta['selected_cells']:
            p = review_dir / 'patches' / (row['cell_uid'] + '.png')
            if base.sha256_file(p) != row['patch_sha256']:
                raise ValueError(f'Review patch changed: {p}')
        return meta
    if review_dir.exists() and any(review_dir.iterdir()):
        raise ValueError('A new review directory must be empty.')
    template = analysis / base.TEMPLATE_NAME
    with template.open() as f:
        rows = list(base.csv.DictReader(f))
    if len({r['cell_uid'] for r in rows}) != len(rows):
        raise ValueError('Duplicate cell IDs in source template.')
    for row in rows:
        row['patch_path'] = str(source_patch(analysis, row))
    if count > len(rows):
        raise ValueError(f'Only {len(rows)} cells are available.')
    ordering = base.Review.__new__(base.Review)
    ordering.rows = rows
    by_uid = {r['cell_uid']: r for r in rows}
    selected = [dict(by_uid[u]) for u in ordering._review_order(seed)[:count]]
    review_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / 'patches').mkdir()
    audit = []
    for row in selected:
        original = Path(row['patch_path'])
        target = review_dir / 'patches' / original.name
        shutil.copy2(original, target)
        row['patch_path'] = str(target)
        row['state'] = ''
        audit.append({'cell_uid': row['cell_uid'], 'patch_sha256': base.sha256_file(target),
                      'source_patch': str(original)})
    with (review_dir / base.TEMPLATE_NAME).open('w', newline='') as f:
        writer = base.csv.DictWriter(f, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    meta = {'source_analysis': str(analysis), 'source_template_sha256': base.sha256_file(template),
            'source_measurements_sha256': base.sha256_file(analysis / 'per_microglia_cell_measurements.csv'),
            'cells': count, 'seed': seed,
            'sampling': 'Seeded round-robin over proposed class and animal; proposals are not human labels.',
            'proposed_class_counts': dict(collections.Counter(r['pseudo_state'] for r in selected)),
            'animal_counts': dict(collections.Counter(r['animal_id'] for r in selected)),
            'region_counts': dict(collections.Counter(r['region'] for r in selected)),
            'selected_cells': audit}
    atomic_json(receipt, meta)
    return meta


class Training:
    def __init__(self, review, python, min_labels):
        self.review, self.python, self.min_labels = review, python, min_labels
        self.root = review.analysis_dir
        self.state_path = self.root / 'training_status.json'
        self.lock = threading.Lock()
        self.state = json.loads(self.state_path.read_text()) if self.state_path.is_file() else {'status': 'waiting_for_labels'}
        if self.state['status'] == 'running':
            threading.Thread(target=self._resume_monitor, daemon=True).start()

    def _resume_monitor(self):
        # A server restart must not strand a detached training process.
        import time
        while self.state['status'] == 'running':
            report = Path(self.state.get('output_dir', str(self.root / 'retrained'))) / 'retraining_report.json'
            if report.is_file():
                with self.lock:
                    self.state.update(status='complete', returncode=0, report=json.loads(report.read_text()),
                                      finished_at=datetime.now(timezone.utc).isoformat())
                    atomic_json(self.state_path, self.state)
                return
            try:
                os.kill(int(self.state['pid']), 0)
            except ProcessLookupError:
                with self.lock:
                    self.state.update(status='failed', returncode=None, finished_at=datetime.now(timezone.utc).isoformat())
                    atomic_json(self.state_path, self.state)
                return
            time.sleep(3)

    def info(self):
        counts = self.review.counts()
        deficits = []
        if counts['labelled'] < self.min_labels:
            deficits.append(f"{self.min_labels - counts['labelled']} more labelled cells")
        for name, n in counts['per_class'].items():
            if n < 12:
                deficits.append(f'{12 - n} more {name} labels')
        if counts['animals'] < 3:
            deficits.append('labels from at least 3 animals')
        return {**self.state, 'ready': not deficits, 'needs': deficits,
                'minimum_labels': self.min_labels, 'automatic_at': len(self.review.rows)}

    def start(self):
        with self.lock:
            if self.state['status'] != 'waiting_for_labels':
                raise ValueError('This training run has already been submitted; see its status.')
            info = self.info()
            if not info['ready']:
                raise ValueError('Still needed: ' + ', '.join(info['needs']))
            snapshot = self.root / 'labels_used_for_training.csv'
            with self.review.lock:
                shutil.copy2(self.review.output_csv, snapshot)
            output = self.root / 'retrained'
            command = [self.python, str(PAPER / 'Fig5/06_retrain_microglia_classes.py'),
                       '--review-dir', str(self.root), '--labels-csv', str(snapshot),
                       '--output-dir', str(output), '--min-labels', str(self.min_labels)]
            log = self.root / 'retraining.log'
            with log.open('ab') as stream:
                process = subprocess.Popen(command, cwd=PAPER, stdin=subprocess.DEVNULL,
                                           stdout=stream, stderr=subprocess.STDOUT, start_new_session=True,
                                           env={**os.environ, 'OMP_NUM_THREADS': '4', 'MKL_NUM_THREADS': '4'})
            self.state = {'status': 'running', 'pid': process.pid, 'command': command,
                          'labels': self.review.counts()['labelled'],
                          'labels_sha256': base.sha256_file(snapshot),
                          'started_at': datetime.now(timezone.utc).isoformat(),
                          'output_dir': str(output), 'log': str(log)}
            atomic_json(self.state_path, self.state)
            threading.Thread(target=self._wait, args=(process,), daemon=True).start()
            return self.info()

    def _wait(self, process):
        code = process.wait()
        with self.lock:
            self.state.update(status='complete' if code == 0 else 'failed', returncode=code,
                              finished_at=datetime.now(timezone.utc).isoformat())
            if code == 0:
                report = Path(self.state.get('output_dir', str(self.root / 'retrained'))) / 'retraining_report.json'
                self.state['report'] = json.loads(report.read_text())
            atomic_json(self.state_path, self.state)


class Handler(base.Handler):
    def send_json(self, payload, status=200):
        if isinstance(payload, dict):
            payload['training'] = self.server.training.info()
            if 'gate' in payload:
                payload['gate'] = {'min_cells': self.server.training.min_labels, 'min_per_class': 12,
                                   'min_classes': 4, 'min_animals': 3}
                payload['gate_met'] = payload['training']['ready']
        super().send_json(payload, status)

    def do_GET(self):
        endpoint = urlparse(self.path).path
        if endpoint == '/':
            self.send_bytes(HTML.encode(), 'text/html; charset=utf-8')
        elif endpoint == '/api/training':
            self.send_json({'ok': True, **self.review.counts()})
        else:
            super().do_GET()

    def do_POST(self):
        with self.server.mutation_lock:
            self._post()

    def _post(self):
        endpoint = urlparse(self.path).path
        if endpoint == '/api/retrain':
            if self.headers.get_content_type() != 'application/json':
                self.error_json(415, 'Content-Type must be application/json')
                return
            try:
                self.server.training.start()
                self.send_json({'ok': True, **self.review.counts()})
            except ValueError as exc:
                self.error_json(422, str(exc))
            return
        if endpoint in ('/api/label', '/api/clear'):
            if self.server.training.state['status'] != 'waiting_for_labels':
                self.error_json(409, 'Review is finished; training uses the saved label snapshot.')
                return
        super().do_POST()
        if endpoint == '/api/label' and self.review.counts()['labelled'] == len(self.review.rows):
            if self.server.training.info()['ready'] and self.server.training.state['status'] == 'waiting_for_labels':
                self.server.training.start()


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure 5 · Microglia review</title>
<style>
:root{color-scheme:dark;--bg:#10171d;--panel:#1a252e;--ink:#edf4f5;--muted:#a4b8c4;--accent:#65d6c4;--line:#354650}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,sans-serif}header{padding:18px 28px;border-bottom:1px solid var(--line);display:flex;justify-content:space-between}h1{font-size:20px;margin:0}small,.muted{color:var(--muted)}.layout{display:grid;grid-template-columns:minmax(0,1fr) 320px;min-height:calc(100vh - 75px)}main{padding:22px;max-width:1000px;width:100%;margin:auto}aside{padding:22px;background:var(--panel);border-left:1px solid var(--line)}.views{display:grid;grid-template-columns:1fr 1fr;gap:14px}.views figure{margin:0}canvas{background:#000;width:100%;aspect-ratio:1;image-rendering:pixelated;border:1px solid var(--line);border-radius:6px}figcaption{color:var(--muted);text-align:center;margin-top:4px;font-size:13px}.states{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin:18px 0 12px}button{font:inherit;cursor:pointer;color:var(--ink);background:#22333f;border:1px solid var(--line);padding:10px 14px;border-radius:6px}button:hover{border-color:var(--accent)}button:disabled{opacity:.4;cursor:default}button.state{text-align:left}button.state b{display:block}button.state small{display:block;font-size:12px}button.selected{border:2px solid var(--accent)}button.proposed{outline:1px dashed #efa857}kbd{color:var(--accent);margin-right:8px}.row{display:flex;gap:8px}.row button{flex:1}input[type=text]{width:100%;padding:9px;margin:5px 0 14px;background:#10171d;color:var(--ink);border:1px solid var(--line);border-radius:5px}label{font-size:13px}.group{border-bottom:1px solid var(--line);padding-bottom:18px;margin-bottom:18px}progress{width:100%;accent-color:var(--accent);height:12px}table{width:100%;font-size:13px}td{padding:4px 0}td:last-child{text-align:right}#train{width:100%;background:#20544b}#status{min-height:28px;color:#ffc27d;margin-top:12px}#training{margin:12px 0;font-size:13px}#meta{font-size:13px;color:var(--muted);margin:12px 0}#suggestion{font-size:13px;color:#efa857}.help{font-size:13px}.success{color:var(--accent)}@media(max-width:900px){.layout{grid-template-columns:1fr}aside{border-left:0;border-top:1px solid var(--line)}main{padding:16px}header{padding:14px 16px}}
</style></head><body><header><h1>Figure 5 · Microglia review</h1><span id="headercount" class="muted"></span></header>
<div class="layout"><main><div class="views"><figure><canvas id="signal" width="96" height="96"></canvas><figcaption>Iba1 signal</figcaption></figure><figure><canvas id="overlay" width="96" height="96"></canvas><figcaption>Same signal · target cell outlined in cyan</figcaption></figure></div>
<div id="meta">Loading cell…</div><div class="help">Classify the outlined cell. Each choice saves immediately and opens the next cell.</div><div class="states" id="states"></div><div class="row"><button id="back">Previous cell</button><button id="skip">Unsure · skip (S)</button></div><div id="status" role="status"></div><div id="suggestion"></div></main>
<aside><div class="group"><label for="reviewer">Your name</label><input id="reviewer" type="text" placeholder="Enter your name to save labels" autocomplete="name"><small>Labels are saved locally after every choice. You can close this page and resume later.</small></div><div class="group"><b id="progress"></b><progress id="bar" max="300" value="0"></progress><div id="animals" class="muted"></div><table id="counts"></table></div><div class="group"><button id="train" disabled>Finish review & retrain</button><div id="training"></div><small>Review 200–300 cells. Retraining starts automatically after all 300 are labelled. It uses your labels and holds out separate animals for validation.</small></div><label><input id="showproposal" type="checkbox"> Show model suggestion</label><p class="help muted">Keys 1–4 select a class. S skips uncertain cells. Previous cell lets you correct a choice before finishing.</p><small id="saved"></small></aside></div><script>
'use strict';const $=id=>document.getElementById(id);let cfg,cur=null,busy=false,history=[],session='',imageReady=false,trainingActive=false;
function status(s){$('status').textContent=s||''}
async function api(path,payload){const r=await fetch(path,payload?{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)}:{});const d=await r.json();if(!r.ok||!d.ok)throw Error(d.error||'Request failed');return d}
function proposal(){document.querySelectorAll('.state').forEach(b=>{b.classList.toggle('proposed',$('showproposal').checked&&cur&&b.dataset.state===cur.pseudo_state)});$('suggestion').textContent=$('showproposal').checked&&cur?`Model suggestion: ${cur.pseudo_state} (not a human label)`:''}
function controls(){document.querySelectorAll('.state').forEach(b=>b.disabled=busy||trainingActive||!cur||!imageReady);$('skip').disabled=busy||!cur;$('back').disabled=busy||!history.length}
function render(c){cur=c;imageReady=false;controls();const im=new Image();im.onload=()=>{for(const name of ['signal','overlay']){const canvas=$(name);canvas.width=im.width;canvas.height=im.height;const ctx=canvas.getContext('2d');ctx.drawImage(im,0,0);const data=ctx.getImageData(0,0,im.width,im.height),src=new Uint8ClampedArray(data.data);for(let y=0;y<im.height;y++)for(let x=0;x<im.width;x++){const i=(y*im.width+x)*4,v=src[i];data.data[i]=data.data[i+1]=data.data[i+2]=v;if(name==='overlay'&&src[i+1]>127){const edge=x===0||y===0||x===im.width-1||y===im.height-1||src[i-4+1]<128||src[i+4+1]<128||src[i-im.width*4+1]<128||src[i+im.width*4+1]<128;if(edge){data.data[i]=65;data.data[i+1]=230;data.data[i+2]=230}}}ctx.putImageData(data,0,0)}imageReady=true;controls()};im.onerror=()=>{status('Image could not load. Reload before classifying.');imageReady=false;controls()};im.src='/api/patch?cell_uid='+encodeURIComponent(c.cell_uid);$('meta').textContent=`${c.cell_uid} · ${c.region}`;document.querySelectorAll('.state').forEach(b=>b.classList.toggle('selected',b.dataset.state===c.state));proposal();status(c.state?`Saved as ${c.state}. Select a class to change it.`:'')}
function stats(d){$('headercount').textContent=`${d.labelled} / ${d.total_cells} saved`;$('progress').textContent=`${d.labelled} of ${d.total_cells} cells`;$('bar').max=d.total_cells;$('bar').value=d.labelled;$('animals').textContent=`Labels from ${d.animals} animals`;$('counts').innerHTML=Object.entries(d.per_class).map(([s,n])=>`<tr><td>${s}</td><td>${n}</td></tr>`).join('');const t=d.training;trainingActive=t.status!=='waiting_for_labels';controls();$('train').disabled=busy||!t.ready||t.status!=='waiting_for_labels';if(t.status==='waiting_for_labels')$('training').textContent=t.ready?'Ready to retrain. Continue to 300, or finish now.':'Still needed: '+t.needs.join(', ')+'.';else if(t.status==='running')$('training').textContent='Retraining in progress. Your label snapshot has been saved.';else if(t.status==='complete'){$('training').textContent=`Retraining complete. Held-out animal balanced accuracy: ${(100*t.report.validation_balanced_accuracy).toFixed(1)}%.`;$('training').className='success'}else $('training').textContent='Retraining failed. Your labels are safe; see retraining.log in the review folder.';$('saved').textContent='Session: '+session}
async function next(after=''){const d=await api('/api/next'+(after?'?after='+encodeURIComponent(after):''));stats(d);if(d.done){cur=null;imageReady=false;$('meta').textContent='All cells in this batch are labelled.';status(d.training.ready?'Retraining will start automatically.':'More class coverage is needed before retraining.');controls()}else render(d.cell)}
async function label(s){if(busy||trainingActive||!cur||!imageReady)return;const reviewer=$('reviewer').value.trim();if(!reviewer){status('Enter your name before saving a label.');$('reviewer').focus();return}busy=true;controls();try{const uid=cur.cell_uid;const d=await api('/api/label',{cell_uid:uid,state:s,reviewer,session});history.push(uid);stats(d);await next(uid)}catch(e){status(e.message)}finally{busy=false;controls()}}
$('reviewer').value=localStorage.getItem('fig5-reviewer')||'';$('reviewer').oninput=()=>localStorage.setItem('fig5-reviewer',$('reviewer').value);$('showproposal').onchange=proposal;$('skip').onclick=async()=>{if(busy||!cur)return;busy=true;controls();try{await next(cur.cell_uid)}catch(e){status(e.message)}finally{busy=false;controls()}};$('back').onclick=async()=>{if(busy||!history.length)return;busy=true;controls();try{const d=await api('/api/cell?cell_uid='+encodeURIComponent(history.pop()));stats(d);render(d.cell)}catch(e){status(e.message)}finally{busy=false;controls()}};$('train').onclick=async()=>{busy=true;controls();try{stats(await api('/api/retrain',{}))}catch(e){status(e.message)}finally{busy=false;controls()}};document.addEventListener('keydown',e=>{if(/input|textarea|button/i.test(e.target.tagName)||e.repeat)return;const n=Number(e.key);if(n>=1&&n<=4){e.preventDefault();label(cfg.states[n-1])}if(e.key.toLowerCase()==='s'){$('skip').click();e.preventDefault()}});
(async()=>{try{cfg=await api('/api/config');session=cfg.analysis_dir.split('/').pop();for(const [i,s] of cfg.states.entries()){const b=document.createElement('button');b.className='state';b.dataset.state=s;b.innerHTML=`<b><kbd>${i+1}</kbd>${s}</b><small>${cfg.state_help[s]}</small>`;b.onclick=()=>label(s);$('states').append(b)}stats(cfg);await next();setInterval(async()=>{try{stats(await api('/api/training'))}catch(e){status('Connection lost. Reconnect before labelling.')}},4000)}catch(e){status(e.message)}})();
</script></body></html>'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--analysis-dir', type=Path, default=PAPER / 'analyses/Fig5/results/human_final_run')
    ap.add_argument('--review-dir', type=Path, required=True)
    ap.add_argument('--count', type=int, default=300)
    ap.add_argument('--min-labels', type=int, default=200)
    ap.add_argument('--seed', type=int, default=20260906)
    ap.add_argument('--host', default='127.0.0.1')
    ap.add_argument('--port', type=int, default=8765)
    ap.add_argument('--retrain-python', default='/home/server/anaconda3/envs/paper_apotome_repro/bin/python')
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args()
    if not 120 <= args.min_labels <= args.count:
        ap.error('Require 120 <= min-labels <= count.')
    meta = prepare(args.analysis_dir.resolve(), args.review_dir.resolve(), args.count, args.seed)
    print(json.dumps({k:v for k,v in meta.items() if k != 'selected_cells'}, indent=2), flush=True)
    if args.prepare_only:
        return
    review = base.Review(args.review_dir.resolve(), args.review_dir.resolve() / base.DEFAULT_OUTPUT_NAME, args.seed)
    server = base.ThreadingHTTPServer((args.host, args.port), Handler)
    server.daemon_threads = True
    server.review = review
    server.mutation_lock = threading.Lock()
    server.training = Training(review, args.retrain_python, args.min_labels)
    atomic_json(review.analysis_dir / 'server.json', {'pid': os.getpid(), 'host': args.host,
                'port': server.server_port, 'url': f'http://{args.host}:{server.server_port}', 'cells': len(review.rows)})
    print(f'Review: http://{args.host}:{server.server_port} · {len(review.rows)} cells', flush=True)
    if review.counts()['labelled'] == len(review.rows) and server.training.info()['ready'] and server.training.state['status'] == 'waiting_for_labels':
        server.training.start()
    try:
        server.serve_forever(0.25)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == '__main__':
    main()
