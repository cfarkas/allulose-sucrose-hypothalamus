"""Check soma-driven segmentation and persisted human mask revisions."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from PIL import Image

PAPER=next(p for p in Path(__file__).resolve().parents if (p/'Fig5/08_review_microglia_masks.py').is_file())
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,PAPER/file);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
segmenter=load('seg_test','Fig5/07_refine_microglia_masks.py')
reviewer=load('review_test','Fig5/08_review_microglia_masks.py')
trainer=load('train_test','Fig5/06_retrain_microglia_classes.py')

class SegmentationTests(unittest.TestCase):
    def test_nuclei_do_not_create_cells_without_iba1(self):
        dapi=np.zeros((96,96),np.int32);dapi[30:35,30:35]=1
        cells,somas,*_=segmenter.segment(np.zeros_like(dapi,float),dapi)
        self.assertFalse(cells.any());self.assertFalse(somas.any())

    def test_signal_centres_seed_distinct_cells_and_unstained_nucleus_is_excluded(self):
        yy,xx=np.indices((96,96));signal=20+180*np.exp(-((yy-35)**2+(xx-30)**2)/18)+160*np.exp(-((yy-35)**2+(xx-65)**2)/18)
        dapi=np.zeros((96,96),np.int32);dapi[34:38,31:34]=1;dapi[34:38,66:69]=2;dapi[70:76,70:76]=3
        cells,somas,*_=segmenter.segment(signal,dapi)
        self.assertGreater(cells[35,30],0);self.assertGreater(cells[35,65],0)
        self.assertNotEqual(cells[35,30],cells[35,65]);self.assertEqual(cells[73,73],0)
        self.assertTrue(np.all(cells[somas>0]==somas[somas>0]))

class MaskReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='fig5-mask-test-');self.root=Path(self.temp.name)
        self.source=self.root/'source';self.source.mkdir();self.batch=self.root/'review'
        uid='fixture__s01__iba1sig_00001';self.uid=uid
        cell=np.zeros((96,96),bool);cell[43:54,43:54]=True;soma=np.zeros_like(cell);soma[46:51,46:51]=True
        patch=np.zeros((96,96,3),np.uint8);patch[...,0]=np.arange(96,dtype=np.uint8)[None,:];patch[...,1]=cell*255
        paths={}
        for kind in ['microglia','soma','nucleus']:
            p=self.source/f'{kind}_patches'/'A'/'S01'/f'{uid}.png';p.parent.mkdir(parents=True)
            Image.fromarray(patch if kind=='microglia' else soma.astype(np.uint8)*255).save(p);paths[kind]=str(p)
        data={'cell_uid':uid,'sample':'A','section':'S01','animal_id':'A','animal_key':'a','patch_path':paths['microglia'],
              'cell_qc_pass':True,'soma_patch_path':paths['soma'],'nucleus_patch_path':paths['nucleus'],
              'patch_origin_y':0,'patch_origin_x':0,'source_height':96,'source_width':96,'soma_centroid_y':48,'soma_centroid_x':48}
        pd.DataFrame([data]).to_csv(self.source/'per_microglia_cell_measurements.csv',index=False)
        pd.DataFrame([{**data,'region':'ARC','state':'','pseudo_state':'Ramified','pseudo_confidence':.8}]).to_csv(self.source/reviewer.base.TEMPLATE_NAME,index=False)
        reviewer.prepare(self.source,self.batch,self.root/'missing_previous',1,123)
        self.review=reviewer.Review(self.batch,self.batch/reviewer.base.DEFAULT_OUTPUT_NAME,123)

    def tearDown(self):self.temp.cleanup()

    def payload(self):
        c=self.review.cell(self.uid)
        return {'cell_uid':self.uid,'state':'Ramified','reviewer':'UNIT_TEST_ONLY','session':'DISPOSABLE_FIXTURE',
                'expected_patch_sha256':c['patch_sha256'],'expected_mask_revision':c['mask_revision'],'cell_pixels':c['cell_pixels'],'soma_pixels':c['soma_pixels']}

    def test_masks_are_saved_and_classifier_reads_corrected_patch(self):
        payload=self.payload();original=self.source/'microglia_patches/A/S01'/f'{self.uid}.png';source_bytes=original.read_bytes()
        raw_before=np.array(Image.open(original))[...,0]
        payload['cell_pixels'].append(48*96+54)
        self.review.save(payload)
        label=self.review.labels[self.uid]
        self.assertEqual(label['mask_revision'],1);self.assertTrue(label['mask_reviewed'])
        self.assertEqual(original.read_bytes(),source_bytes)
        saved=np.array(Image.open(self.batch/'patches'/f'{self.uid}.png'))
        np.testing.assert_array_equal(saved[...,0],raw_before);self.assertEqual(saved[48,54,1],255)
        meta,cells,_=trainer.load_source(self.batch)
        self.assertEqual(cells.iloc[0].patch_path,str(self.batch/'patches'/f'{self.uid}.png'))
        reviewer.prepare(self.source,self.batch,self.root/'missing_previous',1,123)
        with self.assertRaisesRegex(ValueError,'changed'):
            self.review.save(payload)
        updated=self.payload();updated['soma_pixels'].append(45*96+48)
        with self.review.lock:self.review.apply_mask(updated)
        self.assertNotIn(self.uid,self.review.labels)
        self.assertEqual(self.review.cell(self.uid)['patch_sha256'],updated['expected_patch_sha256'])
        with self.assertRaisesRegex(ValueError,'changed'):
            self.review.save(updated)
        self.assertTrue((self.batch/'mask_edits'/f'{self.uid}.npz').is_file())

    def test_invalid_mask_is_rejected_without_modifying_files(self):
        before=(self.batch/'patches'/f'{self.uid}.png').read_bytes()
        payload=self.payload();payload['soma_pixels'].append(1)
        with self.assertRaisesRegex(ValueError,'inside'):self.review.save(payload)
        self.assertEqual(before,(self.batch/'patches'/f'{self.uid}.png').read_bytes())
        self.assertFalse(self.review.output_csv.exists())

if __name__=='__main__':unittest.main()
