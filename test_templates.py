import importlib.util, json, tempfile, unittest
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
s=importlib.util.spec_from_file_location('tool',Path(__file__).with_name('af_server_templates.py'))
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
ROOT=Path(__file__).resolve().parent.parent
class Tests(unittest.TestCase):
 def test_supplied_examples(self):
  examples=ROOT/'upload'
  expected=json.loads(next(examples.glob('*.json')).read_text())
  seq=''.join(m.AA[x] for x in MMCIF2Dict(str(examples/expected[0]['name']))['_entity_poly_seq.mon_id'])
  with tempfile.TemporaryDirectory() as tmp:
   tmp=Path(tmp); config=tmp/'config.json'
   config.write_text(json.dumps({'sequence':seq,'templates':[{'file':str(examples/h['name'])} for h in expected]}))
   reports=m.prepare(config,tmp/'out')
   actual=json.loads((tmp/'out/query_to_hit.json').read_text())
   for a,e,r in zip(actual,expected,reports):
    self.assertEqual(dict(zip(a['queryIndices'],a['templateIndices'])),dict(zip(e['queryIndices'],e['templateIndices'])))
    d=MMCIF2Dict(str(tmp/'out'/a['name']))
    self.assertEqual(len(set(d['_atom_site.label_asym_id'])),1)
    self.assertEqual(len(d['_entity_poly_seq.num']),len(r['template_sequence']))
    self.assertLess(r['mapped_resolved'],r['mapped'])
 def test_manual_and_errors(self):
  pairs,_=m.mapping('ACD','ACD',{'mode':'manual','numbering':'auth','pairs':[[1,'10'],[2,'10A'],[3,'12']]},{'10':0,'10A':1,'12':2})
  self.assertEqual(pairs,[(0,0),(1,1),(2,2)])
  with self.assertRaises(ValueError):m.mapping('ACD','ACD',{'mode':'manual','blocks':[{'query':[1,3],'template':[1,2]}]}, {})
 def test_alignment_deletion(self):
  pairs,_=m.mapping('ACDEFGHIKLMNPQRSTVWY','ACDEFGHIKNPQRSTVWY',{}, {})
  self.assertNotIn(9,[q for q,t in pairs]);self.assertNotIn(10,[q for q,t in pairs])
if __name__=='__main__':unittest.main()
