#!/usr/bin/env python3
"""Prepare single-chain mmCIF templates and AlphaFold Server mapping JSON."""
import argparse, json, re, sys
from pathlib import Path
from Bio import Align
from Bio.PDB import PDBParser, MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1_extended as AA


def sequence(value, base):
    if isinstance(value, dict):
        lines = (base / value['fasta']).read_text().splitlines()
        if sum(x.startswith('>') for x in lines) > 1:
            raise ValueError('Provide a single-record FASTA')
        value = ''.join(x for x in lines if not x.startswith('>'))
    value = ''.join(value.split()).upper()
    if not value or set(value) - set('ACDEFGHIKLMNPQRSTVWYX'):
        raise ValueError('Sequence must contain standard amino-acid letters or X')
    return value


def load(path):
    if path.suffix.lower() not in ('.pdb', '.ent'):
        return MMCIF2Dict(str(path))
    # PDB coordinates alone cannot reliably establish unresolved polymer residues.
    structure = PDBParser(QUIET=True).get_structure('template', str(path))
    io = MMCIFIO(); io.set_structure(structure)
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / 'converted.cif'; io.save(str(out))
        return MMCIF2Dict(str(out))


def extract(d, chain, model):
    labels = d['_atom_site.label_asym_id']
    models = d.get('_atom_site.pdbx_PDB_model_num', ['1'] * len(labels))
    available = sorted(set(labels))
    if chain is None:
        candidates = sorted({labels[i] for i, x in enumerate(d['_atom_site.label_seq_id']) if x not in ('.', '?')})
        if len(candidates) != 1:
            raise ValueError(f'Specify label chain ID; protein candidates: {candidates}')
        chain = candidates[0]
    rows = [i for i, x in enumerate(labels) if x == chain and models[i] == str(model) and d['_atom_site.label_seq_id'][i] not in ('.', '?')]
    if not rows:
        raise ValueError(f'No polymer atoms for label chain {chain}, model {model}; chains: {available}')
    entities = {d['_atom_site.label_entity_id'][i] for i in rows}
    if len(entities) != 1: raise ValueError('Chain spans multiple entities')
    entity = next(iter(entities))
    seqrows = [i for i, e in enumerate(d.get('_entity_poly_seq.entity_id', [])) if e == entity]
    monomers = {}
    for i in seqrows:
        n = int(d['_entity_poly_seq.num'][i]); mon = d['_entity_poly_seq.mon_id'][i]
        if n in monomers and monomers[n] != mon: raise ValueError('Microheterogeneous polymer sequence is unsupported')
        monomers[n] = mon
    if not monomers:
        for i in rows: monomers[int(d['_atom_site.label_seq_id'][i])] = d['_atom_site.label_comp_id'][i]
    if sorted(monomers) != list(range(1, max(monomers)+1)):
        raise ValueError('Polymer sequence must have consecutive label_seq_id values starting at 1')
    seq = ''.join(AA.get(monomers[n].upper(), 'X') for n in sorted(monomers))
    auth = {}; observed = set()
    for i in rows:
        t = int(d['_atom_site.label_seq_id'][i])-1; observed.add(t)
        ins = d.get('_atom_site.pdbx_PDB_ins_code', ['?']*len(labels))[i]
        key = d['_atom_site.auth_seq_id'][i] + ('' if ins in ('.','?') else ins)
        if key in auth and auth[key] != t: raise ValueError('Ambiguous author residue numbering')
        auth[key] = t
    out = {'data_':'template', '_entry.id':['template']}
    for k, v in d.items():
        if k.startswith('_atom_site.'): out[k] = [v[i] for i in rows]
        elif k.startswith('_chem_comp.'): out[k] = v
    out['_atom_site.pdbx_PDB_model_num'] = ['1']*len(rows)
    out['_struct_asym.id']=[chain]; out['_struct_asym.entity_id']=[entity]
    out['_entity.id']=[entity]; out['_entity.type']=['polymer']
    out['_entity_poly.entity_id']=[entity]; out['_entity_poly.type']=['polypeptide(L)']
    out['_entity_poly.pdbx_strand_id']=[d['_atom_site.auth_asym_id'][rows[0]]]
    out['_entity_poly.pdbx_seq_one_letter_code']=[seq]
    out['_entity_poly.pdbx_seq_one_letter_code_can']=[seq]
    out['_entity_poly_seq.entity_id']=[entity]*len(seq)
    out['_entity_poly_seq.num']=[str(n) for n in sorted(monomers)]
    out['_entity_poly_seq.mon_id']=[monomers[n] for n in sorted(monomers)]
    out['_entity_poly_seq.hetero']=['n']*len(seq)
    return out, seq, auth, observed, chain, bool(seqrows)


def mapping(query, template, spec, auth):
    mode = spec.get('mode','auto')
    if mode == 'auto':
        aligner = Align.PairwiseAligner(mode='global')
        aligner.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')
        aligner.open_gap_score=-10; aligner.extend_gap_score=-0.5
        aligner.end_gap_score=0
        aligns = aligner.align(query, template); first = aligns[0]
        pairs = [(int(q),int(t)) for q,t in zip(*first.indices) if q >= 0 and t >= 0]
        try: aligns[1]; ambiguous=True
        except IndexError: ambiguous=False
        return pairs, {'alignment_score':first.score, 'multiple_optimal_alignments':ambiguous}
    if mode != 'manual': raise ValueError('mode must be auto or manual')
    def idx(x):
        if spec.get('numbering','label') == 'auth':
            if str(x) not in auth: raise ValueError(f'Author residue {x} has no coordinates; use label numbering for unresolved residues')
            return auth[str(x)]
        return int(x)-1
    pairs = [(int(q)-1,idx(t)) for q,t in spec.get('pairs',[])]
    for block in spec.get('blocks',[]):
        q0,q1=block['query']; t0,t1=block['template']
        if spec.get('numbering','label') == 'auth':
            t0,t1=idx(t0)+1,idx(t1)+1
        if q1-q0 != t1-t0: raise ValueError('Manual blocks must have equal lengths; split blocks at insertions/deletions')
        pairs.extend((q0+i-1,t0+i-1) for i in range(q1-q0+1))
    return pairs, {}


def prepare(config_path, output):
    base = config_path.parent; cfg = json.loads(config_path.read_text())
    query = sequence(cfg['sequence'],base); specs=cfg['templates']
    if not 1 <= len(specs) <= 4: raise ValueError('AlphaFold Server requires 1–4 templates per query chain')
    prepared=[]; report=[]
    for n,spec in enumerate(specs):
        source=(base/spec['file']).resolve()
        d,seq,auth,observed,chain,full=extract(load(source),spec.get('chain'),spec.get('model',1))
        pairs,info=mapping(query,seq,spec,auth); pairs=sorted(pairs)
        if not pairs: raise ValueError('Empty mapping')
        if len({q for q,t in pairs}) != len(pairs) or len({t for q,t in pairs}) != len(pairs): raise ValueError('Duplicate mapping indices')
        if any(q<0 or q>=len(query) or t<0 or t>=len(seq) for q,t in pairs): raise ValueError('Mapping index out of range')
        if any(a[1]>=b[1] for a,b in zip(pairs,pairs[1:])): raise ValueError('Mapping must preserve sequence order')
        if not any(t in observed for q,t in pairs): raise ValueError('Mapping contains no resolved residues')
        name=f'template_{n}.cif'
        hit={'name':name,'queryIndices':[q for q,t in pairs],'templateIndices':[t for q,t in pairs]}
        report.append(dict(source=str(source),label_chain=chain,template_sequence=seq,full_polymer_sequence_available=full,
            mapped=len(pairs),query_coverage=len(pairs)/len(query),identity=sum(query[q]==seq[t] for q,t in pairs)/len(pairs),
            mapped_resolved=sum(t in observed for q,t in pairs),**info,
            residue_pairs=[{'query_position':q+1,'template_label_seq_id':t+1,'query_aa':query[q],'template_aa':seq[t],'resolved':t in observed} for q,t in pairs]))
        prepared.append((d,hit))
    output.mkdir(parents=True,exist_ok=True)
    paths=[output/x[1]['name'] for x in prepared]+[output/'query_to_hit.json',output/'mapping_report.json']
    if any(p.exists() for p in paths): raise ValueError('Output files already exist; choose a fresh output directory')
    for d,hit in prepared:
        io=MMCIFIO(); io.set_dict(d); io.save(str(output/hit['name']))
    (output/'query_to_hit.json').write_text(json.dumps([h for d,h in prepared],indent=2)+'\n')
    (output/'mapping_report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('config',type=Path); p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    try:
        report=prepare(args.config.resolve(),args.out.resolve())
        for r in report:
            print(f"{r['label_chain']}: mapped {r['mapped']}, identity {r['identity']:.1%}, resolved {r['mapped_resolved']}")
            if not r['full_polymer_sequence_available']: print('WARNING: coordinate-only sequence; unresolved residues were not inferred',file=sys.stderr)
            if r.get('multiple_optimal_alignments'): print('WARNING: alignment has tied optimal solutions; review mapping_report.json',file=sys.stderr)
    except (ValueError,KeyError,OSError) as e: p.exit(2,f'Error: {e}\n')

if __name__=='__main__': main()
