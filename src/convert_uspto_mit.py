#!/usr/bin/env python
"""Convert mapped USPTO-MIT/USPTO-480K files to MARS raw.csv layout."""
import argparse, csv, json, os, sys

SPLITS=(('train','train'),('valid','val'),('test','test'))

def detokenize(s): return ''.join(str(s).strip().split())

def normalize_reaction(rxn):
    rxn=detokenize(rxn)
    if '>>' in rxn:
        a,b=rxn.split('>>',1); return a+'>>'+b
    parts=rxn.split('>')
    if len(parts)==3: return parts[0]+'>>'+parts[2]
    raise ValueError('reaction must be reactants>>product or reactants>reagents>product')

def parse_source(src, source_mode):
    src=detokenize(src)
    if source_mode=='separated' and '>' in src:
        return src.split('>',1)[0]
    return src.replace('>','.') if source_mode=='mixed' else src

def validate_mapped_rxn(rxn):
    try: reactants, product = rxn.split('>>',1)
    except ValueError: return False,'bad_reaction_format'
    from rdkit import Chem
    r=Chem.MolFromSmiles(reactants); p=Chem.MolFromSmiles(product)
    if r is None or p is None: return False,'parse_failure'
    product_maps=[a.GetAtomMapNum() for a in p.GetAtoms() if a.GetAtomicNum()>1]
    if not product_maps or all(m==0 for m in product_maps): return False,'unmapped_product'
    reactant_maps={a.GetAtomMapNum() for a in r.GetAtoms() if a.GetAtomMapNum()!=0}
    missing=sorted(set(m for m in product_maps if m!=0)-reactant_maps)
    if missing: return False,'product_map_missing_in_reactants:'+','.join(map(str,missing[:10]))
    product_unmapped=sum(1 for a in p.GetAtoms() if a.GetAtomicNum()>1 and a.GetAtomMapNum()==0)
    reactant_unmapped=sum(1 for a in r.GetAtoms() if a.GetAtomicNum()>1 and a.GetAtomMapNum()==0)
    if product_unmapped>reactant_unmapped: return False,'unmapped_product_atoms_exceed_reactant'
    return True,''

def read_csv_split(path, split, typed):
    import pandas as pd
    df=pd.read_csv(path)
    rxn_col='rxn_smiles' if 'rxn_smiles' in df.columns else 'reaction_smiles' if 'reaction_smiles' in df.columns else None
    if rxn_col is None: raise ValueError(f'{path} needs rxn_smiles or reaction_smiles')
    for i,row in df.iterrows():
        cls=int(row['class']) if typed and 'class' in df.columns else int(row.get('class',1) if 'class' in df.columns else 1)
        yield {'rxn_smiles':normalize_reaction(row[rxn_col]),'class':cls,'source_id':row.get('source_id',f'{split}:{i}'),'source_split':split,'source_format':'csv'}

def read_txt_split(src_path,tgt_path,split,input_task,source_mode):
    with open(src_path) as fs, open(tgt_path) as ft:
        for i,(src,tgt) in enumerate(zip(fs,ft)):
            src=detokenize(src); tgt=detokenize(tgt)
            if input_task=='forward': rxn=parse_source(src,source_mode)+'>>'+tgt
            else: rxn=tgt+'>>'+src
            yield {'rxn_smiles':rxn,'class':1,'source_id':f'{split}:{i}','source_split':split,'source_format':'src_tgt'}

def find_csv(input_dir, mars_split, mt_split):
    for name in (f'{mars_split}.csv',f'{mt_split}.csv',f'raw_{mars_split}.csv'):
        p=os.path.join(input_dir,name)
        if os.path.exists(p): return p
    return None

def convert(args):
    os.makedirs(args.output_dir,exist_ok=True); manifest={'splits':{},'remapped':False}
    if args.remap_with_rxnmapper:
        manifest['remapped']=True; print('WARNING: remapped results are not directly comparable to original mapped-data benchmarks.', file=sys.stderr)
    for mars_split, mt_split in SPLITS:
        rows=[]; invalid=[]; counts={'input_rows':0,'valid_rows':0,'invalid_rows':0,'unmapped_rows':0,'parse_failures':0}
        csv_path=find_csv(args.input_dir,mars_split,mt_split)
        if csv_path: iterator=read_csv_split(csv_path,mars_split,args.typed)
        else:
            src=os.path.join(args.input_dir,f'src-{mt_split}.txt'); tgt=os.path.join(args.input_dir,f'tgt-{mt_split}.txt')
            if not (os.path.exists(src) and os.path.exists(tgt)): raise FileNotFoundError(f'Missing CSV or {src}/{tgt}')
            iterator=read_txt_split(src,tgt,mars_split,args.input_task,args.source_mode)
        for row in iterator:
            counts['input_rows']+=1
            ok,reason=validate_mapped_rxn(row['rxn_smiles'])
            if ok: rows.append(row); counts['valid_rows']+=1
            else:
                bad=dict(row); bad['invalid_reason']=reason; invalid.append(bad); counts['invalid_rows']+=1
                if reason.startswith('unmapped'): counts['unmapped_rows']+=1
                if reason=='parse_failure': counts['parse_failures']+=1
        frac=counts['invalid_rows']/counts['input_rows'] if counts['input_rows'] else 0
        max_frac=0.0 if args.strict_map and args.max_invalid_fraction is None else (args.max_invalid_fraction if args.max_invalid_fraction is not None else 0.01)
        out_raw=os.path.join(args.output_dir,mars_split,'raw'); os.makedirs(out_raw,exist_ok=True)
        import pandas as pd
        pd.DataFrame(rows,columns=['rxn_smiles','class','source_id','source_split','source_format']).to_csv(os.path.join(out_raw,'raw.csv'),index=False)
        pd.DataFrame(invalid).to_csv(os.path.join(out_raw,'raw_invalid.csv'),index=False)
        counts['invalid_fraction']=frac; manifest['splits'][mars_split]=counts
        with open(os.path.join(out_raw,'conversion_manifest.json'),'w') as f: json.dump(counts,f,indent=2)
        if counts['invalid_rows'] and (args.strict_map or frac>max_frac):
            raise SystemExit(f'{mars_split}: {counts["invalid_rows"]} invalid/unmapped rows. MARS requires atom-mapped reaction SMILES; refusing to continue.')
    with open(os.path.join(args.output_dir,'conversion_manifest.json'),'w') as f: json.dump(manifest,f,indent=2)
    print(json.dumps(manifest,indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input_dir',required=True); p.add_argument('--output_dir',default='data/USPTO480K')
    p.add_argument('--input_task',choices=['forward','retro'],default='forward')
    p.add_argument('--source_mode',choices=['separated','mixed'],default='separated')
    p.add_argument('--typed',action='store_true'); p.add_argument('--strict_map',action='store_true')
    p.add_argument('--max_invalid_fraction',type=float,default=None)
    p.add_argument('--remap_with_rxnmapper',action='store_true',help='Reserved; enables explicit non-comparable remapping mode.')
    args=p.parse_args(); convert(args)
if __name__=='__main__': main()
