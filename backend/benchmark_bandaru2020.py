"""
benchmark_bandaru2020.py - re-run the cross-check against Bandaru et al. 2020 after ANY change.

What it does
  Sends the full XIST transcript to your running backend (/analyze) and asks: of the 28-nt guides that
  the paper tested in the lab, how many does the tool propose?
    - "SS-seed" guides: seed fully facing single-stranded RNA (paper: mean relative expression ~0.35,
      strong knockdown)
    - "other-seed" guides: no seed base facing single-stranded RNA (paper: ~0.78, weak knockdown)
  A useful tool proposes SS-seed guides more often than other-seed guides, and clearly more often
  than the share of all windows it proposes ("chance").

Usage (backend must be running, Cas13 variant LshCas13a):
    python benchmark_bandaru2020.py
    python benchmark_bandaru2020.py --url http://localhost:5000 --tool rnafold
Needs: pars_real_GM12878.tsv.gz (for the XIST sequence) and literature/bandaru2020_xist_guides.tsv
"""
import argparse
import csv
import gzip
import json
import os
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--url', default='http://localhost:5000')
ap.add_argument('--tool', default='rnafold')
args = ap.parse_args()

seq = None
with gzip.open(os.path.join(HERE, 'pars_real_GM12878.tsv.gz'), 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        r = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if r['noncode_id'].split(',')[0] == 'NONHSAT137541.2':
            seq = r['sequence']
            break
if not seq:
    raise SystemExit('XIST (NONHSAT137541.2) not found in pars_real_GM12878.tsv.gz')

guides = []
with open(os.path.join(HERE, 'literature', 'bandaru2020_xist_guides.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if int(r['guide_len']) == 28:
            guides.append((int(r['start_pos']), int(r['seed_ss_nt'])))

print('Sending XIST ({} nt) to {} ... (can take a minute)'.format(len(seq), args.url))
req = urllib.request.Request(args.url + '/analyze', data=json.dumps(
    {'sequence': seq, 'tool': args.tool, 'variant': 'lshcas13a'}).encode(),
    headers={'Content-Type': 'application/json'})
d = json.load(urllib.request.urlopen(req, timeout=900))
proposed = {g['position_start'] for g in d['grnas']}
chance = 100.0 * len(proposed) / (len(seq) - 27)

ss  = [a for a, s in guides if s == 8]
oth = [a for a, s in guides if s == 0]
hit = lambda L: sum(1 for a in L if a in proposed)
print('\nStructure method      : {}'.format(d.get('tool_used')))
print('Guides proposed       : {} ({:.1f}% of all possible windows = chance level)'.format(len(proposed), chance))
print('Paper SS-seed guides  : {} of {} proposed ({:.0f}%)'.format(hit(ss), len(ss), 100.0 * hit(ss) / len(ss)))
print('Paper other-seed      : {} of {} proposed ({:.0f}%)'.format(hit(oth), len(oth), 100.0 * hit(oth) / len(oth)))
print('\nRead it like this: SS-seed should be clearly above both the other-seed share and the chance level.')
print('(60 vs 60 guides, one transcript: small differences are within noise.)')
