"""
load_literature.py - load the guides and single-stranded regions published for XIST in
Bandaru et al. 2020 (Sci Rep 10:11610, Supplementary Information 2, sheets 1 and 2).

Files (in the 'literature' folder next to this script):
    bandaru2020_xist_guides.tsv    146 tested guides (position, length, class, seed bases facing SS RNA)
    bandaru2020_xist_regions.tsv   11 single-stranded regions derived from PARS in the paper

Each guide is checked against the stored XIST sequence (guide must be the reverse complement of the
transcript at the stated position); rows that do not match are skipped and reported.
Only class-level knockdown results are published by the paper, not per-guide values.

Run:  python load_literature.py       (after load_pars_real.py; safe to run again)
"""
import csv
import os
import sys

import mysql.connector
from mysql.connector import Error
from dotenv import load_dotenv

load_dotenv()
HERE      = os.path.dirname(os.path.abspath(__file__))
LIT_DIR   = os.path.join(HERE, 'literature')
NONCODE   = 'NONHSAT137541.2'
CITATION  = 'Bandaru et al. 2020, Sci Rep 10:11610'

DB_CONFIG = {
    'host': os.getenv('DB_HOST', 'localhost'), 'user': os.getenv('DB_USER', 'root'),
    'password': os.getenv('DB_PASSWORD', ''), 'port': int(os.getenv('DB_PORT', '3306')),
    'database': os.getenv('DB_NAME', 'cas13_grna_db'), 'charset': 'utf8mb4',
}
COMP = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


def revcomp(seq):
    return ''.join(COMP[c] for c in reversed(seq))


def main():
    g_file = os.path.join(LIT_DIR, 'bandaru2020_xist_guides.tsv')
    r_file = os.path.join(LIT_DIR, 'bandaru2020_xist_regions.tsv')
    if not (os.path.exists(g_file) and os.path.exists(r_file)):
        print('[SKIP] literature files not found; continuing without literature support')
        return
    conn = None
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        cur = conn.cursor()
        cur.execute('SELECT id, sequence FROM transcripts WHERE noncode_id=%s', (NONCODE,))
        row = cur.fetchone()
        if not row:
            print('[SKIP] {} is not in the database; literature not loaded'.format(NONCODE))
            return
        tid, seq = row
        cur.execute('DELETE FROM literature_guides  WHERE transcript_id=%s AND citation=%s', (tid, CITATION))
        cur.execute('DELETE FROM literature_regions WHERE transcript_id=%s AND citation=%s', (tid, CITATION))

        ok = bad = 0
        with open(g_file) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                a, b, ln = int(r['start_pos']), int(r['end_pos']), int(r['guide_len'])
                guide = r['sequence'].upper().replace('T', 'U')
                if b - a + 1 != ln or revcomp(seq[a - 1:b]) != guide:
                    bad += 1
                    continue
                cur.execute("""INSERT INTO literature_guides
                    (transcript_id, citation, guide_name, start_pos, end_pos, guide_len, target_type, pfs, ss_nt, seed_ss_nt)
                    VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)""",
                    (tid, CITATION, r['guide_name'], a, b, ln, r['target_type'], r['pfs'],
                     int(r['ss_nt']), int(r['seed_ss_nt'])))
                ok += 1
        n_reg = 0
        with open(r_file) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                cur.execute("INSERT INTO literature_regions (transcript_id, citation, start_pos, end_pos, kind) "
                            "VALUES (%s,%s,%s,%s,'SS')", (tid, CITATION, int(r['start_pos']), int(r['end_pos'])))
                n_reg += 1
        conn.commit()
        print('[OK] Literature guides loaded: {} ({} skipped: sequence did not match XIST)'.format(ok, bad))
        print('[OK] Literature SS regions loaded: {}'.format(n_reg))
        if bad:
            sys.exit(1)
    except Error as e:
        print('[ERROR] {}'.format(e))
        if conn:
            conn.rollback()
        sys.exit(1)
    finally:
        if conn and conn.is_connected():
            conn.close()


if __name__ == '__main__':
    main()
