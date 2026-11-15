"""
load_pars_real.py - load REAL experimental PARS scores into MySQL.
Cas13 gRNA Designer

Input file (made by build_pars_database.py, kept in this folder):
    pars_real_GM12878.tsv.gz
It holds, for every matched lncRNA: NONCODE ID, sequence, and the PARS score and
read count at every position. Nothing in this project is simulated.

Run:
    python load_pars_real.py
Safe to run again: existing rows for the same transcript are replaced.
"""
import gzip
import os
import re
import sys

import mysql.connector
from mysql.connector import Error
from dotenv import load_dotenv

load_dotenv()

HERE = os.path.dirname(os.path.abspath(__file__))
PARS_FILE = os.getenv('PARS_FILE', os.path.join(HERE, 'pars_real_GM12878.tsv.gz'))

DB_CONFIG = {
    'host':     os.getenv('DB_HOST', 'localhost'),
    'user':     os.getenv('DB_USER', 'root'),
    'password': os.getenv('DB_PASSWORD', ''),
    'port':     int(os.getenv('DB_PORT', '3306')),
    'database': os.getenv('DB_NAME', 'cas13_grna_db'),
    'charset':  'utf8mb4',
}

# Describes the experiment the file comes from. Shown to users in the tool.
DATASET = {
    'name':              'GSE50676 GM12878 PARS (V1/S1)',
    'species':           'human',
    'method':            'PARS',
    'cell_type':         'GM12878 (lymphoblastoid)',
    'probing_condition': 'renatured (RNA refolded in vitro)',
    'source':            'GEO GSE50676, samples GSM1226157 (S1) and GSM1226158 (V1)',
    'citation':          'Wan et al. 2014, Nature 505:706-709',
    'score_definition':  'log2((V1+5)/(S1+5)); positive = paired, negative = single-stranded',
}


def read_rows(path):
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        for n, line in enumerate(f, start=2):
            parts = line.rstrip('\n').split('\t')
            if len(parts) != len(header):
                raise ValueError("line {}: expected {} columns, found {}".format(n, len(header), len(parts)))
            yield dict(zip(header, parts))


def main():
    print('=' * 60)
    print('  Loading real PARS data')
    print('=' * 60)
    if not os.path.exists(PARS_FILE):
        print("[ERROR] File not found: {}".format(PARS_FILE))
        sys.exit(1)

    conn = None
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        cur = conn.cursor()

        cols = ', '.join(DATASET.keys())
        ph = ', '.join(['%s'] * len(DATASET))
        cur.execute("INSERT INTO probing_datasets ({}) VALUES ({}) "
                    "ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), "
                    "cell_type=VALUES(cell_type), probing_condition=VALUES(probing_condition), "
                    "source=VALUES(source), citation=VALUES(citation), "
                    "score_definition=VALUES(score_definition)".format(cols, ph),
                    tuple(DATASET.values()))
        dataset_id = cur.lastrowid
        print("[OK] Dataset: {}  (id {})".format(DATASET['name'], dataset_id))

        n_tr = n_scores = 0
        for row in read_rows(PARS_FILE):
            ids = row['noncode_id'].split(',')
            seq = row['sequence'].upper().replace('T', 'U')
            length = int(row['length'])
            scores = row['pars_scores'].split(';')
            reads = row['total_reads'].split(';')

            # ---- integrity checks: stop instead of loading something wrong ----
            if not re.fullmatch(r'[ACGU]+', seq):
                raise ValueError("{}: sequence contains letters other than A,C,G,U".format(ids[0]))
            if not (len(seq) == length == len(scores) == len(reads)):
                raise ValueError("{}: length mismatch (sequence {}, scores {}, reads {}, header {})".format(
                    ids[0], len(seq), len(scores), len(reads), length))

            gene = row.get('gene_name') or row['source_id']
            cur.execute("""
                INSERT INTO transcripts (noncode_id, alt_noncode_ids, gene_name, source_id, transcript_len, sequence)
                VALUES (%s,%s,%s,%s,%s,%s)
                ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), alt_noncode_ids=VALUES(alt_noncode_ids),
                    gene_name=VALUES(gene_name), source_id=VALUES(source_id),
                    transcript_len=VALUES(transcript_len), sequence=VALUES(sequence)
            """, (ids[0], ','.join(ids[1:]) or None, gene, row['source_id'], length, seq))
            tid = cur.lastrowid

            cur.execute("DELETE FROM pars_scores WHERE transcript_id=%s AND dataset_id=%s", (tid, dataset_id))
            batch = []
            for i, (sc, rd) in enumerate(zip(scores, reads), start=1):
                if sc != 'NA':
                    batch.append((tid, dataset_id, i, float(sc), int(float(rd))))
            for k in range(0, len(batch), 5000):
                cur.executemany("INSERT INTO pars_scores (transcript_id, dataset_id, position, pars_score, read_count) "
                                "VALUES (%s,%s,%s,%s,%s)", batch[k:k + 5000])
            n_tr += 1
            n_scores += len(batch)

        conn.commit()
        print("[OK] Transcripts loaded : {}".format(n_tr))
        print("[OK] Scored positions   : {}".format(n_scores))

        cur.execute("SELECT noncode_id, gene_name, transcript_len FROM transcripts WHERE noncode_id='NONHSAT137541.2'")
        x = cur.fetchone()
        if x:
            print("[OK] XIST present: {} ({} nt)".format(x[0], x[2]))
    except (Error, ValueError) as e:
        print("[ERROR] {}".format(e))
        if conn:
            conn.rollback()
        sys.exit(1)
    finally:
        if conn and conn.is_connected():
            conn.close()


if __name__ == '__main__':
    main()
