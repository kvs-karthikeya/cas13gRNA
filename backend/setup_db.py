"""
setup_db.py - create the MySQL database and tables.
Cas13 gRNA Designer - real experimental probing data (PARS)

Run:
    python setup_db.py            create anything that is missing
    python setup_db.py --reset    DROP all project tables first (use this once locally
                                  if you set the project up with an older version)

Tables
  probing_datasets  one row per experiment (method, cell type, source, citation)
  transcripts       lncRNAs that have experimental data, with their sequence
  pars_scores       the real score at every position that has enough reads
  literature_*      guides and single-stranded regions published in Bandaru et al. 2020 (XIST)
Positions without reads are simply not stored ("no data").
"""
import os
import sys

import mysql.connector
from mysql.connector import Error
from dotenv import load_dotenv

load_dotenv()

DB_CONFIG = {
    'host':     os.getenv('DB_HOST', 'localhost'),
    'user':     os.getenv('DB_USER', 'root'),
    'password': os.getenv('DB_PASSWORD', ''),
    'port':     int(os.getenv('DB_PORT', '3306')),
    'charset':  'utf8mb4',
}
DB_NAME = os.getenv('DB_NAME', 'cas13_grna_db')

# Drop order matters: tables that point at `transcripts` must go first.
# ss_regions and grna_cache existed in older versions and also point at transcripts.
TABLES_NEWEST_FIRST = ['literature_guides', 'literature_regions', 'grna_cache', 'ss_regions',
                       'pars_scores', 'transcripts', 'probing_datasets']


def main():
    reset = '--reset' in sys.argv
    print('=' * 60)
    print('  Cas13 gRNA Designer - Database Setup')
    print('=' * 60)
    conn = None
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        cur = conn.cursor()
        cur.execute("CREATE DATABASE IF NOT EXISTS `{}` CHARACTER SET utf8mb4 "
                    "COLLATE utf8mb4_unicode_ci".format(DB_NAME))
        cur.execute("USE `{}`".format(DB_NAME))
        print("[OK] Database '{}' ready".format(DB_NAME))

        if reset:
            cur.execute("SET FOREIGN_KEY_CHECKS = 0")      # belt and braces for unusual old layouts
            for t in TABLES_NEWEST_FIRST:
                cur.execute("DROP TABLE IF EXISTS `{}`".format(t))
            cur.execute("SET FOREIGN_KEY_CHECKS = 1")
            print("[OK] Old tables dropped (--reset)")
        else:
            # detect the older layout (it had a column called `condition`)
            cur.execute("SELECT COUNT(*) FROM information_schema.columns "
                        "WHERE table_schema=%s AND table_name='transcripts' AND column_name='condition'",
                        (DB_NAME,))
            if cur.fetchone()[0]:
                print("[ERROR] This database has the OLD table layout.")
                print("        Run once:  python setup_db.py --reset")
                sys.exit(1)

        cur.execute("""
            CREATE TABLE IF NOT EXISTS probing_datasets (
                id                 INT AUTO_INCREMENT PRIMARY KEY,
                name               VARCHAR(200) NOT NULL UNIQUE,
                species            VARCHAR(100) NOT NULL,
                method             VARCHAR(50)  NOT NULL,
                cell_type          VARCHAR(100),
                probing_condition  VARCHAR(100),
                source             VARCHAR(200),
                citation           VARCHAR(300),
                score_definition   VARCHAR(300),
                created_at         TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
        """)
        cur.execute("""
            CREATE TABLE IF NOT EXISTS transcripts (
                id              INT AUTO_INCREMENT PRIMARY KEY,
                noncode_id      VARCHAR(60)  NOT NULL UNIQUE,
                alt_noncode_ids TEXT,
                gene_name       VARCHAR(100),
                source_id       VARCHAR(100),
                species         VARCHAR(100) DEFAULT 'human',
                rna_type        VARCHAR(50)  DEFAULT 'lncRNA',
                transcript_len  INT          NOT NULL,
                sequence        LONGTEXT     NOT NULL,
                created_at      TIMESTAMP    DEFAULT CURRENT_TIMESTAMP,
                INDEX idx_gene_name (gene_name),
                INDEX idx_source_id (source_id)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
        """)
        cur.execute("""
            CREATE TABLE IF NOT EXISTS pars_scores (
                transcript_id  INT    NOT NULL,
                dataset_id     INT    NOT NULL,
                position       INT    NOT NULL,
                pars_score     FLOAT  NOT NULL,
                read_count     INT    NOT NULL,
                PRIMARY KEY (transcript_id, dataset_id, position),
                FOREIGN KEY (transcript_id) REFERENCES transcripts(id) ON DELETE CASCADE,
                FOREIGN KEY (dataset_id)    REFERENCES probing_datasets(id) ON DELETE CASCADE
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
        """)
        cur.execute("""
            CREATE TABLE IF NOT EXISTS literature_guides (
                id             INT AUTO_INCREMENT PRIMARY KEY,
                transcript_id  INT          NOT NULL,
                citation       VARCHAR(200) NOT NULL,
                guide_name     VARCHAR(50),
                start_pos      INT          NOT NULL,
                end_pos        INT          NOT NULL,
                guide_len      INT          NOT NULL,
                target_type    VARCHAR(10),
                pfs            VARCHAR(2),
                ss_nt          INT,
                seed_ss_nt     INT,
                INDEX idx_lit_guide (transcript_id, start_pos),
                FOREIGN KEY (transcript_id) REFERENCES transcripts(id) ON DELETE CASCADE
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
        """)
        cur.execute("""
            CREATE TABLE IF NOT EXISTS literature_regions (
                id             INT AUTO_INCREMENT PRIMARY KEY,
                transcript_id  INT          NOT NULL,
                citation       VARCHAR(200) NOT NULL,
                start_pos      INT          NOT NULL,
                end_pos        INT          NOT NULL,
                kind           VARCHAR(10)  DEFAULT 'SS',
                INDEX idx_lit_region (transcript_id),
                FOREIGN KEY (transcript_id) REFERENCES transcripts(id) ON DELETE CASCADE
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
        """)
        conn.commit()
        print("[OK] Tables ready: probing_datasets, transcripts, pars_scores, literature_guides, literature_regions")
    except Error as e:
        print("[ERROR] MySQL: {}".format(e))
        sys.exit(1)
    finally:
        if conn and conn.is_connected():
            conn.close()


if __name__ == '__main__':
    main()
