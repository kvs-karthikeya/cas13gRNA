"""
Cas13 gRNA Designer - Complete Production Backend
Based on Bandaru et al. 2020
AI: Groq API (free, fast, no GPU needed)
Run: python app_py_final.py
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
import subprocess
import tempfile
import os
import re
import shutil
import json
import logging
import mysql.connector
from mysql.connector import Error
import requests
from dotenv import load_dotenv
from collections import Counter
from linearfold_integration import (predict_structure_smart, predict_windowed_partition,
                                     run_rnaplfold_unpaired)

load_dotenv()

# ─── LOGGING ─────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*", "methods": ["GET","POST","OPTIONS"], "allow_headers": ["Content-Type"]}})

# ─── GROQ CONFIG ──────────────────────────────────────────────────
GROQ_API_KEY = os.getenv('GROQ_API_KEY', '')
GROQ_MODEL   = os.getenv('GROQ_MODEL',   'llama-3.3-70b-versatile')
GROQ_URL     = 'https://api.groq.com/openai/v1/chat/completions'

# ─── DATABASE CONFIG ──────────────────────────────────────────────
DB_CONFIG = {
    'host':               os.getenv('DB_HOST',     'localhost'),
    'user':               os.getenv('DB_USER',     'root'),
    'password':           os.getenv('DB_PASSWORD', ''),
    'database':           os.getenv('DB_NAME',     'cas13_grna_db'),
    'port':         int(os.getenv('DB_PORT',       '3306')),
    'charset':            'utf8mb4',
    'connection_timeout': 10
}

# Criteria for reading REAL PARS scores (score = log2((V1+5)/(S1+5)); negative = single-stranded).
# These are settings, not facts: agree them with your supervisor and keep them fixed
# before looking at results. Defaults: at least a 2-fold difference between the two
# enzymes, at positions with at least 10 reads.
PARS_SS_CUTOFF        = float(os.getenv('PARS_SS_CUTOFF', '-1.0'))   # score <= this  -> single-stranded
PARS_DS_CUTOFF        = float(os.getenv('PARS_DS_CUTOFF', '1.0'))    # score >= this  -> double-stranded
PARS_MIN_READS        = int(os.getenv('PARS_MIN_READS', '10'))       # reads needed to trust a position
PARS_MIN_SEED_COVERED = int(os.getenv('PARS_MIN_SEED_COVERED', '6')) # seed positions with data needed (of 8)
MAX_SEQ_LEN   = int(os.getenv('MAX_SEQUENCE_LENGTH', '20000'))
SHORT_SEQ_LIM = 1000
MED_SEQ_LIM   = 5000

# ─── CAS13 VARIANTS ───────────────────────────────────────────────
CAS13_VARIANTS = {
    'lshcas13a': {
        'name':        'LshCas13a',
        'full_name':   'LshCas13a (Leptotrichia shahii) — Used in Bandaru et al. 2020',
        'spacer_len':  28,
        'seed_start':  10,
        'seed_end':    17,
        'pfs':         'none',
        'description': 'Exact variant used in Bandaru et al. 2020. 28nt spacer, central seed positions 11-18.'
    },
    'lwacas13a': {
        'name':        'LwaCas13a',
        'full_name':   'LwaCas13a (Leptotrichia wadei)',
        'spacer_len':  28,
        'seed_start':  10,
        'seed_end':    17,
        'pfs':         'non_g',
        'description': 'Alternative Cas13a ortholog. 28nt spacer, central seed 11-18. Non-G PFS.'
    },
    'rfxcas13d': {
        'name':        'RfxCas13d (CasRx)',
        'full_name':   'RfxCas13d / CasRx (Ruminococcus flavefaciens XPD3002)',
        'spacer_len':  23,
        'seed_start':  0,
        'seed_end':    5,
        'pfs':         'none',
        'description': 'Compact Cas13d. 23nt spacer, 5-prime seed 1-6. No PFS.'
    }
}


# ─── GROQ HELPER ──────────────────────────────────────────────────
def call_groq(system_prompt, user_prompt, max_tokens=1000):
    """
    Calls Groq API
    Free tier: 30 req/min, 14400 req/day
    Models: llama-3.1-8b-instant (fastest)
             llama-3.3-70b-versatile (best quality)
             gemma2-9b-it (good for science)
    """
    if not GROQ_API_KEY:
        return None, "No Groq API key. Add GROQ_API_KEY to .env file."

    try:
        headers = {
            'Authorization': f'Bearer {GROQ_API_KEY}',
            'Content-Type':  'application/json'
        }
        payload = {
            'model':    GROQ_MODEL,
            'messages': [
                {'role': 'system', 'content': system_prompt},
                {'role': 'user',   'content': user_prompt}
            ],
            'max_tokens':   max_tokens,
            'temperature':  0.3,
            'stream':       False
        }
        resp = requests.post(
            GROQ_URL, headers=headers,
            json=payload, timeout=30
        )
        if resp.status_code == 200:
            data    = resp.json()
            content = data['choices'][0]['message']['content']
            return content, None
        elif resp.status_code == 429:
            return None, "Rate limit reached. Please wait 1 minute and try again."
        else:
            return None, f"Groq API error {resp.status_code}: {resp.text[:200]}"

    except requests.Timeout:
        return None, "Groq API timeout. Please try again."
    except Exception as e:
        return None, f"Groq error: {str(e)}"


# ─── HEALTH CHECK ─────────────────────────────────────────────────
@app.route('/')
def home():
    rnafold_ok      = check_tool('RNAfold',    ['RNAfold',    '--version'])
    rnastructure_ok = check_tool('Fold',       ['Fold',       '--version'])
    linearfold_ok   = check_tool('linearfold', ['linearfold', '--help'])
    groq_ok         = bool(GROQ_API_KEY)
    db_ok, db_stats = check_database_detailed()

    return jsonify({
        "status":              "Cas13 gRNA Designer API is running",
        "version":             "2.0.0",
        "tools": {
            "RNAfold":      "available" if rnafold_ok      else "not found",
            "RNAstructure": "available" if rnastructure_ok else "not found",
            "LinearFold":   "available" if linearfold_ok   else "not found"
        },
        "ai_agent":            "available" if groq_ok else "no API key",
        "ai_provider":         "Groq",
        "ai_model":            GROQ_MODEL,
        "database":            "connected" if db_ok else "not connected",
        "transcripts_loaded":  db_stats.get('transcripts', 0),
        "pars_scores_loaded":  db_stats.get('pars_scores', 0),
        "pars_datasets_loaded": db_stats.get('probing_datasets', 0),
        "variants":            list(CAS13_VARIANTS.keys()),
        "max_seq_len":         MAX_SEQ_LEN
    })


def check_tool(name, cmd):
    try:
        subprocess.run(cmd, capture_output=True, timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_database():
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        conn.close()
        return True
    except:
        return False


def check_database_detailed():
    """Returns (ok, stats_dict) with row counts for health endpoint."""
    try:
        conn   = mysql.connector.connect(**DB_CONFIG)
        cursor = conn.cursor()
        stats  = {}
        for table in ('transcripts', 'pars_scores', 'probing_datasets'):
            try:
                cursor.execute(f'SELECT COUNT(*) FROM {table}')
                stats[table] = cursor.fetchone()[0]
            except:
                stats[table] = 0
        cursor.close()
        conn.close()
        return True, stats
    except:
        return False, {}


def get_db():
    return mysql.connector.connect(**DB_CONFIG)


# ─── ANALYZE ROUTE ────────────────────────────────────────────────
@app.route('/analyze', methods=['POST'])
def analyze():
    data     = request.get_json()
    sequence = data.get('sequence', '').strip().upper()
    tool     = data.get('tool',     'rnafold')
    variant  = data.get('variant',  'lshcas13a')

    if not sequence:
        return jsonify({"error": "No sequence provided"}), 400

    sequence = sequence.replace('T', 'U').replace('t', 'U')

    if re.search(r'[^AUGC]', sequence):
        return jsonify({"error": "Invalid characters. Only A, U, G, C or T allowed."}), 400
    if len(sequence) < 28:
        return jsonify({"error": "Sequence too short. Minimum 28 nucleotides."}), 400
    if len(sequence) > MAX_SEQ_LEN:
        return jsonify({"error": f"Sequence too long. Maximum {MAX_SEQ_LEN} nucleotides."}), 400

    variant_config = CAS13_VARIANTS.get(variant, CAS13_VARIANTS['lshcas13a'])
    n              = len(sequence)

    logger.info(f'[Analyze] {n}nt | tool={tool} | variant={variant_config["name"]}')

    # Smart prediction
    structure, energy, error, method = predict_structure_smart(sequence, tool)
    if error:
        alt = 'rnastructure' if tool == 'rnafold' else 'rnafold'
        structure, energy, error, method = predict_structure_smart(sequence, alt)
        if error:
            return jsonify({"error": f"Structure prediction failed: {error}"}), 500

    tool_used = build_tool_name(tool, method, n)

    # Probabilities
    pair_probs = {}
    if n <= SHORT_SEQ_LIM:
        if tool == 'rnastructure':
            pair_probs, _ = run_rnastructure_partition(sequence)
        else:
            pair_probs, _ = run_rnafold_partition(sequence)
    elif n <= MED_SEQ_LIM:
        pair_probs, _ = predict_windowed_partition(sequence, tool)

    unpaired_probs = calculate_unpaired_probs(sequence, pair_probs)
    if unpaired_probs:
        prob_source = 'partition function ({})'.format('full sequence' if n <= SHORT_SEQ_LIM else 'windowed')
    else:
        # No full base-pair probabilities (long sequence, or the partition step failed):
        # use RNAplfold local accessibility so seed confidence is still available.
        unpaired_probs, plerr = run_rnaplfold_unpaired(sequence)
        if unpaired_probs:
            prob_source = 'RNAplfold local accessibility (window 200, max span 100)'
        else:
            prob_source = None
            logger.warning('  No unpaired probabilities available: {}'.format(plerr))
    regions        = parse_regions(structure)
    regions        = enrich_regions_with_confidence(regions, unpaired_probs)
    grnas          = design_grnas(sequence, structure, regions, unpaired_probs, variant_config)

    logger.info(f'  Found {len(grnas)} gRNAs | method={method}')

    # Experimental check is part of the core result (not only of the AI agent step)
    pars_validation = validate_grnas_with_pars(sequence, grnas, variant_config)
    evidence        = build_evidence(pars_validation)

    return jsonify({
        "pars_validation": pars_validation,
        "evidence":        evidence,
        "sequence":       sequence,
        "structure":      structure,
        "energy":         energy,
        "tool_used":      tool_used,
        "method":         method,
        "variant":        variant_config,
        "regions":        regions,
        "grnas":          grnas,
        "unpaired_probs": unpaired_probs,
        "stats": {
            "length":             n,
            "ss_regions":         len(regions["ss"]),
            "ds_regions":         len(regions["ds"]),
            "ss_nucleotides":     structure.count('.'),
            "ds_nucleotides":     len(structure) - structure.count('.'),
            "grnas_found":        len(grnas),
            "has_probs":          len(unpaired_probs) > 0,
            "variant_name":       variant_config['name'],
            "spacer_length":      variant_config['spacer_len'],
            "prediction_method":  method,
            "probability_source": prob_source
        }
    })


def build_evidence(pv):
    """One clear label per result: is there experimental structure data, or is it prediction only?"""
    if pv.get('status') == 'success':
        counts = pv.get('summary', {}).get('counts', {})
        ds = pv.get('dataset', {})
        return {
            "tier":  "experimental_data_available",
            "label": "Experimental PARS data available for {} - guides are checked against it. "
                     "The structure itself is predicted without using that data.".format(pv['transcript']['gene_name']),
            "dataset": "{} | {} | {}".format(ds.get('method',''), ds.get('cell_type',''), ds.get('probing_condition','')),
            "counts": counts
        }
    return {
        "tier":  "predicted_only",
        "label": "Predicted only - no experimental structure data is stored for this sequence. "
                 "Treat guide quality as a computational estimate that needs lab testing.",
        "dataset": "", "counts": {}
    }


def build_tool_name(tool, method, n):
    base = "RNAstructure (Mathews Lab) — Used in Bandaru et al. 2020" \
           if tool == 'rnastructure' else "ViennaRNA RNAfold"
    if method == 'linearfold':
        return f"LinearFold O(n) — Fast prediction for {n}nt sequence"
    elif method == 'windowed':
        return f"{base} — Windowed (200nt windows)"
    elif method == 'windowed_fallback':
        return f"{base} — Windowed fallback"
    return f"{base} — Full prediction"


# ─── AI AGENT ROUTE ───────────────────────────────────────────────
@app.route('/ai_agent', methods=['POST'])
def ai_agent():
    data = request.get_json()

    sequence  = data.get('sequence',  '')
    structure = data.get('structure', '')
    grnas     = data.get('grnas',     [])
    regions   = data.get('regions',   {})
    variant   = data.get('variant',   {})
    energy    = data.get('energy',    'N/A')
    stats     = data.get('stats',     {})
    method    = data.get('method',    'full')

    if not sequence or not grnas:
        return jsonify({"error": "No sequence or gRNAs provided"}), 400

    pars_validation = validate_grnas_with_pars(sequence, grnas, variant)
    ai_analysis     = get_groq_analysis(
        sequence, structure, grnas, regions,
        variant, energy, stats, pars_validation, method
    )

    return jsonify({
        "pars_validation": pars_validation,
        "ai_analysis":     ai_analysis
    })


# ─── CHAT ROUTE ───────────────────────────────────────────────────
@app.route('/ai_chat', methods=['POST'])
def ai_chat():
    data = request.get_json()

    question     = data.get('question',        '')
    sequence     = data.get('sequence',        '')
    grnas        = data.get('grnas',           [])
    regions      = data.get('regions',         {})
    variant      = data.get('variant',         {})
    pars_val     = data.get('pars_validation', {})
    chat_history = data.get('chat_history',    [])

    if not question:
        return jsonify({"error": "No question provided"}), 400

    if not GROQ_API_KEY:
        return jsonify({"answer": "Groq API key not configured. Add GROQ_API_KEY to .env file."})

    answer = get_groq_chat(
        question, sequence, grnas,
        regions, variant, pars_val, chat_history
    )
    return jsonify({"answer": answer})


# ─── WHAT THE AI IS TOLD ABOUT THE EVIDENCE ───────────────────────
AI_FACTS = """FACTS ABOUT THIS TOOL (obey them):
- Use only the numbers given below. Never invent scores, positions, percentages or experimental results.
- The RNA structure is a COMPUTATIONAL prediction made without using experimental data. Never call
  predicted accessibility "validated" or "confirmed".
- Faculty score = computational ranking (lower is better; score 1 = no paired bases in the 28 nt window).
  SS% = predicted unpaired fraction of the guide's target window. Seed confidence = mean predicted unpaired
  probability over the seed (may be unavailable for some inputs; say so instead of guessing).
- PARS experimental evidence is separate. It is stored in the MySQL database (GEO GSE50676, GM12878,
  RNA refolded in vitro) and checked for every guide. PARS scores are log2((V1+5)/(S1+5)): negative =
  single-stranded, positive = paired. They are NOT probabilities.
- PARS statuses: confirmed / partial / inconclusive / contradicted are outcomes of the experimental check of the
  seed. no_data means the transcript HAS data but too few seed positions have enough reads. Never say PARS
  was not calculated or not stored unless the evidence line says PREDICTED ONLY.
- A seed rich in C can look paired in PARS because C positions score high in this dataset; mention this caution
  when the data says so.
- Calibration: on XIST, structure prediction without experimental constraints separated the lab-tested working
  guides from non-working guides in Bandaru et al. 2020 only weakly (AUC about 0.6, where 0.5 is chance).
  Never present the ranking as proven; recommend lab testing of several top guides.
- "Tested in Bandaru 2020" means that exact guide was tested; only class-level knockdown results are published.
- Guide positions are positions in the sequence the user pasted. If an offset is given, transcript position =
  shown position + offset."""


def describe_validation(pv, top=5):
    """Plain-text summary of the experimental check, built from the backend response."""
    if pv.get('status') != 'success':
        return "EVIDENCE: PREDICTED ONLY. " + str(pv.get('message', 'No experimental structure data for this sequence.'))
    ds, cr = pv.get('dataset', {}), pv.get('criteria', {})
    sm = pv.get('summary', {})
    c  = sm.get('counts', {})
    tr = pv.get('transcript', {})
    lines = [
        "EVIDENCE: experimental PARS data is available for {} ({}).".format(tr.get('gene_name', ''), tr.get('noncode_id', '')),
        "Dataset: {} | {} | {} | {}".format(ds.get('method', ''), ds.get('cell_type', ''),
                                           ds.get('probing_condition', ''), ds.get('citation', '')),
        "Rule: a position counts if it has >= {} reads; confirmed if the seed mean score <= {}; contradicted if >= {}; "
        "at least {} seed positions needed.".format(cr.get('min_reads'), cr.get('ss_cutoff'), cr.get('ds_cutoff'),
                                                     cr.get('min_seed_covered')),
        "Pasted sequence starts at transcript position {} (offset {}).".format(tr.get('offset', 0) + 1, tr.get('offset', 0)),
        "Result over all {} candidates: confirmed {}, partial {}, inconclusive {}, contradicted {}, no_data {}.".format(
            sm.get('total', '?'), c.get('confirmed', 0), c.get('partial', 0), c.get('inconclusive', 0),
            c.get('contradicted', 0), c.get('no_data', 0)),
        "Guides tested in Bandaru 2020: {}; guides inside the paper's single-stranded regions: {}.".format(
            sm.get('literature_tested', 0), sm.get('in_paper_ss_region', 0)),
        "Top candidates:"
    ]
    lit = pv.get('literature', {}) or {}
    for rank, val in list(pv.get('grnas', {}).items())[:top]:
        extra = ""
        li = lit.get(rank) or lit.get(str(rank)) or {}
        if li.get('tested'):
            extra += " | tested in paper: {} (class mean relative expression {})".format(
                li['tested']['seed_class'], li['tested']['class_mean_relative_expression'])
        if li.get('in_paper_ss_region'):
            extra += " | seed inside a paper single-stranded region"
        lines.append("  rank {}: {} - {} (seed positions with data {}/{}){}".format(
            rank, val['status'], val['message'], val.get('seed_positions_scored', '?'),
            val.get('seed_positions_total', '?'), extra))
    return "\n".join(lines)


# ─── GROQ ANALYSIS ────────────────────────────────────────────────
def get_groq_analysis(sequence, structure, grnas, regions,
                      variant, energy, stats, pars_validation, method):
    n        = stats.get('length', len(sequence))
    ss_count = stats.get('ss_nucleotides', 0)
    ds_count = stats.get('ds_nucleotides', 0)

    ss_info = "\n".join([
        "  pos {}-{} ({}nt){}".format(
            r['start_display'], r['end_display'], r['length'],
            f", confidence: {r['confidence']}%" if r.get('confidence') is not None else ""
        )
        for r in regions.get('ss', [])[:10]
    ])

    grna_info = "\n".join([
        "  Rank {}: {} | pos {}-{} | score: {} | SS%: {}%".format(
            g['rank'], g['grna'],
            g['position_start'], g['position_end'],
            g['faculty_score'], g['ss_percent']
        )
        for g in grnas[:5]
    ])

    pars_info = describe_validation(pars_validation, top=5)

    system_prompt = ("""You are an expert RNA biologist specializing in CRISPR-Cas13 systems with deep knowledge of
Bandaru et al. 2020 on structure-based gRNA design. You write precise, honest, scientific analysis.
Always respond with valid JSON only. No markdown. No extra text.

""" + AI_FACTS)

    user_prompt = f"""Analyze these Cas13 gRNA design results and respond with JSON.

Variant: {variant.get('name', 'LshCas13a')}
Sequence length: {n}nt
Prediction method: {method}
Energy: {energy} kcal/mol (sum of window energies when windowed)
SS nucleotides: {ss_count} ({round(ss_count/n*100,1) if n>0 else 0}%)
DS nucleotides: {ds_count} ({round(ds_count/n*100,1) if n>0 else 0}%)
SS regions: {stats.get('ss_regions',0)}
Total gRNAs: {stats.get('grnas_found',0)}

Structure (first 80nt): {structure[:80]}...

SS Regions:
{ss_info if ss_info else 'None found'}

Top gRNAs:
{grna_info if grna_info else 'None designed'}

SCORING: Lower score = better. Score 1 = perfect. Score 18+ = poor.
a = DS-complementing nucleotides (want low)
b = spacer length

Probability source: {stats.get('probability_source') or 'not available'}

Experimental evidence:
{pars_info}

Respond ONLY with this JSON structure:
{{
  "structure_analysis": "2-3 sentence scientific analysis of SS/DS ratio and best target regions",
  "grna_recommendation": "2-3 sentence evaluation of top gRNAs with recommendation. Be honest if quality is poor.",
  "pars_interpretation": "2-3 sentences on PARS validation results or what experimental validation would show",
  "summary": "2-3 sentences overall assessment and specific next steps for the researcher"
}}"""

    content, error = call_groq(system_prompt, user_prompt, max_tokens=800)

    if error:
        return {
            "structure_analysis":  f"AI error: {error}",
            "grna_recommendation": "",
            "pars_interpretation": "",
            "summary":             ""
        }

    try:
        text = content.strip().replace('```json','').replace('```','').strip()
        return json.loads(text)
    except json.JSONDecodeError:
        return {
            "structure_analysis":  content[:400] if content else "Analysis unavailable",
            "grna_recommendation": "",
            "pars_interpretation": "",
            "summary":             ""
        }


# ─── GROQ CHAT ────────────────────────────────────────────────────
def get_groq_chat(question, sequence, grnas, regions,
                  variant, pars_validation, history):
    n = len(sequence)

    grna_summary = "\n".join([
        f"Rank {g['rank']}: {g['grna']} | pos {g['position_start']}-{g['position_end']} | score: {g['faculty_score']} | SS%: {g['ss_percent']}%"
        for g in grnas[:10]
    ])

    ss_summary = ", ".join([
        f"pos {r['start_display']}-{r['end_display']} ({r['length']}nt)"
        for r in regions.get('ss', [])[:8]
    ])

    pars_summary = describe_validation(pars_validation, top=10) if pars_validation else "EVIDENCE: not available."

    history_text = ""
    if history:
        history_text = "\nRecent conversation:\n"
        for msg in history[-4:]:
            role = "User" if msg['role']=='user' else "Assistant"
            history_text += f"{role}: {msg['content']}\n"

    system_prompt = ("""You are an expert Cas13 gRNA design assistant.
Answer questions accurately using only the data provided.
Be specific about positions, sequences, and scores. Explain WHY a value is missing instead of just saying so.
When asked about the table: Score = computational ranking; SS% = predicted unpaired fraction; Seed Conf = predicted
seed accessibility; PFS = variant compatibility; PARS = separate experimental evidence.
Keep answers to 3-6 sentences. No markdown.

""" + AI_FACTS)

    user_prompt = f"""Context:
Variant: {variant.get('name','LshCas13a')}
Sequence: {n}nt
SS regions: {ss_summary if ss_summary else 'None'}
Total gRNAs: {len(grnas)}
Scoring: LOWER = BETTER (score 1=perfect, 18+=poor)

gRNA candidates:
{grna_summary if grna_summary else 'None designed'}

Experimental evidence:
{pars_summary}
{history_text}
Question: {question}"""

    content, error = call_groq(system_prompt, user_prompt, max_tokens=400)

    if error:
        return f"AI error: {error}"
    return content or "Could not generate response."


# ─── PARS VALIDATION ──────────────────────────────────────────────
def _seed_offsets(variant):
    """0-based seed offsets inside the spacer window for the chosen Cas13 variant.
    Falls back to LshCas13a (positions 11-18) when the variant is unknown."""
    name = (variant or {}).get('name', '')
    for vc in CAS13_VARIANTS.values():
        if vc['name'] == name:
            return vc['seed_start'], vc['seed_end']
    vc = CAS13_VARIANTS['lshcas13a']
    return vc['seed_start'], vc['seed_end']


# Bandaru et al. 2020, Supplementary Table 4: mean XIST expression relative to no-insert control
# (lower = stronger knockdown) for guides grouped by how many central-seed bases face single-stranded RNA.
# Only these class means are published, not per-guide values.
LIT_CITATION   = 'Bandaru et al. 2020, Sci Rep 10:11610'
LIT_SEED_CLASS = {
    8: ('seed fully single-stranded',   '0.35 +/- 0.12 (n=60 guides)'),
    0: ('seed not single-stranded',     '0.78 +/- 0.10 (n=60 guides)'),
}


def literature_support(cursor, tid, offset, grnas, seed_off_start, seed_off_end):
    """Which of our guides were tested in the paper, or sit inside the paper's single-stranded regions."""
    try:
        cursor.execute("SELECT guide_name, start_pos, guide_len, target_type, seed_ss_nt "
                       "FROM literature_guides WHERE transcript_id=%s", (tid,))
        tested = {(r['start_pos'], r['guide_len']): r for r in cursor.fetchall()}
        cursor.execute("SELECT start_pos, end_pos FROM literature_regions WHERE transcript_id=%s", (tid,))
        regions = [(r['start_pos'], r['end_pos']) for r in cursor.fetchall()]
    except Error:
        return {}
    if not tested and not regions:
        return {}
    out = {}
    for g in grnas:
        ps = g['position_start'] + offset
        wl = g['position_end'] - g['position_start'] + 1
        item = {}
        t = tested.get((ps, wl))
        if t:
            cls, mean = LIT_SEED_CLASS.get(t['seed_ss_nt'], ('seed partly single-stranded', 'class mean not published'))
            item['tested'] = {"guide": t['guide_name'], "paper_target_type": t['target_type'],
                              "seed_ss_nt": t['seed_ss_nt'], "seed_class": cls,
                              "class_mean_relative_expression": mean, "citation": LIT_CITATION}
        if regions and all(any(x <= p <= y for x, y in regions)
                           for p in range(ps + seed_off_start, ps + seed_off_end + 1)):
            item['in_paper_ss_region'] = True
        if item:
            out[g['rank']] = item
    return out


def validate_grnas_with_pars(sequence, grnas, variant=None):
    """Check EVERY designed gRNA against REAL experimental PARS scores (read from MySQL).

    Only positions with at least PARS_MIN_READS reads are used. The seed (variant-specific positions,
    same as in gRNA design) is judged by the MEAN PARS score of its informative positions
    (score = log2((V1+5)/(S1+5)); negative = single-stranded, positive = paired):
      no_data       fewer than PARS_MIN_SEED_COVERED seed positions have enough reads
      confirmed     mean <= SS cutoff        (clearly single-stranded)
      partial       SS cutoff < mean < 0     (leans single-stranded)
      inconclusive  0 <= mean < DS cutoff
      contradicted  mean >= DS cutoff        (clearly paired)
    A single noisy position can no longer veto a seed; this matters because C positions score
    high across XIST, so C-rich seeds are flagged with a caution.
    """
    try:
        conn   = get_db()
        cursor = conn.cursor(dictionary=True)
        transcript = find_transcript(cursor, sequence)

        if not transcript:
            return {
                "status":  "no_match",
                "message": "This sequence is not one of the lncRNAs with experimental PARS data in the "
                           "database. Computational prediction only.",
                "grnas":   {}
            }

        tid    = transcript['id']
        offset = transcript['offset']          # 0 when the whole transcript was pasted
        seed_off_start, seed_off_end = _seed_offsets(variant)

        cursor.execute("""
            SELECT d.id, d.name, d.method, d.cell_type, d.probing_condition, d.source, d.citation,
                   COUNT(*) AS n
            FROM pars_scores s JOIN probing_datasets d ON d.id = s.dataset_id
            WHERE s.transcript_id = %s GROUP BY d.id ORDER BY n DESC LIMIT 1
        """, (tid,))
        ds = cursor.fetchone()
        if not ds:
            return {"status": "no_match", "grnas": {},
                    "message": "No experimental PARS scores are stored for this transcript. "
                               "Computational prediction only."}

        cursor.execute(
            "SELECT position, pars_score, read_count FROM pars_scores WHERE transcript_id=%s AND dataset_id=%s",
            (tid, ds['id'])
        )
        scores = {r['position']: r['pars_score'] for r in cursor.fetchall() if r['read_count'] >= PARS_MIN_READS}

        results = {}
        for grna in grnas:
            ps = grna['position_start'] + offset
            pe = grna['position_end'] + offset
            seed_pos  = list(range(ps + seed_off_start, ps + seed_off_end + 1))
            seed_vals = [scores[p] for p in seed_pos if p in scores]
            win_vals  = [scores[p] for p in range(ps, pe + 1) if p in scores]
            covered   = len(seed_vals)
            n_ss = sum(1 for v in seed_vals if v <= PARS_SS_CUTOFF)
            n_ds = sum(1 for v in seed_vals if v >= PARS_DS_CUTOFF)
            sa = round(sum(seed_vals) / covered, 3) if covered else None
            wa = round(sum(win_vals) / len(win_vals), 3) if win_vals else None

            seed_seq = sequence[grna['position_start'] - 1 + seed_off_start:
                                grna['position_start'] - 1 + seed_off_end + 1]
            c_rich = len(seed_seq) > 0 and seed_seq.count('C') / len(seed_seq) >= 0.5

            if covered < PARS_MIN_SEED_COVERED:
                status, conf = "no_data", "Unknown"
                msg = ("Data exists for this transcript, but only {} of {} seed positions have enough reads "
                       "(>= {}); need {}").format(covered, len(seed_pos), PARS_MIN_READS, PARS_MIN_SEED_COVERED)
            elif sa >= PARS_DS_CUTOFF:
                status, conf = "contradicted", "Low"
                msg = "Seed looks double-stranded in the experiment (mean score {} over {} positions)".format(sa, covered)
            elif sa <= PARS_SS_CUTOFF:
                status, conf = "confirmed", "High"
                msg = "Seed looks single-stranded in the experiment (mean score {} over {} positions)".format(sa, covered)
            elif sa < 0:
                status, conf = "partial", "Medium"
                msg = "Seed leans single-stranded (mean score {} over {} positions)".format(sa, covered)
            else:
                status, conf = "inconclusive", "Unknown"
                msg = "Seed is neither clearly single- nor double-stranded (mean score {} over {} positions)".format(sa, covered)
            if c_rich and status in ("contradicted", "inconclusive"):
                msg += ". Caution: C-rich seed; C positions score high in this dataset"

            results[grna['rank']] = {
                "status": status, "message": msg,
                "pars_confirmed": status == "confirmed",
                "seed_pars_avg": sa, "window_pars_avg": wa, "confidence": conf,
                "seed_positions_scored": covered, "seed_positions_total": len(seed_pos),
                "seed_positions_ss": n_ss, "seed_positions_ds": n_ds, "seed_c_rich": c_rich
            }

        literature = literature_support(cursor, tid, offset, grnas, seed_off_start, seed_off_end)
        cursor.close()
        conn.close()
        label = transcript['gene_name'] or transcript['noncode_id']
        return {
            "status":     "success",
            "transcript": {"id": tid, "gene_name": label, "noncode_id": transcript['noncode_id'],
                           "transcript_len": transcript['transcript_len'], "offset": offset},
            "dataset":    {k: ds[k] for k in ('name', 'method', 'cell_type', 'probing_condition', 'source', 'citation')},
            "criteria":   {"min_reads": PARS_MIN_READS, "ss_cutoff": PARS_SS_CUTOFF, "ds_cutoff": PARS_DS_CUTOFF,
                           "min_seed_covered": PARS_MIN_SEED_COVERED},
            "message":    "Checked against experimental {} data for {} ({}, {})".format(
                              ds['method'], label, ds['cell_type'], ds['probing_condition']),
            "summary":    {"total": len(grnas), "counts": dict(Counter(v['status'] for v in results.values())),
                           "literature_tested": sum(1 for v in literature.values() if 'tested' in v),
                           "in_paper_ss_region": sum(1 for v in literature.values() if v.get('in_paper_ss_region'))},
            "grnas":      results,
            "literature": literature
        }

    except Error as e:
        return {"status": "db_error", "message": "Database error: {}".format(e), "grnas": {}}
    except Exception as e:
        return {"status": "error", "message": "Validation error: {}".format(e), "grnas": {}}


def find_transcript(cursor, sequence):
    """Find the stored transcript that contains EXACTLY this sequence.

    Returns a dict with 'offset' (0-based start of the pasted sequence inside the stored
    transcript) or None. There is deliberately no 'closest match' fallback: an unrelated RNA
    must never be validated against another transcript's experimental data.
    """
    seq = sequence.upper().replace('T', 'U')
    cols = "id, gene_name, noncode_id, source_id, transcript_len"

    cursor.execute("SELECT {} FROM transcripts WHERE sequence = %s LIMIT 1".format(cols), (seq,))
    t = cursor.fetchone()
    if t:
        t['offset'] = 0
        return t

    # pasted sequence is a fragment of a stored transcript
    cursor.execute("SELECT {}, LOCATE(%s, sequence) AS pos FROM transcripts "
                   "WHERE LOCATE(%s, sequence) > 0 LIMIT 2".format(cols), (seq, seq))
    rows = cursor.fetchall()
    if len(rows) == 1:
        t = rows[0]
        t['offset'] = t.pop('pos') - 1
        return t
    return None          # not found, or ambiguous (found in several transcripts)


# ─── STRUCTURE PREDICTION ─────────────────────────────────────────
def run_rnafold_partition(sequence):
    tmp_dir = tempfile.mkdtemp()
    try:
        sp = os.path.join(tmp_dir, 'seq.fa')
        with open(sp, 'w') as f:
            f.write(f">seq\n{sequence}\n")
        r = subprocess.run(['RNAfold','--partfunc','--noPS',sp],
                           capture_output=True, text=True, timeout=300, cwd=tmp_dir)
        if r.returncode != 0:
            return {}, "RNAfold partition error"
        dp = os.path.join(tmp_dir, 'seq_dp.ps')
        if not os.path.exists(dp):
            for fn in os.listdir(tmp_dir):
                if fn.endswith('_dp.ps'):
                    dp = os.path.join(tmp_dir, fn); break
        if not os.path.exists(dp):
            return {}, "Dot plot not found"
        return parse_rnafold_dotplot(dp), None
    except Exception as e:
        return {}, str(e)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def parse_rnafold_dotplot(dp_file):
    pp = {}
    try:
        with open(dp_file) as f:
            content = f.read()
        for m in re.findall(r'(\d+)\s+(\d+)\s+([0-9.]+)\s+ubox', content):
            i = int(m[0])-1; j = int(m[1])-1; p = float(m[2])**2
            pp[(i,j)] = p; pp[(j,i)] = p
    except: pass
    return pp


def run_rnastructure_partition(sequence):
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.seq', delete=False) as f:
            f.write(f"; seq\nsequence\n{sequence.replace('U','T')}\n1\n")
            sp = f.name
        pp  = sp.replace('.seq', '.pfs')
        pro = sp.replace('.seq', '_probs.txt')
        r1  = subprocess.run(['partition', sp, pp], capture_output=True, text=True, timeout=300)
        os.unlink(sp)
        if r1.returncode != 0: return {}, "partition error"
        r2 = subprocess.run(['ProbabilityPlot', pp, pro, '--text'], capture_output=True, text=True, timeout=120)
        if os.path.exists(pp): os.unlink(pp)
        if r2.returncode != 0: return {}, "ProbabilityPlot error"
        result = parse_rnastructure_probs(pro)
        if os.path.exists(pro): os.unlink(pro)
        return result, None
    except Exception as e:
        return {}, str(e)


def parse_rnastructure_probs(path):
    pp = {}
    try:
        with open(path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        i = int(parts[0])-1; j = int(parts[1])-1
                        p = 10**(-float(parts[2]))
                        pp[(i,j)] = p; pp[(j,i)] = p
                    except: continue
    except: pass
    return pp


# ─── REGION PARSING ───────────────────────────────────────────────
def parse_regions(structure):
    ss = []; ds = []; i = 0
    while i < len(structure):
        if structure[i] == '.':
            j = i
            while j < len(structure) and structure[j] == '.': j += 1
            ss.append({"start":i,"end":j-1,"start_display":i+1,"end_display":j,"length":j-i})
            i = j
        else:
            j = i
            while j < len(structure) and structure[j] != '.': j += 1
            ds.append({"start":i,"end":j-1,"start_display":i+1,"end_display":j,"length":j-i})
            i = j
    return {"ss": ss, "ds": ds}


def calculate_unpaired_probs(sequence, pair_probs):
    if not pair_probs: return {}
    n = len(sequence)
    return {
        str(i): round(max(0.0, min(1.0, 1.0 - sum(
            pair_probs.get((i,j),0) for j in range(n) if j != i
        ))), 4)
        for i in range(n)
    }


def enrich_regions_with_confidence(regions, unpaired_probs):
    for r in regions["ss"]:
        if unpaired_probs:
            probs = [unpaired_probs.get(str(i),0.5) for i in range(r["start"],r["end"]+1)]
            r["confidence"] = round(sum(probs)/len(probs)*100,1)
            r["confidence_label"] = get_conf_label(r["confidence"])
        else:
            r["confidence"] = None; r["confidence_label"] = "No data"
    for r in regions["ds"]:
        if unpaired_probs:
            probs = [1.0-unpaired_probs.get(str(i),0.5) for i in range(r["start"],r["end"]+1)]
            r["confidence"] = round(sum(probs)/len(probs)*100,1)
            r["confidence_label"] = get_conf_label(r["confidence"])
        else:
            r["confidence"] = None; r["confidence_label"] = "No data"
    return regions


def get_conf_label(c):
    if c>=90: return "Very High"
    if c>=70: return "High"
    if c>=50: return "Medium"
    if c>=30: return "Low"
    return "Very Low"


# ─── gRNA DESIGN ──────────────────────────────────────────────────
def design_grnas(sequence, structure, regions, unpaired_probs, vc):
    results = []
    wl = vc['spacer_len']; ss = vc['seed_start']; se = vc['seed_end']; pfs = vc['pfs']
    n  = len(sequence)

    for i in range(n - wl + 1):
        sw   = structure[i:i+wl]
        seed = sw[ss:se+1]
        if not all(c=='.' for c in seed): continue

        target = sequence[i:i+wl]
        total_ss = sw.count('.')
        ss_pct = round(total_ss/wl*100)
        grna = reverse_complement(target)

        if is_self_complementary(grna): continue
        if forms_dimer(grna):           continue

        pfs_ok = True; pfs_note = "N/A"
        if pfs == 'non_g':
            fp = i + wl
            if fp < n:
                fn = sequence[fp]
                if fn == 'G': pfs_ok = False; pfs_note = f"G-PFS at pos {fp+1}"
                else: pfs_note = f"Non-G PFS: {fn} ✓"
            else: pfs_note = "End of sequence"
        if not pfs_ok: continue

        a  = sw.count('(') + sw.count(')')
        b  = wl
        fs = round(1 - (a*(-1)) - ((28-b)*(-0.5)), 2)

        sc = None
        if unpaired_probs:
            sp = [unpaired_probs.get(str(i+ss+k),0.5) for k in range(se-ss+1)]
            sc = round(sum(sp)/len(sp)*100,1)

        wc = None
        if unpaired_probs:
            wp = [unpaired_probs.get(str(i+k),0.5) for k in range(wl)]
            wc = round(sum(wp)/len(wp)*100,1)

        sas = i+ss; sae = i+se
        tl  = "Multiple SS regions"
        for r in regions["ss"]:
            if r["start"] <= sas and r["end"] >= sae:
                cs = f" | confidence: {r['confidence']}%" if r.get('confidence') is not None else ""
                tl = f"Loop pos {r['start_display']}-{r['end_display']} ({r['length']}nt){cs}"
                break

        results.append({
            "rank":0,"grna":grna,"target":target,
            "position_start":i+1,"position_end":i+wl,
            "ss_percent":ss_pct,"ds_nucleotides":a,
            "faculty_score":fs,"seed_confidence":sc,
            "window_confidence":wc,"struct_window":sw,
            "target_loop":tl,"pfs_note":pfs_note
        })

    results.sort(key=lambda x: x["faculty_score"])
    for idx,r in enumerate(results): r["rank"] = idx+1
    return results


def reverse_complement(seq):
    c = {'A':'U','U':'A','G':'C','C':'G'}
    return ''.join(c.get(b,b) for b in reversed(seq))

def is_self_complementary(seq):
    h = len(seq)//2
    f = seq[:h]; s = seq[len(seq)-h:]
    rc = reverse_complement(f)
    return sum(1 for a,b in zip(rc,s) if a==b)/h > 0.70

def forms_dimer(seq):
    el=6; fi=seq[:el]; th=seq[-el:]
    rc=reverse_complement(th)
    return sum(1 for a,b in zip(fi,rc) if a==b)/el > 0.65


# ─── RUN ──────────────────────────────────────────────────────────
if __name__ == '__main__':
    port  = int(os.getenv('FLASK_PORT', '5000'))
    debug = os.getenv('FLASK_DEBUG','False').lower() == 'true'

    print("=" * 60)
    print("  Cas13 gRNA Designer v2.0")
    print("  Based on Bandaru et al. 2020")
    print("  AI: Groq API (free, fast, no GPU needed)")
    print("=" * 60)

    rf  = check_tool('RNAfold',    ['RNAfold',   '--version'])
    rna = check_tool('Fold',       ['Fold',       '--version'])
    lf  = check_tool('linearfold', ['linearfold', '--help'])
    db  = check_database()

    print(f"\n  RNAfold     : {'✓' if rf  else '✗'}")
    print(f"  RNAstructure: {'✓' if rna else '✗'}")
    print(f"  LinearFold  : {'✓' if lf  else '✗ (optional)'}")
    print(f"  Groq AI     : {'✓' if GROQ_API_KEY else '✗ Add GROQ_API_KEY to .env'}")
    print(f"  MySQL DB    : {'✓' if db  else '✗ Check DB_PASSWORD in .env'}")
    print(f"  Model       : {GROQ_MODEL}")
    print(f"  Max seq     : {MAX_SEQ_LEN}nt")
    print(f"\n  Speed routing:")
    print(f"    <1000nt  → Full prediction")
    print(f"    1-5knt   → Windowed (200nt)")
    print(f"    >5000nt  → LinearFold O(n)")
    print(f"\n  Server: http://localhost:{port}\n")

    app.run(host='0.0.0.0', port=port, debug=debug)