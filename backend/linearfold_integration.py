"""
LinearFold Integration Module
Cas13 gRNA Designer v2.0
Smart structure prediction routing:
  < 1000 nt  → Full RNAfold or RNAstructure
  1000-5000nt → Windowed prediction (200nt windows)
  > 5000 nt  → LinearFold O(n) linear time

If LinearFold is not installed falls back to windowed.
"""

import subprocess
import tempfile
import os
import re
import shutil
import shlex

SHORT_SEQ_LIM = 1000
MED_SEQ_LIM   = 5000
WINDOW_SIZE   = 200
STEP_SIZE     = 100

# LinearFold: '-V' selects the Vienna (thermodynamic) energy model, so the reported number is a
# real free energy in kcal/mol. Without it LinearFold prints a model SCORE (e.g. +800), not kcal/mol.
LINEARFOLD_FLAGS = shlex.split(os.getenv('LINEARFOLD_FLAGS', '-V'))

# Local accessibility (RNAplfold) used when full base-pair probabilities are not available (> 5000 nt)
PLFOLD_WINDOW = int(os.getenv('PLFOLD_WINDOW', '200'))
PLFOLD_SPAN   = int(os.getenv('PLFOLD_SPAN',   '100'))


# ─── SMART ROUTER ──────────────────────────────────────────────────
def predict_structure_smart(sequence, tool='rnafold'):
    """
    Routes structure prediction based on sequence length.

    Returns:
        (structure, energy, error, method)
        method one of: 'full', 'windowed', 'linearfold', 'windowed_fallback'
    """
    n = len(sequence)

    if n < SHORT_SEQ_LIM:
        # Short: full prediction
        if tool == 'rnastructure':
            struct, energy, err = run_rnastructure_fold(sequence, timeout=300)
        else:
            struct, energy, err = run_rnafold(sequence, timeout=300)

        if err and tool != 'rnafold':
            struct, energy, err = run_rnafold(sequence, timeout=300)
        if err:
            return None, None, err, 'full'
        return struct, energy, None, 'full'

    elif n <= MED_SEQ_LIM:
        # Medium: windowed prediction
        struct, energy, err = predict_windowed_structure(sequence, tool)
        if err:
            return None, None, err, 'windowed'
        return struct, energy, None, 'windowed'

    else:
        # Long: try LinearFold first, fall back to windowed
        struct, energy, err = run_linearfold(sequence)
        if err:
            # Fallback to windowed
            struct, energy, err2 = predict_windowed_structure(sequence, tool)
            if err2:
                return None, None, f"LinearFold: {err} | Windowed: {err2}", 'windowed_fallback'
            return struct, energy, None, 'windowed_fallback'
        return struct, energy, None, 'linearfold'


# ─── RNAFOLD ───────────────────────────────────────────────────────
def run_rnafold(sequence, timeout=300):
    """Run ViennaRNA RNAfold for MFE structure."""
    tmp_dir = tempfile.mkdtemp()
    try:
        input_fa = os.path.join(tmp_dir, 'seq.fa')
        with open(input_fa, 'w') as f:
            f.write(f'>seq\n{sequence}\n')

        result = subprocess.run(
            ['RNAfold', '--noPS', input_fa],
            capture_output=True, text=True,
            timeout=timeout, cwd=tmp_dir
        )
        if result.returncode != 0:
            return None, None, f'RNAfold error: {result.stderr[:200]}'

        lines = result.stdout.strip().split('\n')
        # Lines: >header, sequence, structure (energy)
        struct_line = None
        for line in lines:
            if line and line[0] in '.()':
                struct_line = line
                break
            elif line and (line[0] in '.()' or ' ' in line):
                parts = line.split()
                if parts and parts[0] and parts[0][0] in '.()':
                    struct_line = line
                    break

        if not struct_line:
            # Try last non-empty line
            for line in reversed(lines):
                if line.strip() and line.strip()[0] in '.()':
                    struct_line = line
                    break

        if not struct_line:
            return None, None, 'Could not parse RNAfold output'

        # Parse: ".((...)). (-5.60)" or just ".((...))."
        m = re.match(r'([.()[\]{}]+)\s*(?:\(\s*([+-]?\d+\.?\d*)\s*\))?', struct_line.strip())
        if not m:
            return None, None, f'Could not parse structure line: {struct_line}'

        struct = m.group(1)
        energy = float(m.group(2)) if m.group(2) else 0.0

        # Pad or trim to match sequence length
        if len(struct) < len(sequence):
            struct = struct + '.' * (len(sequence) - len(struct))
        elif len(struct) > len(sequence):
            struct = struct[:len(sequence)]

        return struct, energy, None

    except FileNotFoundError:
        return None, None, 'RNAfold not found. Install ViennaRNA package.'
    except subprocess.TimeoutExpired:
        return None, None, f'RNAfold timeout ({timeout}s) for {len(sequence)}nt sequence.'
    except Exception as e:
        return None, None, str(e)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ─── RNASTRUCTURE ──────────────────────────────────────────────────
def run_rnastructure_fold(sequence, timeout=600):
    """Run RNAstructure Fold command."""
    tmp_dir = tempfile.mkdtemp()
    try:
        seq_file    = os.path.join(tmp_dir, 'seq.seq')
        ct_file     = os.path.join(tmp_dir, 'out.ct')
        dot_file    = os.path.join(tmp_dir, 'out.dot')

        # Write sequence in RNAstructure format (DNA)
        dna_seq = sequence.replace('U', 'T')
        with open(seq_file, 'w') as f:
            f.write(f'; seq\nsequence\n{dna_seq}\n1\n')

        result = subprocess.run(
            ['Fold', seq_file, ct_file],
            capture_output=True, text=True,
            timeout=timeout, cwd=tmp_dir
        )
        if result.returncode != 0:
            return None, None, f'Fold error: {result.stderr[:200]}'

        # Convert CT to dot-bracket
        r2 = subprocess.run(
            ['ct2dot', ct_file, '1', dot_file],
            capture_output=True, text=True,
            timeout=60, cwd=tmp_dir
        )
        if r2.returncode != 0:
            # Try parsing CT directly
            struct, energy = parse_ct_file(ct_file, len(sequence))
            if struct:
                return struct, energy, None
            return None, None, f'ct2dot error: {r2.stderr[:200]}'

        # Parse dot-bracket
        struct, energy = parse_dot_file(dot_file)
        if struct:
            return struct, energy, None
        return None, None, 'Could not parse RNAstructure output'

    except FileNotFoundError:
        return None, None, 'RNAstructure Fold not found.'
    except subprocess.TimeoutExpired:
        return None, None, f'RNAstructure timeout ({timeout}s).'
    except Exception as e:
        return None, None, str(e)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def parse_dot_file(path):
    """Parse RNAstructure dot-bracket output file."""
    try:
        with open(path) as f:
            lines = f.readlines()
        struct = ''
        energy = 0.0
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # Try to extract energy from header
                m = re.search(r'ENERGY\s*=\s*([+-]?\d+\.?\d*)', line)
                if m:
                    energy = float(m.group(1))
            elif line and line[0] in '.()':
                struct += line
        return (struct if struct else None), energy
    except Exception:
        return None, 0.0


def parse_ct_file(path, n):
    """Parse CT format file to get dot-bracket structure."""
    try:
        pairs = {}
        energy = 0.0
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    # First line has length and energy
                    try:
                        length = int(parts[0])
                        m = re.search(r'([+-]?\d+\.?\d+)', line)
                        if m:
                            energy = float(m.group(1))
                        continue
                    except ValueError:
                        pass
                if len(parts) >= 5:
                    try:
                        pos    = int(parts[0])
                        paired = int(parts[4])
                        pairs[pos] = paired
                    except ValueError:
                        continue

        struct = []
        for i in range(1, n+1):
            j = pairs.get(i, 0)
            if j == 0:
                struct.append('.')
            elif j > i:
                struct.append('(')
            else:
                struct.append(')')
        return ''.join(struct), energy
    except Exception:
        return None, 0.0


# ─── LINEARFOLD ────────────────────────────────────────────────────
def run_linearfold(sequence, timeout=120):
    """
    Run LinearFold O(n) for very long sequences.
    LinearFold reads sequence from stdin.
    Output: structure energy
    """
    try:
        result = subprocess.run(
            ['linearfold'] + LINEARFOLD_FLAGS,
            input=sequence,
            capture_output=True, text=True,
            timeout=timeout
        )
        if result.returncode != 0 and not result.stdout:
            return None, None, f'LinearFold error: {result.stderr[:200]}'

        output = result.stdout.strip()
        lines  = output.split('\n')
        for line in lines:
            line = line.strip()
            m = re.match(r'([.()[\]{}]+)\s*(?:\(?\s*([+-]?\d+\.?\d*)\s*\)?)?', line)
            if m and m.group(1):
                struct = m.group(1)
                energy = float(m.group(2)) if m.group(2) else 0.0
                # Adjust length
                if len(struct) < len(sequence):
                    struct = struct + '.' * (len(sequence) - len(struct))
                elif len(struct) > len(sequence):
                    struct = struct[:len(sequence)]
                return struct, energy, None

        return None, None, 'Could not parse LinearFold output'

    except FileNotFoundError:
        return None, None, 'LinearFold not installed. Using windowed fallback.'
    except subprocess.TimeoutExpired:
        return None, None, f'LinearFold timeout ({timeout}s).'
    except Exception as e:
        return None, None, str(e)


# ─── WINDOWED STRUCTURE ────────────────────────────────────────────
def predict_windowed_structure(sequence, tool='rnafold'):
    """
    Windowed structure prediction for medium-length sequences.
    Folds consecutive non-overlapping 200 nt windows and stitches the results.
    (Earlier versions folded every window twice and discarded the first pass.)

    If folding fails for ANY window the whole call fails with an error. A failed window is
    never replaced by 'all unpaired' dots, because that would invent single-stranded regions.
    Energy returned = sum of the window free energies.
    """
    final, energy, failed, total = _combine_windows(sequence, tool)
    if total == 0:
        return None, None, 'No windows could be processed'
    if failed:
        return None, None, 'Folding failed for {} of {} windows'.format(failed, total)
    return final, round(energy, 2), None


def _combine_windows(sequence, tool):
    """
    Fold non-overlapping WINDOW_SIZE windows and stitch the results.
    Returns (structure, summed energy, failed windows, total windows).
    """
    n      = len(sequence)
    result = []
    energy = 0.0
    failed = 0
    total  = 0

    for start in range(0, n, WINDOW_SIZE):
        end  = min(start + WINDOW_SIZE, n)
        wseq = sequence[start:end]

        if len(wseq) < 10:
            result.append('.' * len(wseq))      # fragment too short to fold (< 10 nt); cannot pair meaningfully
            continue

        total += 1
        if tool == 'rnastructure':
            ws, we, err = run_rnastructure_fold(wseq, timeout=60)
            if err or not ws:
                ws, we, err = run_rnafold(wseq, timeout=60)
        else:
            ws, we, err = run_rnafold(wseq, timeout=60)
            if err or not ws:
                ws, we, err = run_rnastructure_fold(wseq, timeout=60)

        if err or not ws or len(ws) != len(wseq):
            failed += 1
            result.append('.' * len(wseq))      # placeholder only; the caller rejects the whole result
            continue

        energy += we if we is not None else 0.0
        result.append(_balance_brackets(ws))

    combined = ''.join(result)
    if len(combined) < n:
        combined += '.' * (n - len(combined))
    return combined[:n], energy, failed, total


def _balance_brackets(struct):
    """
    Ensure dot-bracket string is properly balanced.
    Any unmatched bracket becomes a dot.
    """
    result = list(struct)
    stack  = []
    for i, ch in enumerate(result):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if stack:
                stack.pop()
            else:
                result[i] = '.'
    # Unmatched opens
    for i in stack:
        result[i] = '.'
    return ''.join(result)


# ─── WINDOWED PARTITION ────────────────────────────────────────────
def predict_windowed_partition(sequence, tool='rnafold'):
    """
    Windowed partition function for base pair probabilities.
    Used for medium-length sequences (1000-5000nt).

    Returns:
        (pair_probs dict, error)
        pair_probs: {(i,j): probability}
    """
    n         = len(sequence)
    all_probs = {}

    for start in range(0, n, STEP_SIZE):
        end    = min(start + WINDOW_SIZE, n)
        wseq   = sequence[start:end]

        if len(wseq) < 28:
            break

        if tool == 'rnastructure':
            probs, err = _rnastructure_partition_window(wseq)
        else:
            probs, err = _rnafold_partition_window(wseq)

        if err or not probs:
            continue

        # Map window coordinates to global coordinates
        for (i, j), p in probs.items():
            gi = i + start
            gj = j + start
            if gi < n and gj < n:
                # Average if position already has a probability
                if (gi, gj) in all_probs:
                    all_probs[(gi, gj)] = (all_probs[(gi, gj)] + p) / 2
                    all_probs[(gj, gi)] = all_probs[(gi, gj)]
                else:
                    all_probs[(gi, gj)] = p
                    all_probs[(gj, gi)] = p

    return all_probs, None


def _rnafold_partition_window(sequence):
    """Run RNAfold --partfunc on a window."""
    tmp_dir = tempfile.mkdtemp()
    try:
        sp = os.path.join(tmp_dir, 'seq.fa')
        with open(sp, 'w') as f:
            f.write(f'>seq\n{sequence}\n')

        r = subprocess.run(
            ['RNAfold', '--partfunc', '--noPS', sp],
            capture_output=True, text=True,
            timeout=120, cwd=tmp_dir
        )
        if r.returncode != 0:
            return {}, 'RNAfold partition error'

        # Find dot plot file
        dp_file = None
        for fn in os.listdir(tmp_dir):
            if fn.endswith('_dp.ps'):
                dp_file = os.path.join(tmp_dir, fn)
                break

        if not dp_file or not os.path.exists(dp_file):
            return {}, 'Dot plot not found'

        probs = _parse_dotplot(dp_file)
        return probs, None

    except Exception as e:
        return {}, str(e)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _parse_dotplot(dp_file):
    """Parse RNAfold dot plot PostScript file."""
    pp = {}
    try:
        with open(dp_file) as f:
            content = f.read()
        for m in re.findall(r'(\d+)\s+(\d+)\s+([0-9.]+)\s+ubox', content):
            i = int(m[0]) - 1
            j = int(m[1]) - 1
            p = float(m[2]) ** 2
            pp[(i, j)] = p
            pp[(j, i)] = p
    except Exception:
        pass
    return pp


def _rnastructure_partition_window(sequence):
    """Run RNAstructure partition + ProbabilityPlot on a window."""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.seq', delete=False) as f:
            dna = sequence.replace('U', 'T')
            f.write(f'; seq\nsequence\n{dna}\n1\n')
            seq_file = f.name

        pfs_file  = seq_file.replace('.seq', '.pfs')
        prob_file = seq_file.replace('.seq', '_probs.txt')

        r1 = subprocess.run(
            ['partition', seq_file, pfs_file],
            capture_output=True, text=True, timeout=120
        )
        if os.path.exists(seq_file):
            os.unlink(seq_file)

        if r1.returncode != 0:
            return {}, 'partition error'

        r2 = subprocess.run(
            ['ProbabilityPlot', pfs_file, prob_file, '--text'],
            capture_output=True, text=True, timeout=60
        )
        if os.path.exists(pfs_file):
            os.unlink(pfs_file)

        if r2.returncode != 0:
            return {}, 'ProbabilityPlot error'

        probs = _parse_rnastructure_probs(prob_file)
        if os.path.exists(prob_file):
            os.unlink(prob_file)
        return probs, None

    except Exception as e:
        return {}, str(e)


def _parse_rnastructure_probs(path):
    """Parse RNAstructure ProbabilityPlot text output."""
    pp = {}
    try:
        with open(path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        i = int(parts[0]) - 1
                        j = int(parts[1]) - 1
                        p = 10 ** (-float(parts[2]))
                        pp[(i, j)] = p
                        pp[(j, i)] = p
                    except Exception:
                        continue
    except Exception:
        pass
    return pp


# ─── LOCAL ACCESSIBILITY (RNAplfold) ───────────────────────────────
def run_rnaplfold_unpaired(sequence, window=None, span=None, timeout=600):
    """
    Per-nucleotide probability of being unpaired from RNAplfold (local folding, ViennaRNA).
    Works for any length in seconds-to-minutes, so seed confidence is available above 5000 nt too.
    Returns ({str(0-based index): probability}, error or None).
    Settings: window 200 / max base-pair span 100 (standard values for accessibility analysis).
    """
    window = window or PLFOLD_WINDOW
    span   = span   or PLFOLD_SPAN
    tmpdir = tempfile.mkdtemp(prefix='plfold_')
    try:
        with open(os.path.join(tmpdir, 'seq.fa'), 'w') as f:
            f.write('>seq\n' + sequence + '\n')
        with open(os.path.join(tmpdir, 'seq.fa')) as fin:
            result = subprocess.run(
                ['RNAplfold', '-W', str(window), '-L', str(span), '-u', '1'],
                stdin=fin, cwd=tmpdir, capture_output=True, text=True, timeout=timeout
            )
        out = os.path.join(tmpdir, 'seq_lunp')
        if not os.path.exists(out):
            return {}, 'RNAplfold produced no output: {}'.format(result.stderr[:150])
        probs = {}
        with open(out) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[1] != 'NA':
                    probs[str(int(parts[0]) - 1)] = round(max(0.0, min(1.0, float(parts[1]))), 4)
        if len(probs) < 0.9 * len(sequence):
            return {}, 'RNAplfold output incomplete ({} of {} positions)'.format(len(probs), len(sequence))
        return probs, None
    except FileNotFoundError:
        return {}, 'RNAplfold not installed'
    except subprocess.TimeoutExpired:
        return {}, 'RNAplfold timeout ({}s)'.format(timeout)
    except Exception as e:
        return {}, str(e)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
