"""
pipeline_utils.py
-----------------
Shared conversion functions used by 01_convert.py and 02_esm_fold.py.
Single source of truth for the lyric → AA mapping.
Do not modify the LETTER_TO_AA table or AA_BY_ABUNDANCE ordering —
any change here invalidates all existing bar_index_snapshot.json files.
"""

import math
import random
import re

# ---------------------------------------------------------------------------
# Frequency tables
# ---------------------------------------------------------------------------

AA_BY_ABUNDANCE = [
    ("L", 9.66), ("A", 8.25), ("G", 7.07), ("V", 6.87), ("E", 6.75),
    ("S", 6.56), ("I", 5.96), ("K", 5.84), ("R", 5.53), ("D", 5.45),
    ("T", 5.34), ("P", 4.70), ("N", 4.06), ("Q", 3.93), ("F", 3.86),
    ("H", 2.27), ("Y", 2.92), ("M", 2.42), ("C", 1.37), ("W", 1.08),
]

AA_BY_RARITY = list(reversed(AA_BY_ABUNDANCE))
CANONICAL_AAS = set(aa for aa, _ in AA_BY_ABUNDANCE)
ALL_AAS = [aa for aa, _ in AA_BY_ABUNDANCE]

LETTER_TO_AA = {
    "E": "L", "T": "A", "A": "G", "I": "V", "N": "E", "S": "S",
    "R": "I", "H": "T", "L": "P", "D": "K", "C": "R", "M": "D",
    "F": "F", "P": "Q", "G": "N", "W": "H", "Y": "C", "V": "W",
    "K": "Y", "Q": "M",
}

NON_STANDARD = {"Z": 0.09, "J": 0.16, "X": 0.23, "B": 1.48, "U": 2.73, "O": 7.64}
ENGLISH_MAX_PCT = 12.49


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def clean(bar: str) -> str:
    """Strip non-alpha chars, uppercase."""
    return re.sub(r"[^A-Za-z]", "", bar).upper()


def build_prob_vector(char_eng_pct: float, lambda_val: float) -> list:
    char_norm = 1.0 - (char_eng_pct / ENGLISH_MAX_PCT)
    n = len(AA_BY_RARITY)
    distances = [-lambda_val * abs(char_norm - (1.0 - i / (n - 1))) for i in range(n)]
    max_d = max(distances)
    exps = [math.exp(d - max_d) for d in distances]
    total = sum(exps)
    return [e / total for e in exps]


def sample_softmax(prob_vector: list, rng: random.Random) -> str:
    r = rng.random()
    cumulative = 0.0
    for i, p in enumerate(prob_vector):
        cumulative += p
        if r <= cumulative:
            return AA_BY_RARITY[i][0]
    return AA_BY_RARITY[-1][0]


# ---------------------------------------------------------------------------
# Converters — (cleaned_str, lambda, rng) -> (sequence, mask)
# ---------------------------------------------------------------------------

def convert_concordance(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask, cache = [], [], {}
    for letter in cleaned:
        if letter in LETTER_TO_AA:
            seq.append(LETTER_TO_AA[letter])
        elif letter in NON_STANDARD:
            if letter not in cache:
                cache[letter] = build_prob_vector(NON_STANDARD[letter], lambda_val)
            mask.append(len(seq))
            seq.append(sample_softmax(cache[letter], rng))
    return "".join(seq), mask


def convert_alanine(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask = [], []
    for letter in cleaned:
        if letter in LETTER_TO_AA:
            seq.append(LETTER_TO_AA[letter])
        elif letter in NON_STANDARD:
            mask.append(len(seq))
            seq.append("A")
    return "".join(seq), mask


def convert_random(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask = [], []
    for letter in cleaned:
        if letter in LETTER_TO_AA:
            seq.append(LETTER_TO_AA[letter])
        elif letter in NON_STANDARD:
            mask.append(len(seq))
            seq.append(rng.choice(ALL_AAS))
    return "".join(seq), mask


def convert_native(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask, cache = [], [], {}
    for letter in cleaned:
        if letter in CANONICAL_AAS:
            seq.append(letter)
        elif letter in NON_STANDARD:
            if letter not in cache:
                cache[letter] = build_prob_vector(NON_STANDARD[letter], lambda_val)
            mask.append(len(seq))
            seq.append(sample_softmax(cache[letter], rng))
    return "".join(seq), mask


def convert_native_alanine(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask = [], []
    for letter in cleaned:
        if letter in CANONICAL_AAS:
            seq.append(letter)
        elif letter in NON_STANDARD:
            mask.append(len(seq))
            seq.append("A")
    return "".join(seq), mask


CONVERTERS = {
    "concordance": convert_concordance,
    "alanine": convert_alanine,
    "random": convert_random,
    "native": convert_native,
    "native_alanine": convert_native_alanine,
}
