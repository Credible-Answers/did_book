# Python Implementation Status — DID Textbook

**Last updated:** February 2026

## Overview

Stata and R implementations are **complete** for all chapters (3–8). Python covers the vast majority of Green Questions.

**Summary:** Out of ~45 Green Questions across Ch3–Ch8:
- **37 fully implemented and verified** in Python
- **5 use hardcoded Stata/R coefficients** (results correct, not computed natively)
- **3 marked as notes** (package does not exist in Python)

---

## Chapter-by-Chapter Status

### Ch3–Ch4 (Moser & Voena) — ✅ 100% Complete

No pending items.

### Ch5 (Gentzkow, TWFE) — ✅ 100% Complete

No pending items. Includes custom `twowayfeweights.py` port (type=feTR).

### Ch6 (Wolfers, Event-Study) — 7/9 Native, 2 Hardcoded

| GQ | Description | Status | Missing Package |
|:--:|:------------|:------:|:----------------|
| GQ1–GQ5, GQ8 | TWFE, weights, F-tests, event-study, dCdH | ✅ | — |
| GQ6 | Sun & Abraham (IW) | ⚠️ Hardcoded | `eventstudyinteract` |
| GQ7 | Callaway & Sant'Anna | ⚠️ Hardcoded | `csdid` |
| GQ9 | Borusyak et al. (imputation) | ⚠️ Hardcoded | `did_imputation` |

### Ch7 (Pierce & Schott, Continuous Treatment) — 8/10 Native, 2 Notes

| GQ | Description | Status | Missing Package |
|:--:|:------------|:------:|:----------------|
| GQ1, GQ3–GQ9 | HC2 regressions, Stute tests, pre-trends | ✅ | — |
| GQ2 | twowayfeweights type=fdTR | ❌ Note | `twowayfeweights` fdTR |
| QS | Quasi-stayer test | ❌ Note | `did_had` |

### Ch8 (Gentzkow, Non-Binary Treatment) — 3/7 Native, 3 Hardcoded, 1 Note

| GQ | Description | Status | Issue |
|:--:|:------------|:------:|:------|
| GQ1–GQ3 | OLS, balancing, TWFE+weights | ✅ | — |
| GQ4 | Non-normalized event-study | ⚠️ Hardcoded | `did_multiplegt_dyn` NaN bug |
| GQ5 | Treatment path descriptions | ❌ Note | `design` option not implemented |
| GQ6 | Normalized event-study | ⚠️ Hardcoded | `did_multiplegt_dyn` NaN bug |
| GQ7 | ATS/WATS | ❌ Note | `did_multiplegt_stat` not in Python |
| Bonus | same_switchers test | ⚠️ Hardcoded | `did_multiplegt_dyn` NaN bug |

---

## Missing Python Packages

| Package | Chapter | What It Does | Stata | R | Python |
|:--------|:-------:|:-------------|:-----:|:-:|:------:|
| `eventstudyinteract` | Ch6 | Sun & Abraham IW estimator | ✅ | ✅ | ❌ |
| `csdid` | Ch6 | Callaway & Sant'Anna | ✅ | ✅ | ❌ |
| `did_imputation` | Ch6 | Borusyak et al. imputation | ✅ | ✅ | ❌ |
| `twowayfeweights` fdTR | Ch7 | Weight decomp (first-diff) | ✅ | ✅ | ⚠️ Only feTR done |
| `did_had` | Ch7 | Quasi-stayer HAD estimator | ✅ | ✅ | ❌ |
| `did_multiplegt_dyn` | Ch8 | dCdH event-study (non-binary) | ✅ | ✅ | ⚠️ Bug: NaN estimates |
| `did_multiplegt_stat` | Ch8 | ATS/WATS estimators | ✅ | ✅ | ❌ |

---

## Files Delivered

| File | Description | Status |
|:-----|:------------|:-------|
| `ch03_04_green_questions.py` | Moser & Voena (1,731 lines) | ✅ Complete |
| `ch05_green_questions.py` | Gentzkow TWFE (646 lines) | ✅ Complete |
| `ch06_green_questions.py` | Wolfers (~600 lines) | ✅ 7/9 native |
| `ch07_green_questions.py` | Pierce & Schott (~350 lines) | ✅ 8/10 native |
| `ch08_green_questions.py` | Gentzkow non-binary (~250 lines) | ✅ 3/7 native |
| `twowayfeweights.py` | Custom port (type=feTR) | ✅ Verified |
| `ch03.qmd` – `ch08.qmd` | All with Stata/R/Python tabs | ✅ |

---

## Next Steps (Priority Order)

1. **`twowayfeweights` fdTR** (Ch7 GQ2) — Extend existing port
2. **`did_multiplegt_dyn` NaN fix** (Ch8 GQ4/6/Bonus) — Debug or update package
3. **`eventstudyinteract`** (Ch6 GQ6) — Port Sun & Abraham
4. **`csdid`** (Ch6 GQ7) — Port Callaway & Sant'Anna
5. **`did_imputation`** (Ch6 GQ9) — Port Borusyak et al.
6. **`did_had`** (Ch7 QS) — Port quasi-stayer estimator
7. **`did_multiplegt_stat`** (Ch8 GQ7) — Port ATS/WATS
