# XRF Data Analyzer — Development Status

GitHub: https://github.com/M-Song-ChE/XRF-analyzer  
Last updated: 2026-05-17

---

## Overview

A single-file tkinter GUI app that loads XRF spot-map CSV files, converts mass% to atomic fraction, and computes averaged composition with standard deviation for a user-selected subset of elements.

---

## Implemented Features

### Data Processing
- Load multiple CSV files (each file managed independently)
- Auto-detect CSV format by searching for the X, Y header row — tolerant of variable metadata rows at the top
- mass% → atomic fraction conversion per spot: `moles_i = mass%_i / atomic_mass_i`, normalized over selected elements only
- Section filtering:
  - Sections where the moles subtotal of **all selected elements is exactly 0.0** are marked "Undefined" and excluded from the average
  - Sections where at least one selected element is > 0 are always included (e.g. Pt=0 & Ni>0 → included)
- Mean and sample standard deviation (ddof=1) computed across all valid spots

### Periodic Table Panel
- Full 118-element interactive periodic table
- Color-coded by category (alkali, alkaline-earth, transition, post-transition, metalloid, nonmetal, halogen, noble-gas, lanthanide, actinide)
- Three button states: included (vivid color), detected-but-excluded (dim + groove), undetected (flat/inactive)
- Click to toggle elements in/out of composition calculation
- **Zoom in/out**: toolbar `−` / `+` buttons and mouse scroll wheel (range: 30%–300% of natural size)
- Zoom anchors to horizontal center of the table
- Auto-resizes when the pane is resized (zoom factor preserved)
- Default zoom: 50% of natural auto-fit size

### Per-Element Stats Table (middle panel)
- Shows all detected elements
- Columns: Element / Z / Mean at% (renorm) / ±σ / N spots / Mean mass%
- Included elements: green background + bold; excluded elements: gray text with values in parentheses
- Column header click to sort

### Per-File Composition Table (bottom panel)
- One row per loaded file — filename, composition (alloy notation), Total/Defined/Undefined section counts, per-element at% ± σ
- Alloy notation with Unicode subscript digits: `Pt₄₅.₃ Ni₅₄.₇`
- **Export CSV** button
- Click a file in the sidebar to view that file's stats alone; "View All Files Combined" to aggregate
- All columns uniform width with default stretch behavior

### Layout
- Left sidebar: file list, element selection controls, category legend
- Right area: 3-pane vertical split (periodic table / element stats / per-file table)
- Each pane resizable by dragging the sash
- Initial sash positions set automatically on startup (~1/4 for PT, ~55% for element table, remainder for per-file table)

---

## File Structure

```
XRF data extractor/
├── xrf_analyzer.py      # Main app (single file, ~1000 lines)
├── DEV_STATUS.md        # This file
├── PtNi_Rec_2.csv       # Sample data
├── .gitignore
└── .venv/               # Python virtual environment (numpy, tkinter)
```

---

## Calculation Summary

```
Per spot:
  moles_i  = mass%_i / atomic_mass_i   (selected elements only)
  subtotal = sum(moles_i)

  if subtotal == 0.0  →  Undefined (excluded)
  else                →  at%_i = moles_i / subtotal

Per file / all files:
  mean_at%_i = mean(at%_i across valid spots)
  σ_i        = std(at%_i across valid spots, ddof=1)
```

---

## Commit History

| Hash | Description |
|------|-------------|
| `fa3eba2` | Add DEV_STATUS.md |
| `dd5a438` | Remove alloy expression box; restore per-element table; halve PT default size |
| `99fd439` | Fix layout: cap PT pane height, center PT grid, per-file column stretch |
| `566f81e` | Add PT zoom controls, adjust table font sizes, uniform column widths |
| `fddcd55` | Initial commit |

---

## Known Issues / Future Ideas

- [ ] Composition column text may be clipped for long alloy strings (consider tooltip)
- [ ] Focus handling after file removal could be smoother
- [ ] Sort column name mismatch between heading label and internal column ID (minor)
- [ ] Future: spot-coordinate heatmap visualization
- [ ] Future: composition comparison chart across multiple files
