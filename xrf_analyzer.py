#!/usr/bin/env python3
"""XRF Data Analyzer — mass fraction → atomic fraction with interactive periodic table."""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import csv
import os
from collections import defaultdict

import numpy as np

# ─── Subscript helper ─────────────────────────────────────────────────────────

_SUB = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')

def to_sub(s):
    """Convert digit characters in s to Unicode subscripts; keep '.' and others."""
    return s.translate(_SUB)

def alloy_notation(incl_elems, mean_af, decimals=1):
    """Return compact alloy string with subscript numbers, e.g. Pt₄₅.₃Ni₅₄.₇"""
    parts = []
    for e in incl_elems:
        if e in mean_af:
            num = f"{mean_af[e]*100:.{decimals}f}"
            parts.append(e + to_sub(num))
    return "".join(parts)

# ─── Element data ─────────────────────────────────────────────────────────────

ATOMIC_NUMBERS = {
    'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
    'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
    'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,
    'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,
    'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,
    'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,
    'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,
    'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,
    'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,
    'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,
    'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,
    'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,
    'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,'Rf':104,'Db':105,
    'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
    'Nh':113,'Fl':114,'Mc':115,'Lv':116,'Ts':117,'Og':118,
}

ATOMIC_MASSES = {
    'H':1.008,'He':4.003,'Li':6.941,'Be':9.012,'B':10.811,'C':12.011,
    'N':14.007,'O':15.999,'F':18.998,'Ne':20.180,'Na':22.990,'Mg':24.305,
    'Al':26.982,'Si':28.086,'P':30.974,'S':32.065,'Cl':35.453,'Ar':39.948,
    'K':39.098,'Ca':40.078,'Sc':44.956,'Ti':47.867,'V':50.942,'Cr':51.996,
    'Mn':54.938,'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.38,
    'Ga':69.723,'Ge':72.630,'As':74.922,'Se':78.971,'Br':79.904,'Kr':83.798,
    'Rb':85.468,'Sr':87.620,'Y':88.906,'Zr':91.224,'Nb':92.906,'Mo':95.960,
    'Tc':98.000,'Ru':101.07,'Rh':102.906,'Pd':106.42,'Ag':107.868,'Cd':112.411,
    'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.600,'I':126.904,'Xe':131.293,
    'Cs':132.905,'Ba':137.327,'La':138.905,'Ce':140.116,'Pr':140.908,'Nd':144.242,
    'Pm':145.000,'Sm':150.360,'Eu':151.964,'Gd':157.250,'Tb':158.925,'Dy':162.500,
    'Ho':164.930,'Er':167.259,'Tm':168.934,'Yb':173.045,'Lu':174.967,
    'Hf':178.490,'Ta':180.948,'W':183.840,'Re':186.207,'Os':190.230,'Ir':192.217,
    'Pt':195.084,'Au':196.967,'Hg':200.592,'Tl':204.380,'Pb':207.200,'Bi':208.980,
    'Po':209.000,'At':210.000,'Rn':222.000,'Fr':223.000,'Ra':226.000,'Ac':227.000,
    'Th':232.038,'Pa':231.036,'U':238.029,'Np':237.000,'Pu':244.000,'Am':243.000,
    'Cm':247.000,'Bk':247.000,'Cf':251.000,'Es':252.000,'Fm':257.000,'Md':258.000,
    'No':259.000,'Lr':262.000,'Rf':267.000,'Db':268.000,'Sg':271.000,'Bh':270.000,
    'Hs':277.000,'Mt':276.000,'Ds':281.000,'Rg':280.000,'Cn':285.000,
    'Nh':284.000,'Fl':289.000,'Mc':288.000,'Lv':293.000,'Ts':294.000,'Og':294.000,
}

# (symbol, display_row, display_col) — 18-column, rows 0-6: periods 1-7,
# row 7: spacer, rows 8-9: lanthanides/actinides
PERIODIC_TABLE_LAYOUT = [
    ('H',0,0),('He',0,17),
    ('Li',1,0),('Be',1,1),
    ('B',1,12),('C',1,13),('N',1,14),('O',1,15),('F',1,16),('Ne',1,17),
    ('Na',2,0),('Mg',2,1),
    ('Al',2,12),('Si',2,13),('P',2,14),('S',2,15),('Cl',2,16),('Ar',2,17),
    ('K',3,0),('Ca',3,1),
    ('Sc',3,2),('Ti',3,3),('V',3,4),('Cr',3,5),('Mn',3,6),
    ('Fe',3,7),('Co',3,8),('Ni',3,9),('Cu',3,10),('Zn',3,11),
    ('Ga',3,12),('Ge',3,13),('As',3,14),('Se',3,15),('Br',3,16),('Kr',3,17),
    ('Rb',4,0),('Sr',4,1),
    ('Y',4,2),('Zr',4,3),('Nb',4,4),('Mo',4,5),('Tc',4,6),
    ('Ru',4,7),('Rh',4,8),('Pd',4,9),('Ag',4,10),('Cd',4,11),
    ('In',4,12),('Sn',4,13),('Sb',4,14),('Te',4,15),('I',4,16),('Xe',4,17),
    ('Cs',5,0),('Ba',5,1),
    ('Hf',5,3),('Ta',5,4),('W',5,5),('Re',5,6),
    ('Os',5,7),('Ir',5,8),('Pt',5,9),('Au',5,10),('Hg',5,11),
    ('Tl',5,12),('Pb',5,13),('Bi',5,14),('Po',5,15),('At',5,16),('Rn',5,17),
    ('Fr',6,0),('Ra',6,1),
    ('Rf',6,3),('Db',6,4),('Sg',6,5),('Bh',6,6),
    ('Hs',6,7),('Mt',6,8),('Ds',6,9),('Rg',6,10),('Cn',6,11),
    ('Nh',6,12),('Fl',6,13),('Mc',6,14),('Lv',6,15),('Ts',6,16),('Og',6,17),
    ('La',8,2),('Ce',8,3),('Pr',8,4),('Nd',8,5),('Pm',8,6),
    ('Sm',8,7),('Eu',8,8),('Gd',8,9),('Tb',8,10),('Dy',8,11),
    ('Ho',8,12),('Er',8,13),('Tm',8,14),('Yb',8,15),('Lu',8,16),
    ('Ac',9,2),('Th',9,3),('Pa',9,4),('U',9,5),('Np',9,6),
    ('Pu',9,7),('Am',9,8),('Cm',9,9),('Bk',9,10),('Cf',9,11),
    ('Es',9,12),('Fm',9,13),('Md',9,14),('No',9,15),('Lr',9,16),
]

CATEGORY = {
    'H':'nonmetal','He':'noble',
    'Li':'alkali','Be':'alkaline','B':'metalloid','C':'nonmetal',
    'N':'nonmetal','O':'nonmetal','F':'halogen','Ne':'noble',
    'Na':'alkali','Mg':'alkaline','Al':'post_trans','Si':'metalloid',
    'P':'nonmetal','S':'nonmetal','Cl':'halogen','Ar':'noble',
    'K':'alkali','Ca':'alkaline',
    'Sc':'transition','Ti':'transition','V':'transition','Cr':'transition',
    'Mn':'transition','Fe':'transition','Co':'transition','Ni':'transition',
    'Cu':'transition','Zn':'transition',
    'Ga':'post_trans','Ge':'metalloid','As':'metalloid','Se':'nonmetal',
    'Br':'halogen','Kr':'noble',
    'Rb':'alkali','Sr':'alkaline',
    'Y':'transition','Zr':'transition','Nb':'transition','Mo':'transition',
    'Tc':'transition','Ru':'transition','Rh':'transition','Pd':'transition',
    'Ag':'transition','Cd':'transition',
    'In':'post_trans','Sn':'post_trans','Sb':'metalloid','Te':'metalloid',
    'I':'halogen','Xe':'noble',
    'Cs':'alkali','Ba':'alkaline',
    'La':'lanthanide','Ce':'lanthanide','Pr':'lanthanide','Nd':'lanthanide',
    'Pm':'lanthanide','Sm':'lanthanide','Eu':'lanthanide','Gd':'lanthanide',
    'Tb':'lanthanide','Dy':'lanthanide','Ho':'lanthanide','Er':'lanthanide',
    'Tm':'lanthanide','Yb':'lanthanide','Lu':'lanthanide',
    'Hf':'transition','Ta':'transition','W':'transition','Re':'transition',
    'Os':'transition','Ir':'transition','Pt':'transition','Au':'transition',
    'Hg':'transition','Tl':'post_trans','Pb':'post_trans','Bi':'post_trans',
    'Po':'metalloid','At':'halogen','Rn':'noble',
    'Fr':'alkali','Ra':'alkaline',
    'Ac':'actinide','Th':'actinide','Pa':'actinide','U':'actinide',
    'Np':'actinide','Pu':'actinide','Am':'actinide','Cm':'actinide',
    'Bk':'actinide','Cf':'actinide','Es':'actinide','Fm':'actinide',
    'Md':'actinide','No':'actinide','Lr':'actinide',
    'Rf':'transition','Db':'transition','Sg':'transition','Bh':'transition',
    'Hs':'transition','Mt':'transition','Ds':'transition','Rg':'transition',
    'Cn':'transition','Nh':'post_trans','Fl':'post_trans','Mc':'post_trans',
    'Lv':'post_trans','Ts':'halogen','Og':'noble',
}

DIM_COLOR = {
    'alkali':    '#ffd5d5', 'alkaline':  '#ffe8cc', 'transition': '#cce5ff',
    'post_trans':'#ccf5e5', 'metalloid': '#e5f5cc', 'nonmetal':   '#fffbcc',
    'halogen':   '#ffe0ee', 'noble':     '#ecdcff', 'lanthanide': '#fce0ff',
    'actinide':  '#ffd8d8',
}

ACTIVE_COLOR = {
    'alkali':    ('#ff5252', 'white'),
    'alkaline':  ('#ff8f00', 'black'),
    'transition':('#1565c0', 'white'),
    'post_trans':('#00897b', 'white'),
    'metalloid': ('#558b2f', 'white'),
    'nonmetal':  ('#f9a825', 'black'),
    'halogen':   ('#ad1457', 'white'),
    'noble':     ('#6a1b9a', 'white'),
    'lanthanide':('#6a1b9a', 'white'),
    'actinide':  ('#b71c1c', 'white'),
}

THRESHOLD = 0.001   # 0.1 at% to count as "detected"


# ─── CSV parsing ──────────────────────────────────────────────────────────────

def parse_csv(path):
    """Return (elements, data_rows) where data_rows = [(x, y, [mass_fracs...]), ...]."""
    with open(path, newline='', encoding='utf-8-sig') as f:
        rows = list(csv.reader(f))

    hdr_idx = None
    for i, row in enumerate(rows):
        s = [c.strip() for c in row]
        if len(s) >= 2 and s[0] == 'X' and s[1] == 'Y':
            hdr_idx = i
            break
    if hdr_idx is None:
        raise ValueError(f"Cannot find 'X, Y' header row in {os.path.basename(path)}")

    elements = [c.strip() for c in rows[hdr_idx][2:] if c.strip()]
    data = []
    for row in rows[hdr_idx + 2:]:
        if not row or len(row) < 2:
            continue
        try:
            x, y = float(row[0].strip()), float(row[1].strip())
            mf = [float(row[i+2].strip()) if i+2 < len(row) else 0.0
                  for i in range(len(elements))]
            data.append((x, y, mf))
        except (ValueError, IndexError):
            continue
    return elements, data


# ─── Computation ──────────────────────────────────────────────────────────────

def mass_to_atomic(elements, mass_fracs):
    """Convert mass% list → atomic fraction list (sum=1)."""
    moles = [mf / ATOMIC_MASSES.get(e, 0.0) if ATOMIC_MASSES.get(e, 0.0) > 0 else 0.0
             for e, mf in zip(elements, mass_fracs)]
    total = sum(moles)
    if total <= 0:
        return [0.0] * len(elements)
    return [m / total for m in moles]


def _all_elements_in(files_subset):
    """Return list of all element symbols sorted by atomic number."""
    return sorted(
        {e for elems, _ in files_subset.values() for e in elems},
        key=lambda e: ATOMIC_NUMBERS.get(e, 999)
    )


def compute_stats_for(files_subset, included=None):
    """
    Compute mean ± std of atomic fractions for a given subset of files.

    If included is None:
        Use all elements, no renormalization, no section filtering.
    If included is a set:
        For each section, compute moles only for included elements.
        - subtotal (sum of moles) > 0  → valid section: renormalize and include.
        - subtotal == 0 (all selected elements absent) → mathematically undefined
          (0/0), excluded from mean/std; counted separately as n_undefined.
        Pt=0 Ni>0 or Pt>0 Ni=0 sections ARE included (composition is defined).

    Returns: (sorted_elements, mean_af, std_af, n_valid, n_total)
        n_total   = total sections in files_subset
        n_valid   = sections where composition was defined (subtotal > 0)
        n_total - n_valid = "no metal" / undefined sections
    """
    all_elems = _all_elements_in(files_subset)

    if included is not None:
        work_elems = sorted(included, key=lambda e: ATOMIC_NUMBERS.get(e, 999))
    else:
        work_elems = all_elems

    bucket = defaultdict(list)
    n_total = 0
    n_valid = 0

    for elements, rows in files_subset.values():
        for _, _, mf in rows:
            n_total += 1
            # Compute moles directly from mass% (not from the all-element
            # normalised atomic fraction) so the subtotal check reflects
            # actual Pt+Ni mole content.
            moles_sel = {}
            for e in work_elems:
                idx = elements.index(e) if e in elements else -1
                raw_mass = mf[idx] if idx >= 0 else 0.0
                am = ATOMIC_MASSES.get(e, 0.0)
                moles_sel[e] = raw_mass / am if am > 0 else 0.0

            subtotal = sum(moles_sel.values())

            if included is not None:
                # Strictly exclude only sections where ALL selected elements
                # have zero mass% (subtotal exactly 0).  A tiny but non-zero
                # value means some metal is present and composition is defined.
                if subtotal == 0.0:
                    continue  # undefined: 0/0 → skip, counted as n_total - n_valid
                n_valid += 1
                for e in work_elems:
                    bucket[e].append(moles_sel[e] / subtotal)
            else:
                n_valid += 1
                # For the raw (non-filtered) case use the all-element
                # normalised atomic fraction as before.
                af_all = mass_to_atomic(elements, mf)
                af_dict = dict(zip(elements, af_all))
                for e in work_elems:
                    bucket[e].append(af_dict.get(e, 0.0))

    mean_af, std_af = {}, {}
    for e in work_elems:
        arr = np.array(bucket[e]) if bucket[e] else np.array([0.0])
        mean_af[e] = float(arr.mean())
        std_af[e] = float(arr.std(ddof=1)) if len(arr) > 1 else 0.0

    return work_elems, mean_af, std_af, n_valid, n_total


def detected_elements(files_subset):
    """Return set of elements with mean at% > THRESHOLD in the given file subset."""
    if not files_subset:
        return set()
    elems, mean_af, _, _, _ = compute_stats_for(files_subset)
    return {e for e in elems if mean_af[e] > THRESHOLD}


# ─── Application ──────────────────────────────────────────────────────────────

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("XRF Data Analyzer")
        self.configure(bg='#eceff1')

        self.files_data = {}      # path -> (elements, rows)
        self.file_order = []      # insertion-order list of paths

        # Which file is "focused" in the listbox (None = show all files combined)
        self.focused_path = None

        # Set of element symbols the user has toggled ON for composition calculation
        # Starts as all detected elements; user clicks to toggle
        self.included = set()

        self._btn_widgets = {}    # symbol -> (frame, an_label, sym_label, category)
        self._sort_col = None
        self._sort_rev = False

        self._build_ui()
        self.minsize(1250, 820)
        self.geometry("1450x900")

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        self._setup_styles()

        # Horizontal pane: sidebar | right content (user can drag to resize)
        h_pane = ttk.PanedWindow(self, orient='horizontal')
        h_pane.pack(fill='both', expand=True, padx=8, pady=8)

        sidebar = ttk.Frame(h_pane, style='Side.TFrame', width=270)
        sidebar.pack_propagate(False)
        h_pane.add(sidebar, weight=0)
        self._build_sidebar(sidebar)

        right = ttk.Frame(h_pane, style='Outer.TFrame')
        h_pane.add(right, weight=1)
        self._build_right(right)

    def _setup_styles(self):
        s = ttk.Style()
        s.theme_use('clam')
        s.configure('Outer.TFrame', background='#eceff1')
        s.configure('Side.TFrame',  background='#ffffff')
        s.configure('TLabel',       background='#eceff1', font=('Segoe UI', 11))
        s.configure('SideH.TLabel', background='#ffffff', font=('Segoe UI', 12, 'bold'),
                    foreground='#1565c0')
        s.configure('TButton',      font=('Segoe UI', 11), padding=5)
        s.configure('Treeview',     font=('Segoe UI', 22), rowheight=52)
        s.configure('Treeview.Heading', font=('Segoe UI', 22, 'bold'))

    def _build_sidebar(self, parent):
        tk.Label(parent, text="XRF Data Analyzer", font=('Segoe UI', 14, 'bold'),
                 fg='#1565c0', bg='white').pack(anchor='w', padx=12, pady=(12, 2))
        tk.Label(parent, text="mass% → atomic fraction", font=('Segoe UI', 10),
                 fg='#607d8b', bg='white').pack(anchor='w', padx=12, pady=(0, 8))

        ttk.Separator(parent).pack(fill='x', padx=8)

        # ── File list ──
        ttk.Label(parent, text="Samples (CSV files)", style='SideH.TLabel').pack(
            anchor='w', padx=12, pady=(8, 2))

        lb_frame = tk.Frame(parent, bg='white')
        lb_frame.pack(fill='both', expand=True, padx=12, pady=(0, 4))
        vsb = ttk.Scrollbar(lb_frame)
        vsb.pack(side='right', fill='y')
        self.lb = tk.Listbox(lb_frame, yscrollcommand=vsb.set, selectmode='single',
                             font=('Segoe UI', 11), height=7, activestyle='none',
                             bg='#f5f5f5', relief='solid', borderwidth=1,
                             highlightthickness=0)
        self.lb.pack(side='left', fill='both', expand=True)
        vsb.config(command=self.lb.yview)
        self.lb.bind('<<ListboxSelect>>', self._on_file_select)

        file_btn = tk.Frame(parent, bg='white')
        file_btn.pack(fill='x', padx=12, pady=(0, 4))
        ttk.Button(file_btn, text="Add Files",
                   command=self._add_files).pack(side='left', expand=True, fill='x', padx=(0, 2))
        ttk.Button(file_btn, text="Remove",
                   command=self._remove_file).pack(side='left', expand=True, fill='x', padx=(2, 0))

        # "View all" button to deselect focus
        self.view_all_btn = ttk.Button(parent, text="View All Files Combined",
                                        command=self._view_all)
        self.view_all_btn.pack(fill='x', padx=12, pady=(0, 6))

        ttk.Separator(parent).pack(fill='x', padx=8)

        # ── Status ──
        ttk.Label(parent, text="Current View", style='SideH.TLabel').pack(
            anchor='w', padx=12, pady=(8, 2))
        self.view_var = tk.StringVar(value="—")
        tk.Label(parent, textvariable=self.view_var, font=('Segoe UI', 11),
                 bg='white', fg='#e65100', wraplength=240,
                 justify='left').pack(anchor='w', padx=12)

        self.summary_var = tk.StringVar(value="No data loaded.")
        tk.Label(parent, textvariable=self.summary_var, font=('Segoe UI', 11),
                 bg='white', fg='#37474f', wraplength=240,
                 justify='left').pack(anchor='w', padx=12, pady=(4, 0))

        ttk.Separator(parent).pack(fill='x', padx=8, pady=(8, 4))

        # ── Element selection controls ──
        ttk.Label(parent, text="Element Selection", style='SideH.TLabel').pack(
            anchor='w', padx=12, pady=(4, 2))
        tk.Label(parent, text="Click elements on the table to\ntoggle them in/out of the\ncomposition calculation.",
                 font=('Segoe UI', 10), bg='white', fg='#607d8b',
                 justify='left').pack(anchor='w', padx=12)

        sel_btn = tk.Frame(parent, bg='white')
        sel_btn.pack(fill='x', padx=12, pady=6)
        ttk.Button(sel_btn, text="Select All",
                   command=self._select_all).pack(side='left', expand=True, fill='x', padx=(0, 2))
        ttk.Button(sel_btn, text="Clear All",
                   command=self._clear_all).pack(side='left', expand=True, fill='x', padx=(2, 0))

        self.sel_status_var = tk.StringVar(value="")
        tk.Label(parent, textvariable=self.sel_status_var, font=('Segoe UI', 10),
                 bg='white', fg='#1565c0', wraplength=240,
                 justify='left').pack(anchor='w', padx=12)

        ttk.Separator(parent).pack(fill='x', padx=8, pady=(8, 4))

        # ── Category legend ──
        ttk.Label(parent, text="Category Colors", style='SideH.TLabel').pack(
            anchor='w', padx=12, pady=(4, 2))
        for label, cat in [
            ('Alkali metals',    'alkali'),
            ('Alkaline earth',   'alkaline'),
            ('Transition metals','transition'),
            ('Post-transition',  'post_trans'),
            ('Metalloids',       'metalloid'),
            ('Nonmetals',        'nonmetal'),
            ('Halogens',         'halogen'),
            ('Noble gases',      'noble'),
            ('Lanthanides',      'lanthanide'),
            ('Actinides',        'actinide'),
        ]:
            bg, _ = ACTIVE_COLOR.get(cat, ('#aaaaaa', 'white'))
            f = tk.Frame(parent, bg='white')
            f.pack(anchor='w', padx=12, pady=1, fill='x')
            tk.Label(f, text="  ", bg=bg, width=3, relief='raised').pack(side='left', padx=(0, 6))
            tk.Label(f, text=label, bg='white', font=('Segoe UI', 10),
                     fg='#37474f').pack(side='left')

    def _build_right(self, parent):
        # Vertical pane: periodic table | main results | per-file table
        self._v_pane = ttk.PanedWindow(parent, orient='vertical')
        self._v_pane.pack(fill='both', expand=True)

        pt_frame = ttk.Frame(self._v_pane, style='Outer.TFrame')
        self._v_pane.add(pt_frame, weight=1)
        self._build_periodic_table(pt_frame)

        mid_frame = ttk.Frame(self._v_pane, style='Outer.TFrame')
        self._v_pane.add(mid_frame, weight=3)
        self._build_mid_results(mid_frame)

        bot_frame = ttk.Frame(self._v_pane, style='Outer.TFrame')
        self._v_pane.add(bot_frame, weight=2)
        self._build_file_panel(bot_frame)

        # Set initial sash positions once the window is fully drawn
        self.after(200, self._set_initial_sash)

    def _set_initial_sash(self):
        """Pin sash positions so the PT pane doesn't swallow the tables."""
        total = self._v_pane.winfo_height()
        if total < 200:
            self.after(100, self._set_initial_sash)
            return
        pt_h  = min(260, total // 4)
        mid_h = int((total - pt_h) * 0.60)
        self._v_pane.sashpos(0, pt_h)
        self._v_pane.sashpos(1, pt_h + mid_h)

    def _build_periodic_table(self, parent):
        self._pt_outer  = parent        # used by resize handler
        self._pt_BW     = 56            # current cell width
        self._pt_BH     = 44            # current cell height
        self._pt_PAD    = 1
        self._pt_zoom   = 1.0           # manual zoom multiplier
        self._pt_all_frames = []        # (frame, 'cell'|'spacer') for bulk resize

        # Zoom toolbar
        zoom_bar = tk.Frame(parent, bg='#eceff1')
        zoom_bar.pack(fill='x', padx=4, pady=(2, 0))
        tk.Label(zoom_bar, text="Periodic Table", font=('Segoe UI', 10, 'bold'),
                 bg='#eceff1', fg='#37474f').pack(side='left', padx=(2, 6))
        tk.Button(zoom_bar, text="−", font=('Segoe UI', 11, 'bold'),
                  command=lambda: self._pt_zoom_by(-0.15), relief='flat',
                  bg='#cfd8dc', width=2, cursor='hand2').pack(side='left')
        tk.Button(zoom_bar, text="+", font=('Segoe UI', 11, 'bold'),
                  command=lambda: self._pt_zoom_by(+0.15), relief='flat',
                  bg='#cfd8dc', width=2, cursor='hand2').pack(side='left', padx=(2, 0))

        # Wrapper fills the pane; _pt_frame is centered inside it
        self._pt_wrapper = tk.Frame(parent, bg='#eceff1')
        self._pt_wrapper.pack(fill='both', expand=True, padx=2, pady=2)

        self._pt_frame = tk.Frame(self._pt_wrapper, bg='#eceff1')

        BW, BH, PAD = self._pt_BW, self._pt_BH, self._pt_PAD

        # Lanthanide / Actinide placeholder cells
        self._pt_placeholders = []
        for row, col, txt in [(5, 2, '57–71'), (6, 2, '89–103')]:
            f = tk.Frame(self._pt_frame, width=BW, height=BH,
                         bg='#e0e0e0', relief='flat', bd=1)
            f.grid(row=row, column=col, padx=PAD, pady=PAD, sticky='nsew')
            f.grid_propagate(False)
            lbl = tk.Label(f, text=txt, font=('Segoe UI', 7), bg='#e0e0e0', fg='#9e9e9e')
            lbl.place(relx=0.5, rely=0.5, anchor='center')
            self._pt_placeholders.append(f)
            self._pt_all_frames.append((f, 'cell'))

        # Row 7 spacer
        self._pt_spacers = []
        for col in range(18):
            f = tk.Frame(self._pt_frame, width=BW, height=8, bg='#eceff1')
            f.grid(row=7, column=col, padx=PAD)
            self._pt_spacers.append(f)
            self._pt_all_frames.append((f, 'spacer'))

        # Series labels
        for row, txt in [(8, 'Lantha-\nnides'), (9, 'Acti-\nnides')]:
            tk.Label(self._pt_frame, text=txt, font=('Segoe UI', 8),
                     bg='#eceff1', fg='#78909c', justify='right',
                     anchor='e').grid(row=row, column=0, columnspan=2,
                                      sticky='e', padx=(0, 4))

        # Element buttons
        for sym, disp_row, disp_col in PERIODIC_TABLE_LAYOUT:
            cat    = CATEGORY.get(sym, 'transition')
            an     = ATOMIC_NUMBERS.get(sym, '')
            dim_bg = DIM_COLOR.get(cat, '#e0e0e0')

            fr = tk.Frame(self._pt_frame, width=BW, height=BH,
                          bg=dim_bg, relief='flat', bd=1, cursor='hand2')
            fr.grid(row=disp_row, column=disp_col, padx=PAD, pady=PAD, sticky='nsew')
            fr.grid_propagate(False)

            an_lbl = tk.Label(fr, text=str(an), font=('Segoe UI', 7),
                               bg=dim_bg, fg='#bdbdbd')
            an_lbl.place(x=2, y=1)

            sym_lbl = tk.Label(fr, text=sym, font=('Segoe UI', 11, 'bold'),
                                bg=dim_bg, fg='#bdbdbd')
            sym_lbl.place(relx=0.5, rely=0.62, anchor='center')

            for w in (fr, an_lbl, sym_lbl):
                w.bind('<Button-1>', lambda e, s=sym: self._on_element_click(s))
                w.bind('<Enter>',    lambda e, s=sym: self._on_hover(s, True))
                w.bind('<Leave>',    lambda e, s=sym: self._on_hover(s, False))

            self._btn_widgets[sym] = (fr, an_lbl, sym_lbl, cat)
            self._pt_all_frames.append((fr, 'cell'))

        for c in range(18):
            self._pt_frame.grid_columnconfigure(c, minsize=BW + 2 * PAD)

        # Initial centered placement; kept in sync by _apply_pt_size
        self._pt_frame.place(relx=0.5, y=0, anchor='n')

        # Bind resize on wrapper; scroll-wheel zoom on grid
        self._pt_wrapper.bind('<Configure>', self._on_pt_resize, add='+')
        self._pt_frame.bind('<MouseWheel>', self._on_pt_scroll, add='+')

    def _on_pt_resize(self, event):
        """Recompute cell size whenever the periodic-table pane is resized."""
        if event.widget is not self._pt_wrapper:
            return
        avail_w = max(event.width - 8, 18 * 40)
        base_BW = max(40, avail_w // 19)
        BW = max(28, int(base_BW * self._pt_zoom))
        BH = max(22, int(BW * 44 / 56))
        if abs(BW - self._pt_BW) < 2:
            return
        self._pt_BW, self._pt_BH = BW, BH
        self._apply_pt_size(BW, BH)

    def _on_pt_scroll(self, event):
        """Scroll-wheel zoom on the periodic table."""
        self._pt_zoom_by(0.1 if event.delta > 0 else -0.1)

    def _pt_zoom_by(self, delta):
        """Adjust zoom factor and redraw periodic table cells."""
        self._pt_zoom = max(0.3, min(3.0, self._pt_zoom + delta))
        avail_w = max(self._pt_wrapper.winfo_width() - 8, 18 * 40)
        base_BW = max(40, avail_w // 19)
        BW = max(28, int(base_BW * self._pt_zoom))
        BH = max(22, int(BW * 44 / 56))
        self._pt_BW, self._pt_BH = BW, BH
        self._apply_pt_size(BW, BH)

    def _apply_pt_size(self, BW, BH):
        PAD = self._pt_PAD
        sym_fs = max(8,  BW // 5)
        an_fs  = max(6,  BW // 9)

        for c in range(18):
            self._pt_frame.grid_columnconfigure(c, minsize=BW + 2 * PAD)

        for f, kind in self._pt_all_frames:
            if kind == 'spacer':
                f.configure(width=BW)
            else:
                f.configure(width=BW, height=BH)

        for sym, (fr, an_lbl, sym_lbl, cat) in self._btn_widgets.items():
            fr.configure(width=BW, height=BH)
            an_lbl.configure(font=('Segoe UI', an_fs))
            sym_lbl.configure(font=('Segoe UI', sym_fs, 'bold'))

        # Re-center grid after size change
        self._pt_frame.place(relx=0.5, y=0, anchor='n')

    def _build_mid_results(self, parent):
        """Alloy expression box + per-element stats table."""
        hdr = tk.Frame(parent, bg='#eceff1')
        hdr.pack(fill='x', padx=4, pady=(4, 2))
        self.result_title = tk.Label(hdr,
            text="Load CSV files, then select elements to compute composition.",
            font=('Segoe UI', 10, 'bold'), fg='#1565c0', bg='#eceff1', anchor='w')
        self.result_title.pack(side='left', anchor='w')

        # Alloy expression: symbol bold-large, subscript numbers below
        alloy_outer = tk.Frame(parent, bg='#1565c0', bd=2, relief='flat')
        alloy_outer.pack(fill='x', padx=4, pady=(0, 4))
        self.alloy_text = tk.Text(
            alloy_outer, height=2, font=('Segoe UI', 13),
            fg='#1a237e', bg='#e8eaf6', relief='flat',
            wrap='word', padx=10, pady=6, state='disabled', cursor='arrow')
        self.alloy_text.pack(fill='x')
        # 'sym' = element symbol (bold large); 'sub' = subscript number (smaller, low)
        self.alloy_text.tag_configure('sym',
            font=('Segoe UI', 15, 'bold'), foreground='#1a237e')
        self.alloy_text.tag_configure('sub',
            font=('Segoe UI', 10), foreground='#b71c1c', offset=-3)
        self.alloy_text.tag_configure('err',
            font=('Segoe UI', 9), foreground='#9e9e9e', offset=0)

        # Per-element stats table
        tbl = tk.Frame(parent, bg='#eceff1')
        tbl.pack(fill='both', expand=True, padx=4, pady=(0, 4))

        cols = ('Element', 'Z', 'Mean at% (renorm)', '± σ', 'N spots', 'Mean mass%')
        self.tree = ttk.Treeview(tbl, columns=cols, show='headings')
        for col, w, anc in zip(cols,
                               [130, 80, 220, 150, 130, 200],
                               ['center','center','center','center','center','center']):
            self.tree.heading(col, text=col, command=lambda c=col: self._sort_tree(c))
            self.tree.column(col, width=w, anchor=anc)

        vsb2 = ttk.Scrollbar(tbl, orient='vertical', command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb2.set)
        self.tree.pack(side='left', fill='both', expand=True)
        vsb2.pack(side='right', fill='y')
        self.tree.tag_configure('incl', background='#e8f5e9',
                                font=('Segoe UI', 22, 'bold'))
        self.tree.tag_configure('excl', foreground='#bdbdbd')

    def _build_file_panel(self, parent):
        """Per-file breakdown table + export button (in its own resizable pane)."""
        file_hdr = tk.Frame(parent, bg='#eceff1')
        file_hdr.pack(fill='x', padx=4, pady=(4, 2))
        tk.Label(file_hdr,
                 text="Per-file composition (selected elements, renormalized)",
                 font=('Segoe UI', 14, 'bold'), fg='#1565c0',
                 bg='#eceff1').pack(side='left', anchor='w')
        ttk.Button(file_hdr, text="Export CSV",
                   command=self._export_file_table).pack(side='right', padx=(4, 0))

        self._file_tbl_frame = tk.Frame(parent, bg='#eceff1')
        self._file_tbl_frame.pack(fill='both', expand=True, padx=4, pady=(0, 4))
        self.file_tree = None
        self._file_table_data = []

    # ── File management ───────────────────────────────────────────────────────

    def _add_files(self):
        paths = filedialog.askopenfilenames(
            title="Select XRF CSV Files",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])
        added = 0
        for p in paths:
            if p in self.files_data:
                continue
            try:
                elements, rows = parse_csv(p)
                self.files_data[p] = (elements, rows)
                self.file_order.append(p)
                self.lb.insert(tk.END, os.path.basename(p))
                added += 1
            except Exception as ex:
                messagebox.showerror("Load Error",
                    f"Could not load {os.path.basename(p)}:\n{ex}")
        if added:
            self._on_data_changed()

    def _remove_file(self):
        sel = self.lb.curselection()
        if not sel:
            return
        idx = sel[0]
        name = self.lb.get(idx)
        path = next((p for p in self.files_data if os.path.basename(p) == name), None)
        if path:
            del self.files_data[path]
            self.file_order.remove(path)
            if self.focused_path == path:
                self.focused_path = None
        self.lb.delete(idx)
        self._on_data_changed()

    def _on_file_select(self, event=None):
        sel = self.lb.curselection()
        if not sel:
            return
        idx = sel[0]
        name = self.lb.get(idx)
        path = next((p for p in self.files_data if os.path.basename(p) == name), None)
        if path:
            self.focused_path = path
            self._on_data_changed(reset_included=False)

    def _view_all(self):
        self.lb.selection_clear(0, tk.END)
        self.focused_path = None
        self._on_data_changed(reset_included=False)

    # ── Element selection ─────────────────────────────────────────────────────

    def _on_element_click(self, sym):
        active = self._current_active()
        if sym not in active:
            return
        if sym in self.included:
            self.included.discard(sym)
        else:
            self.included.add(sym)
        self._refresh_all_buttons()
        self._update_results()
        self._update_sel_status()

    def _select_all(self):
        self.included = set(self._current_active())
        self._refresh_all_buttons()
        self._update_results()
        self._update_sel_status()

    def _clear_all(self):
        self.included.clear()
        self._refresh_all_buttons()
        self._update_results()
        self._update_sel_status()

    # ── Data change pipeline ──────────────────────────────────────────────────

    def _current_files(self):
        """Return the file subset to compute stats on (focused or all)."""
        if self.focused_path and self.focused_path in self.files_data:
            return {self.focused_path: self.files_data[self.focused_path]}
        return self.files_data

    def _current_active(self):
        """Set of elements detected (mean > THRESHOLD) in the current file subset."""
        return detected_elements(self._current_files())

    def _on_data_changed(self, reset_included=True):
        if not self.files_data:
            self.focused_path = None
            self.included.clear()
            self._refresh_all_buttons()
            self._update_results()
            self.summary_var.set("No data loaded.")
            self.view_var.set("—")
            return

        active = self._current_active()

        if reset_included:
            # Auto-select all detected elements on fresh load
            self.included = set(active)
        else:
            # Keep existing selection, but drop elements not in this view's active set
            self.included &= active
            # Also add newly detected elements that weren't previously visible
            # (so switching files doesn't leave the user with nothing selected)
            if not self.included:
                self.included = set(active)

        # Update view label
        if self.focused_path:
            self.view_var.set(f"Sample: {os.path.basename(self.focused_path)}")
        else:
            self.view_var.set(f"All {len(self.files_data)} file(s) combined")

        # Update summary
        elems_sorted = sorted(active, key=lambda e: ATOMIC_NUMBERS.get(e, 999))
        subset = self._current_files()
        n_spots = sum(len(rows) for _, rows in subset.values())
        self.summary_var.set(
            f"Spots: {n_spots}\n"
            f"Detected: {len(active)} element(s)\n"
            f"{', '.join(elems_sorted)}"
        )

        self._refresh_all_buttons()
        self._update_results()
        self._update_sel_status()

    def _update_sel_status(self):
        incl = sorted(self.included, key=lambda e: ATOMIC_NUMBERS.get(e, 999))
        if incl:
            self.sel_status_var.set(f"Selected ({len(incl)}): {', '.join(incl)}")
        else:
            self.sel_status_var.set("No elements selected.")

    # ── Button rendering ──────────────────────────────────────────────────────

    def _apply_button_state(self, sym):
        fr, an_lbl, sym_lbl, cat = self._btn_widgets[sym]
        active = self._current_active()

        if sym in self.included:
            # Bright vivid — included in calculation
            bg, fg = ACTIVE_COLOR.get(cat, ('#555', 'white'))
            relief = 'raised'
            an_fg  = fg
        elif sym in active:
            # Detected but excluded — use dim color, slightly darker text
            bg     = DIM_COLOR.get(cat, '#e0e0e0')
            fg     = '#888888'
            an_fg  = '#aaaaaa'
            relief = 'groove'
        else:
            # Not detected at all
            bg     = DIM_COLOR.get(cat, '#e0e0e0')
            fg     = '#c0c0c0'
            an_fg  = '#c0c0c0'
            relief = 'flat'

        fr.configure(bg=bg, relief=relief)
        an_lbl.configure(bg=bg, fg=an_fg)
        sym_lbl.configure(bg=bg, fg=fg)

    def _refresh_all_buttons(self):
        for sym in self._btn_widgets:
            self._apply_button_state(sym)

    def _on_hover(self, sym, entering):
        active = self._current_active()
        if sym not in active:
            return
        fr, an_lbl, sym_lbl, cat = self._btn_widgets[sym]
        if entering:
            base_bg, fg = ACTIVE_COLOR.get(cat, ('#777', 'white'))
            hover_bg = self._lighten(base_bg)
            fr.configure(bg=hover_bg)
            an_lbl.configure(bg=hover_bg)
            sym_lbl.configure(bg=hover_bg)
        else:
            self._apply_button_state(sym)

    @staticmethod
    def _lighten(hex_color):
        r = min(255, int(hex_color[1:3], 16) + 40)
        g = min(255, int(hex_color[3:5], 16) + 40)
        b = min(255, int(hex_color[5:7], 16) + 40)
        return f'#{r:02x}{g:02x}{b:02x}'

    # ── Results panel ─────────────────────────────────────────────────────────

    def _update_results(self):
        for row in self.tree.get_children():
            self.tree.delete(row)
        self.alloy_text.configure(state='normal')
        self.alloy_text.delete('1.0', 'end')
        self.alloy_text.configure(state='disabled')

        subset = self._current_files()
        if not subset:
            self.result_title.configure(
                text="Load CSV files, then select elements to compute composition.")
            self._rebuild_file_table([])
            return

        active = self._current_active()

        if not self.included:
            self.result_title.configure(
                text="No elements selected — click elements on the periodic table.")
            self._rebuild_file_table([])
            return

        # Stats for included elements only — sections where all selected
        # elements are 0 are excluded (undefined composition).
        incl_elems, mean_af, std_af, n_valid, n_total = compute_stats_for(
            subset, included=self.included)

        # Raw stats for ALL detected elements (used for excluded-element rows)
        all_elems, mean_af_raw, std_af_raw, _, _ = compute_stats_for(subset)

        # Mean mass% from raw data
        mass_bucket = defaultdict(list)
        for elems_f, rows in subset.values():
            for _, _, mf in rows:
                mf_dict = dict(zip(elems_f, mf))
                for e in active:
                    mass_bucket[e].append(mf_dict.get(e, 0.0))

        # Title — show valid / total section counts
        view_label = (f"Sample: {os.path.basename(self.focused_path)}"
                      if self.focused_path else
                      f"All {len(self.files_data)} file(s) combined")
        n_undef = n_total - n_valid
        undef_str = f"  |  {n_undef} section(s) undefined (Pt=Ni=0)" if n_undef else ""
        self.result_title.configure(
            text=(f"{view_label}  |  {len(self.included)} element(s) selected  |  "
                  f"{n_valid}/{n_total} valid sections{undef_str}"))

        # ── Alloy expression — e.g.  Pt₄₅.₃(±2.1) Ni₅₄.₇(±1.8) ──
        self.alloy_text.configure(state='normal')
        for e in incl_elems:
            pct   = mean_af[e] * 100
            sigma = std_af[e] * 100
            self.alloy_text.insert('end', e, 'sym')
            self.alloy_text.insert('end', to_sub(f"{pct:.1f}"), 'sub')
            self.alloy_text.insert('end', f"(±{sigma:.1f}) ", 'err')
        self.alloy_text.configure(state='disabled')

        # ── Per-element table ──
        for e in all_elems:
            if mean_af_raw.get(e, 0) <= 0 and e not in active:
                continue
            is_incl = e in self.included

            if is_incl:
                pct_col   = f"{mean_af.get(e, 0)*100:.2f}"
                sigma_col = f"{std_af.get(e, 0)*100:.2f}"
            else:
                pct_col   = f"({mean_af_raw.get(e, 0)*100:.2f})"
                sigma_col = f"({std_af_raw.get(e, 0)*100:.2f})"

            mvals = mass_bucket.get(e, [])
            mean_mass = float(np.mean(mvals)) if mvals else 0.0

            tags = ('incl',) if is_incl else ('excl',)
            self.tree.insert('', 'end',
                values=(e, ATOMIC_NUMBERS.get(e, ''),
                        pct_col, sigma_col,
                        f"{n_valid}/{n_total}" if is_incl else f"—/{n_total}",
                        f"{mean_mass:.2f}"),
                tags=tags)

        # ── Per-file breakdown table ──
        self._rebuild_file_table(incl_elems)

    @staticmethod
    def _alloy_str(incl_elems, f_mean):
        """Build compact alloy notation with subscript digits: Pt₄₅.₃Ni₅₄.₇"""
        return alloy_notation(incl_elems, f_mean, decimals=1)

    def _rebuild_file_table(self, incl_elems):
        """Rebuild the per-file Treeview with Composition column + one column per element."""
        if self.file_tree is not None:
            self.file_tree.destroy()
            self.file_tree = None
        for w in self._file_tbl_frame.winfo_children():
            w.destroy()
        self._file_table_data = []

        if not self.files_data or not incl_elems:
            tk.Label(self._file_tbl_frame,
                     text="(No data / no elements selected)",
                     font=('Segoe UI', 9), fg='#bdbdbd',
                     bg='#eceff1').pack(anchor='w')
            return

        # Columns: File | Composition | Total | Defined | Undefined | El1 mean | El1 σ | ...
        elem_cols = []
        for e in incl_elems:
            elem_cols += [f"{e} at%", f"{e} σ"]

        display_cols = ('File', 'Composition',
                        'Total sec.', 'Defined', 'Undefined') + tuple(elem_cols)
        self.file_tree = ttk.Treeview(
            self._file_tbl_frame, columns=display_cols, show='headings',
            height=min(len(self.files_data) + 1, 5))

        UCOL = 140  # uniform column width
        self.file_tree.heading('File', text='File')
        self.file_tree.column('File', width=UCOL, anchor='w')
        self.file_tree.heading('Composition', text='Composition')
        self.file_tree.column('Composition', width=UCOL, anchor='w')
        for col in ('Total sec.', 'Defined', 'Undefined'):
            self.file_tree.heading(col, text=col)
            self.file_tree.column(col, width=UCOL, anchor='center')
        for e in incl_elems:
            self.file_tree.heading(f"{e} at%", text=f"{e} at%")
            self.file_tree.column(f"{e} at%", width=UCOL, anchor='center')
            self.file_tree.heading(f"{e} σ", text=f"{e} σ")
            self.file_tree.column(f"{e} σ", width=UCOL, anchor='center')

        hsb = ttk.Scrollbar(self._file_tbl_frame, orient='horizontal',
                             command=self.file_tree.xview)
        vsb = ttk.Scrollbar(self._file_tbl_frame, orient='vertical',
                             command=self.file_tree.yview)
        self.file_tree.configure(xscrollcommand=hsb.set, yscrollcommand=vsb.set)
        self.file_tree.grid(row=0, column=0, sticky='nsew')
        vsb.grid(row=0, column=1, sticky='ns')
        hsb.grid(row=1, column=0, sticky='ew')
        self._file_tbl_frame.grid_rowconfigure(0, weight=1)
        self._file_tbl_frame.grid_columnconfigure(0, weight=1)

        self.file_tree.tag_configure('focused',
            background='#fff3e0', font=('Segoe UI', 22, 'bold'))
        self.file_tree.tag_configure('has_undef',
            foreground='#e65100')  # orange text if any undefined sections

        # Store header for CSV export
        self._file_table_data = [display_cols]

        for path in self.file_order:
            if path not in self.files_data:
                continue
            single = {path: self.files_data[path]}
            _, f_mean, f_std, f_valid, f_total = compute_stats_for(
                single, included=set(incl_elems))
            f_undef = f_total - f_valid

            alloy = self._alloy_str(incl_elems, f_mean)
            row = [os.path.basename(path), alloy,
                   str(f_total), str(f_valid), str(f_undef)]
            for e in incl_elems:
                if e in f_mean:
                    row.append(f"{f_mean[e]*100:.2f}")
                    row.append(f"{f_std[e]*100:.2f}")
                else:
                    row += ["—", "—"]

            self._file_table_data.append(tuple(row))
            tags = []
            if path == self.focused_path:
                tags.append('focused')
            if f_undef > 0:
                tags.append('has_undef')
            self.file_tree.insert('', 'end', values=tuple(row), tags=tuple(tags))

    def _export_file_table(self):
        """Export the per-file composition table to a CSV file."""
        if len(self._file_table_data) <= 1:
            messagebox.showinfo("Export", "No data to export.")
            return

        path = filedialog.asksaveasfilename(
            title="Save Per-file Composition",
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])
        if not path:
            return

        import csv as _csv
        try:
            with open(path, 'w', newline='', encoding='utf-8-sig') as f:
                writer = _csv.writer(f)
                for row in self._file_table_data:
                    writer.writerow(row)
            messagebox.showinfo("Export", f"Saved to:\n{path}")
        except Exception as ex:
            messagebox.showerror("Export Error", str(ex))

    def _sort_tree(self, col):
        rows = [(self.tree.item(r, 'values'), r) for r in self.tree.get_children()]
        col_idx = ('Element', 'Z', 'Mean at%', '± σ (at%)', 'N spots', 'Mean mass%').index(col)
        rev = (self._sort_col == col) and (not self._sort_rev)
        self._sort_col = col
        self._sort_rev = rev

        def key(item):
            v = str(item[0][col_idx]).strip('()')
            try:
                return (0, float(v))
            except ValueError:
                return (1, v)

        rows.sort(key=key, reverse=rev)
        for i, (_, row) in enumerate(rows):
            self.tree.move(row, '', i)


if __name__ == '__main__':
    app = App()
    app.mainloop()
