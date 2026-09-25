# -*- coding: utf-8 -*-
"""
portico2d_gui.py — Interface gráfica (Tkinter + matplotlib) do analisador de
pórticos, vigas e treliças planas. O cálculo fica todo em portico2d_core.py.

Requisitos:  pip install openseespy matplotlib
Executar:    python portico2d_gui.py
"""
import math
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, colorchooser

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import MultipleLocator

import numpy as np

import portico2d_core as core
import combinacoes_nbr8800 as nbr
import importar_dxf as dxf
import verificacao_nbr8800 as vf

PALETA = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b",
          "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", "#003f5c", "#a05195"]

MODOS = [("sel", "⬚ Selecionar"), ("mover", "✥ Mover"), ("no", "● Nó"), ("barra", "╱ Barra"),
         ("apoio", "△ Apoio"), ("cno", "↓ Carga nó"), ("cbarra", "⇊ Carga barra"),
         ("del", "✖ Excluir")]

LIGACOES = {"Rígida – Rígida (pórtico)": (False, False),
            "Rótula – Rígida": (True, False),
            "Rígida – Rótula": (False, True),
            "Rótula – Rótula (treliça)": (True, True)}

PRESETS_APOIO = {
    "Articulado fixo (2º gênero)": ("fixo", "fixo", "livre"),
    "Articulado móvel — rola em X (1º gênero)": ("livre", "fixo", "livre"),
    "Articulado móvel — rola em Y": ("fixo", "livre", "livre"),
    "Engaste (3º gênero)": ("fixo", "fixo", "fixo"),
    "Engaste deslizante em X": ("livre", "fixo", "fixo"),
    "Personalizado / molas": None,
}

DIAGRAMAS = ["Estrutura", "Deformada", "Normal (N)", "Cortante (V)", "Momento (M)", "Reações"]


def fnum(s):
    return float(str(s).strip().replace(",", "."))


def fmt(v, nd=2):
    if abs(v) < 0.5 * 10 ** -nd:
        v = 0.0
    return f"{v:.{nd}f}".replace(".", ",")


# =============================================================================
# DIÁLOGOS
# =============================================================================
def desenha_secao(cv, secao, beta, cor="#1f77b4"):
    """Desenha a seção (girada de β) num tk.Canvas, com eixos locais y (no plano) e z."""
    cv.delete("all")
    W, H = int(cv["width"]), int(cv["height"])
    try:
        pols = secao.poligonos(beta)
    except Exception:
        return
    ext = max(max(abs(v) for p in pols for pt in p for v in pt), 1e-6)
    esc = 0.36 * min(W, H) / ext
    cx, cy = W / 2, H / 2
    # eixos
    L = 0.46 * min(W, H)
    cv.create_line(cx, cy, cx, cy - L, fill="#2e7d32", arrow="last", width=1.5)
    cv.create_text(cx + 9, cy - L + 4, text="y", fill="#2e7d32", font=("Segoe UI", 9, "bold"))
    cv.create_line(cx, cy, cx + L, cy, fill="#1565c0", arrow="last", width=1.5, dash=(4, 2))
    cv.create_text(cx + L - 2, cy + 10, text="z", fill="#1565c0", font=("Segoe UI", 9, "bold"))
    for p in pols:
        pts = []
        for z, y in p:
            pts += [cx + z * esc, cy - y * esc]
        cv.create_polygon(pts, fill=cor, outline="#222", width=1.2, stipple="gray50")
        cv.create_polygon(pts, fill="", outline="#222", width=1.2)
    cv.create_oval(cx - 3, cy - 3, cx + 3, cy + 3, fill="#d62728", outline="")
    cv.create_text(6, H - 6, anchor="sw", text=f"β = {beta:g}°   (y = no plano do pórtico)",
                   fill="#555", font=("Segoe UI", 8))


class BarraDialog(tk.Toplevel):
    """Propriedades de uma ou várias barras: seção, ligação, giro β e peso próprio."""

    def __init__(self, parent, modelo, barras):
        super().__init__(parent)
        self.transient(parent)
        self.resizable(False, False)
        self.modelo, self.result = modelo, None
        b0 = barras[0]
        self.title(f"Barra {b0.id}" if len(barras) == 1 else f"{len(barras)} barras")
        frm = ttk.Frame(self, padding=12)
        frm.pack()
        if len(barras) == 1:
            L = modelo.geometria(b0)[0]
            ttk.Label(frm, text=f"Nós {b0.ni} → {b0.nj}    L = {fmt(L, 3)} m", foreground="#444").grid(
                row=0, column=0, columnspan=2, sticky="w", pady=(0, 6))
        self.v_sec = tk.StringVar(value=b0.secao)
        ttk.Label(frm, text="Seção").grid(row=1, column=0, sticky="w")
        c = ttk.Combobox(frm, textvariable=self.v_sec, values=list(modelo.secoes), state="readonly", width=26)
        c.grid(row=1, column=1, sticky="we", pady=2)
        lig = [k for k, v in LIGACOES.items() if v == (b0.rot_i, b0.rot_j)][0]
        self.v_lig = tk.StringVar(value=lig)
        ttk.Label(frm, text="Ligação (i – j)").grid(row=2, column=0, sticky="w")
        ttk.Combobox(frm, textvariable=self.v_lig, values=list(LIGACOES), state="readonly", width=26).grid(
            row=2, column=1, sticky="we", pady=2)
        ttk.Label(frm, text="Giro do perfil β (°)").grid(row=3, column=0, sticky="w")
        fb = ttk.Frame(frm)
        fb.grid(row=3, column=1, sticky="w")
        self.v_beta = tk.StringVar(value=f"{b0.beta:g}")
        ttk.Combobox(fb, textvariable=self.v_beta, values=["0", "90", "180", "270"], width=6).pack(side="left")
        ttk.Button(fb, text="↻ +90°", width=7, command=self._gira).pack(side="left", padx=4)
        self.v_pp = tk.BooleanVar(value=b0.pp)
        ttk.Checkbutton(frm, text="Considerar peso próprio desta barra", variable=self.v_pp).grid(
            row=4, column=0, columnspan=2, sticky="w", pady=(4, 0))
        self.lbl = ttk.Label(frm, text="", foreground="#0b5394")
        self.lbl.grid(row=5, column=0, columnspan=2, sticky="w", pady=(6, 0))
        self.cv = tk.Canvas(frm, width=210, height=210, bg="white", highlightthickness=1, highlightbackground="#ccc")
        self.cv.grid(row=0, column=2, rowspan=7, padx=(12, 0))
        bf = ttk.Frame(frm)
        bf.grid(row=6, column=0, columnspan=2, sticky="w", pady=(10, 0))
        ttk.Button(bf, text="OK", command=self._ok).pack(side="left", padx=4)
        ttk.Button(bf, text="Cancelar", command=self.destroy).pack(side="left")
        for v in (self.v_sec, self.v_beta):
            v.trace_add("write", lambda *a: self._preview())
        self._preview()
        self.bind("<Return>", lambda e: self._ok())
        self.grab_set()
        self.wait_window()

    def _beta(self):
        return fnum(self.v_beta.get() or 0) % 360

    def _gira(self):
        try:
            self.v_beta.set(f"{(self._beta() + 90) % 360:g}")
        except ValueError:
            self.v_beta.set("0")

    def _preview(self):
        s = self.modelo.secoes.get(self.v_sec.get())
        try:
            beta = self._beta()
        except ValueError:
            return
        if s:
            desenha_secao(self.cv, s, beta, s.cor)
            self.lbl.configure(text=f"I no plano = {s.I_plano(beta):.1f} cm⁴")

    def _ok(self):
        try:
            beta = self._beta()
        except ValueError:
            messagebox.showerror("Barra", "β inválido.", parent=self)
            return
        self.result = dict(secao=self.v_sec.get(), lig=LIGACOES[self.v_lig.get()], beta=beta, pp=self.v_pp.get())
        self.destroy()


class FormDialog(tk.Toplevel):
    """Formulário genérico. campos: [{'key','label','type':float|int|str|combo|check,'default','options'}]"""

    def __init__(self, parent, titulo, campos, texto=None):
        super().__init__(parent)
        self.title(titulo)
        self.transient(parent)
        self.resizable(False, False)
        self.result = None
        self.vars = {}
        frm = ttk.Frame(self, padding=12)
        frm.pack(fill="both", expand=True)
        r = 0
        if texto:
            ttk.Label(frm, text=texto, wraplength=340, foreground="#444", justify="left").grid(
                row=r, column=0, columnspan=2, sticky="w", pady=(0, 8))
            r += 1
        first = None
        for c in campos:
            k, t, d = c["key"], c.get("type", "float"), c.get("default", "")
            if t == "check":
                v = tk.BooleanVar(value=bool(d))
                ttk.Checkbutton(frm, text=c["label"], variable=v).grid(row=r, column=0, columnspan=2, sticky="w", pady=2)
            else:
                ttk.Label(frm, text=c["label"]).grid(row=r, column=0, sticky="w", padx=(0, 10), pady=2)
                v = tk.StringVar(value=str(d).replace(".", ",") if t == "float" else str(d))
                if t == "combo":
                    w = ttk.Combobox(frm, textvariable=v, values=c["options"], state="readonly", width=30)
                else:
                    w = ttk.Entry(frm, textvariable=v, width=32)
                    first = first or w
                w.grid(row=r, column=1, sticky="we", pady=2)
            self.vars[k] = (v, t)
            r += 1
        bf = ttk.Frame(frm)
        bf.grid(row=r, column=0, columnspan=2, pady=(10, 0), sticky="e")
        ttk.Button(bf, text="OK", command=self._ok).pack(side="left", padx=4)
        ttk.Button(bf, text="Cancelar", command=self.destroy).pack(side="left")
        self.bind("<Return>", lambda e: self._ok())
        self.bind("<Escape>", lambda e: self.destroy())
        if first:
            first.focus_set()
            first.select_range(0, "end")
        self.grab_set()
        self.wait_window()

    def _ok(self):
        out = {}
        for k, (v, t) in self.vars.items():
            try:
                if t == "float":
                    out[k] = fnum(v.get())
                elif t == "int":
                    out[k] = int(v.get())
                else:
                    out[k] = v.get()
            except ValueError:
                messagebox.showerror("Valor inválido", f"Valor inválido no campo '{k}'.", parent=self)
                return
        self.result = out
        self.destroy()


class SecaoDialog(tk.Toplevel):
    CAMPOS = {"W": [], "I": [("d", "d (mm)"), ("bf", "bf (mm)"), ("tw", "tw (mm)"), ("tf", "tf (mm)")],
              "L": [("b", "b — aba (mm)"), ("t", "t (mm)")],
              "U": [("bw", "bw — alma (mm)"), ("bf", "bf — mesa (mm)"), ("t", "t (mm)")],
              "Ue": [("bw", "bw — alma (mm)"), ("bf", "bf — mesa (mm)"), ("D", "D — enrijecedor (mm)"), ("t", "t (mm)")],
              "GEN": [("A", "A (cm²)"), ("I", "Ix (cm⁴)"), ("Iy", "Iy (cm⁴)")]}
    CATS = {"L": core.CATALOGO_L, "U": core.CATALOGO_U, "Ue": core.CATALOGO_UE}

    def __init__(self, parent, modelo, secao=None, cor=None, nome=None):
        super().__init__(parent)
        self.title("Seção transversal" + (f" — layer '{nome}'" if nome else ""))
        self.transient(parent)
        self.resizable(False, False)
        self.modelo, self.original, self.result = modelo, secao, None
        s = secao or core.Secao("", "W", "ASTM A572 Gr50", {"catalogo": "W200x15,0"}, cor=cor or PALETA[0])
        self.cor = s.cor
        self.nome_manual = secao is not None
        if secao is None and nome:
            s.nome, self.nome_manual = nome, True
        frm = ttk.Frame(self, padding=12)
        frm.pack(fill="both", expand=True)
        self.frm = frm
        fam_labels = [f"{k} — {v}" for k, v in core.FAMILIAS.items()]
        self.v_fam = tk.StringVar(value=f"{s.familia} — {core.FAMILIAS[s.familia]}")
        ttk.Label(frm, text="Família").grid(row=0, column=0, sticky="w")
        cb = ttk.Combobox(frm, textvariable=self.v_fam, values=fam_labels, state="readonly", width=44)
        cb.grid(row=0, column=1, sticky="we", pady=2)
        cb.bind("<<ComboboxSelected>>", lambda e: self._troca_familia(True))
        ttk.Label(frm, text="Material").grid(row=1, column=0, sticky="w")
        self.v_mat = tk.StringVar(value=s.material)
        ttk.Combobox(frm, textvariable=self.v_mat, values=list(modelo.materiais), state="readonly").grid(
            row=1, column=1, sticky="we", pady=2)
        self.dyn = ttk.Frame(frm)
        self.dyn.grid(row=2, column=0, columnspan=2, sticky="we", pady=6)
        ttk.Label(frm, text="Quantidade").grid(row=3, column=0, sticky="w")
        self.v_qtd = tk.StringVar(value="2 — composto (costas c/ costas)" if s.qtd == 2 else "1 — simples")
        self.cb_qtd = ttk.Combobox(frm, textvariable=self.v_qtd, state="readonly",
                                   values=["1 — simples", "2 — composto (costas c/ costas)"])
        self.cb_qtd.grid(row=3, column=1, sticky="we", pady=2)
        ttk.Label(frm, text="Afastamento entre perfis (mm)").grid(row=4, column=0, sticky="w")
        self.v_gap = tk.StringVar(value=str(s.dims.get("gap", 0.0)).replace(".", ","))
        self.en_gap = ttk.Entry(frm, textvariable=self.v_gap, width=10)
        self.en_gap.grid(row=4, column=1, sticky="w", pady=2)
        self.cv = tk.Canvas(frm, width=230, height=230, bg="white", highlightthickness=1, highlightbackground="#ccc")
        self.cv.grid(row=0, column=2, rowspan=8, padx=(12, 0), sticky="n")
        ttk.Label(frm, text="Nome / camada").grid(row=5, column=0, sticky="w")
        self.v_nome = tk.StringVar(value=s.nome)
        en = ttk.Entry(frm, textvariable=self.v_nome)
        en.grid(row=5, column=1, sticky="we", pady=2)
        en.bind("<Key>", lambda e: setattr(self, "nome_manual", True))
        ttk.Label(frm, text="Cor da camada").grid(row=6, column=0, sticky="w")
        self.btn_cor = tk.Button(frm, text="      ", bg=self.cor, command=self._escolhe_cor, relief="groove")
        self.btn_cor.grid(row=6, column=1, sticky="w", pady=2)
        self.lbl_props = ttk.Label(frm, text="", foreground="#0b5394", wraplength=420, justify="left")
        self.lbl_props.grid(row=7, column=0, columnspan=2, sticky="w", pady=(8, 0))
        self.fr_over = ttk.LabelFrame(frm, text="Propriedades para dimensionamento — opcional (vazio = calcular)",
                                      padding=6)
        self.fr_over.grid(row=8, column=0, columnspan=3, sticky="we", pady=(8, 0))
        self.v_over = {}
        for j, (k, lab) in enumerate((("r", "r raio (mm)"), ("Zx", "Zx (cm³)"), ("Zy", "Zy (cm³)"),
                                      ("J", "J = It (cm⁴)"), ("Cw", "Cw (cm⁶)"))):
            ttk.Label(self.fr_over, text=lab).grid(row=0, column=2 * j, sticky="w", padx=(0, 3))
            v = tk.StringVar(value=str(s.dims.get(k, "")).replace(".", ","))
            ttk.Entry(self.fr_over, textvariable=v, width=8).grid(row=0, column=2 * j + 1, padx=(0, 8))
            self.v_over[k] = v
        ttk.Label(self.fr_over, text="Use os valores do catálogo do fabricante para reproduzir exatamente o memorial "
                                     "(sem eles: Zx/Zy com raios derivados da área, J sem raios — conservador).",
                  foreground="#666", wraplength=640).grid(row=1, column=0, columnspan=10, sticky="w", pady=(4, 0))
        bf = ttk.Frame(frm)
        bf.grid(row=9, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(bf, text="OK", command=self._ok).pack(side="left", padx=4)
        ttk.Button(bf, text="Cancelar", command=self.destroy).pack(side="left")
        for v in (self.v_qtd, self.v_mat, self.v_gap):
            v.trace_add("write", lambda *a: self._atualiza())
        self.vdims = {}
        self._troca_familia(False, s.dims)
        self.grab_set()
        self.wait_window()

    def fam(self):
        return self.v_fam.get().split(" — ")[0]

    def _escolhe_cor(self):
        c = colorchooser.askcolor(self.cor, parent=self)[1]
        if c:
            self.cor = c
            self.btn_cor.configure(bg=c)
            self._atualiza()

    def _troca_familia(self, mudou, dims=None):
        for w in self.dyn.winfo_children():
            w.destroy()
        f = self.fam()
        if mudou:
            self.v_mat.set(core.material_padrao_da_familia(f))
        self.vdims = {}
        self.v_preset = tk.StringVar()
        r = 0
        if f == "W":
            ttk.Label(self.dyn, text="Perfil (catálogo)").grid(row=0, column=0, sticky="w")
            self.v_preset.set((dims or {}).get("catalogo", "W200x15,0"))
            c = ttk.Combobox(self.dyn, textvariable=self.v_preset, values=list(core.CATALOGO_W), state="readonly", width=20)
            c.grid(row=0, column=1, sticky="w")
            c.bind("<<ComboboxSelected>>", lambda e: self._atualiza())
            r = 1
        elif f in self.CATS:
            ttk.Label(self.dyn, text="Padrão (preenche)").grid(row=0, column=0, sticky="w")
            c = ttk.Combobox(self.dyn, textvariable=self.v_preset, values=list(self.CATS[f]), state="readonly", width=24)
            c.grid(row=0, column=1, sticky="w")
            c.bind("<<ComboboxSelected>>", lambda e: self._aplica_preset())
            r = 1
        defaults = {"I": dict(d=300, bf=150, tw=6.3, tf=9.5), "L": dict(b=50.8, t=6.35),
                    "U": dict(bw=100, bf=50, t=2.65), "Ue": dict(bw=150, bf=60, D=20, t=2.65),
                    "GEN": dict(A=10.0, I=100.0, Iy=100.0)}.get(f, {})
        if dims and not mudou:
            defaults.update({k: v for k, v in dims.items() if k != "catalogo"})
        for k, lab in self.CAMPOS[f]:
            ttk.Label(self.dyn, text=lab).grid(row=r, column=0, sticky="w", pady=1)
            v = tk.StringVar(value=str(defaults.get(k, "")).replace(".", ","))
            v.trace_add("write", lambda *a: self._atualiza())
            ttk.Entry(self.dyn, textvariable=v, width=12).grid(row=r, column=1, sticky="w", pady=1)
            self.vdims[k] = v
            r += 1
        self._atualiza()

    def _aplica_preset(self):
        f = self.fam()
        vals = self.CATS[f][self.v_preset.get()]
        keys = [k for k, _ in self.CAMPOS[f]]
        for k, v in zip(keys, vals):
            self.vdims[k].set(str(v).replace(".", ","))
        if not self.nome_manual:
            self.v_nome.set(self.v_preset.get())
        self._atualiza()

    def _monta(self):
        f = self.fam()
        if f == "W":
            dims = {"catalogo": self.v_preset.get()}
        else:
            dims = {k: fnum(v.get()) for k, v in self.vdims.items()}
        if hasattr(self, "fr_over"):
            if f in ("W", "I"):
                self.fr_over.grid()
                for k, v in self.v_over.items():
                    if v.get().strip():
                        dims[k] = fnum(v.get())
            else:
                self.fr_over.grid_remove()
        comp = f in core.FAMILIAS_COMPOSTAS
        self.cb_qtd.configure(state="readonly" if comp else "disabled")
        if not comp:
            self.v_qtd.set("1 — simples")
        qtd = 2 if self.v_qtd.get().startswith("2") else 1
        self.en_gap.configure(state="normal" if qtd == 2 else "disabled")
        if qtd == 2:
            dims["gap"] = fnum(self.v_gap.get() or 0)
        nome = self.v_nome.get().strip()
        if not self.nome_manual or not nome:
            if f == "W":
                base = dims["catalogo"]
            elif f == "I":
                base = "I {d:g}x{bf:g}x{tw:g}x{tf:g}".format(**dims)
            elif f == "L":
                base = "L {b:g}x{t:g}".format(**dims)
            elif f == "U":
                base = "U {bw:g}x{bf:g}x{t:g}".format(**dims)
            elif f == "Ue":
                base = "Ue {bw:g}x{bf:g}x{D:g}x{t:g}".format(**dims)
            else:
                base = "GEN A{A:g} I{I:g}".format(**dims)
            nome = ("2x " if qtd == 2 else "") + base.replace(".", ",")
            self.v_nome.set(nome)
        return core.Secao(nome, f, self.v_mat.get(), dims, "x", qtd, self.cor)

    def _atualiza(self):
        try:
            s = self._monta()
            p = s.propriedades()
            txt = (f"Perfil simples: A = {p['A1']:.2f} cm²   Ix = {p['Ix1']:.1f} cm⁴   Iy = {p['Iy1']:.1f} cm⁴\n"
                   f"Seção: A = {p['A']:.2f} cm²   Ix = {p['Ix']:.1f} cm⁴   Iy = {p['Iy']:.1f} cm⁴\n"
                   f"Peso ≈ {p['A'] * 1e-4 * 7850:.2f} kg/m   —   {p['obs']}\n"
                   "Flexão no plano: Ix com β = 0°, Iy com β = 90° (β é definido na barra).")
            self.lbl_props.configure(text=txt)
            desenha_secao(self.cv, s, 0.0, self.cor)
        except (ValueError, KeyError, ZeroDivisionError, TypeError):
            self.lbl_props.configure(text="(preencha as dimensões)")

    def _ok(self):
        try:
            s = self._monta()
            p = s.propriedades()
            assert p["A"] > 0 and p["Ix"] > 0 and p["Iy"] > 0
        except Exception:
            messagebox.showerror("Seção", "Dimensões inválidas.", parent=self)
            return
        if s.nome in self.modelo.secoes and (self.original is None or self.original.nome != s.nome):
            messagebox.showerror("Seção", "Já existe uma seção com esse nome.", parent=self)
            return
        self.result = s
        self.destroy()


class ApoioDialog(tk.Toplevel):
    def __init__(self, parent, nid, apoio=None):
        super().__init__(parent)
        self.title(f"Apoio no nó {nid}")
        self.transient(parent)
        self.resizable(False, False)
        self.result, self.remover = None, False
        a = apoio or core.Apoio(nid)
        frm = ttk.Frame(self, padding=12)
        frm.pack()
        ttk.Label(frm, text="Tipo").grid(row=0, column=0, sticky="w")
        self.v_pre = tk.StringVar(value="Personalizado / molas")
        for k, v in PRESETS_APOIO.items():
            if v == (a.ux, a.uy, a.rz):
                self.v_pre.set(k)
        cb = ttk.Combobox(frm, textvariable=self.v_pre, values=list(PRESETS_APOIO), state="readonly", width=40)
        cb.grid(row=0, column=1, columnspan=2, sticky="we", pady=(0, 8))
        cb.bind("<<ComboboxSelected>>", self._preset)
        self.vt, self.vk = [], []
        for r, (lab, t, k, un) in enumerate((("ux (horizontal)", a.ux, a.kx, "kN/m"),
                                            ("uy (vertical)", a.uy, a.ky, "kN/m"),
                                            ("rz (rotação)", a.rz, a.kr, "kN·m/rad")), start=1):
            ttk.Label(frm, text=lab).grid(row=r, column=0, sticky="w")
            vt = tk.StringVar(value=t)
            ttk.Combobox(frm, textvariable=vt, values=core.TIPOS_GL, state="readonly", width=8).grid(row=r, column=1, pady=2)
            vk = tk.StringVar(value=str(k).replace(".", ","))
            f2 = ttk.Frame(frm)
            f2.grid(row=r, column=2, sticky="w")
            ttk.Entry(f2, textvariable=vk, width=10).pack(side="left")
            ttk.Label(f2, text=f" k ({un}) se mola").pack(side="left")
            self.vt.append(vt)
            self.vk.append(vk)
        bf = ttk.Frame(frm)
        bf.grid(row=5, column=0, columnspan=3, sticky="e", pady=(10, 0))
        if apoio:
            ttk.Button(bf, text="Remover apoio", command=self._rem).pack(side="left", padx=4)
        ttk.Button(bf, text="OK", command=lambda: self._ok(nid)).pack(side="left", padx=4)
        ttk.Button(bf, text="Cancelar", command=self.destroy).pack(side="left")
        self.grab_set()
        self.wait_window()

    def _preset(self, e=None):
        v = PRESETS_APOIO[self.v_pre.get()]
        if v:
            for var, t in zip(self.vt, v):
                var.set(t)

    def _rem(self):
        self.remover = True
        self.destroy()

    def _ok(self, nid):
        try:
            k = [fnum(v.get() or 0) for v in self.vk]
        except ValueError:
            messagebox.showerror("Apoio", "Rigidez inválida.", parent=self)
            return
        t = [v.get() for v in self.vt]
        for ti, ki in zip(t, k):
            if ti == "mola" and ki <= 0:
                messagebox.showerror("Apoio", "Mola exige rigidez > 0.", parent=self)
                return
        if all(x == "livre" for x in t):
            messagebox.showerror("Apoio", "Todos os graus livres: use 'Remover apoio'.", parent=self)
            return
        self.result = core.Apoio(nid, t[0], t[1], t[2], k[0], k[1], k[2])
        self.destroy()


class CombDialog(tk.Toplevel):
    def __init__(self, parent, modelo, comb=None):
        super().__init__(parent)
        self.title("Combinação de carregamentos")
        self.transient(parent)
        self.result = None
        frm = ttk.Frame(self, padding=12)
        frm.pack()
        ttk.Label(frm, text="Nome").grid(row=0, column=0, sticky="w")
        self.v_nome = tk.StringVar(value=comb.nome if comb else f"C{len(modelo.combinacoes) + 1}")
        ttk.Entry(frm, textvariable=self.v_nome, width=24).grid(row=0, column=1, sticky="w")
        ttk.Label(frm, text="Fatores γ por caso (0 = não entra):", foreground="#444").grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(8, 2))
        self.vf = {}
        for r, c in enumerate(modelo.casos, start=2):
            ttk.Label(frm, text=c).grid(row=r, column=0, sticky="w")
            v = tk.StringVar(value=str((comb.fatores.get(c, 0.0) if comb else 1.0)).replace(".", ","))
            ttk.Entry(frm, textvariable=v, width=8).grid(row=r, column=1, sticky="w", pady=1)
            self.vf[c] = v
        bf = ttk.Frame(frm)
        bf.grid(row=99, column=0, columnspan=2, sticky="e", pady=(10, 0))
        ttk.Button(bf, text="OK", command=self._ok).pack(side="left", padx=4)
        ttk.Button(bf, text="Cancelar", command=self.destroy).pack(side="left")
        self.grab_set()
        self.wait_window()

    def _ok(self):
        try:
            f = {c: fnum(v.get() or 0) for c, v in self.vf.items()}
        except ValueError:
            messagebox.showerror("Combinação", "Fator inválido.", parent=self)
            return
        nome = self.v_nome.get().strip()
        if not nome:
            return
        self.result = core.Combinacao(nome, f)
        self.destroy()

class ProjetoBarraDialog(tk.Toplevel):
    """Parâmetros de dimensionamento de uma ou várias barras. Cada comprimento pode ser dado
    como K × L (L = comprimento da barra) OU diretamente em metros (prevalece o valor em metros)."""
    LINHAS = [("kp", "Lp", "Flambagem por flexão NO PLANO"),
              ("kf", "Lf", "Flambagem por flexão FORA DO PLANO"),
              ("kz", "Lz", "Flambagem por torção (KzLz)"),
              ("klb", "Lb", "Contenção lateral p/ FLT (Lb)")]

    def __init__(self, parent, modelo, barras):
        super().__init__(parent)
        self.title("Parâmetros de dimensionamento" + (f" — barra {barras[0].id}" if len(barras) == 1
                                                     else f" — {len(barras)} barras"))
        self.transient(parent)
        self.resizable(False, False)
        self.result = None
        p = vf.pproj(barras[0])
        L = modelo.geometria(barras[0])[0]
        self.L = L
        frm = ttk.Frame(self, padding=12)
        frm.pack()
        ttk.Label(frm, text=(f"Barra {barras[0].id}: L = {fmt(L, 3)} m" if len(barras) == 1 else
                             f"Aplicar às {len(barras)} barras selecionadas (valores da barra {barras[0].id}, "
                             f"L = {fmt(L, 3)} m)"), foreground="#444").grid(row=0, column=0, columnspan=4,
                                                                            sticky="w", pady=(0, 6))
        for j, t in enumerate(("", "K (× L)", "ou comprimento (m)", "Efetivo")):
            ttk.Label(frm, text=t, font=("Segoe UI", 9, "bold")).grid(row=1, column=j, sticky="w", padx=4)
        self.v, self.lbl_ef = {}, {}
        for i, (kk, lk, lab) in enumerate(self.LINHAS, start=2):
            ttk.Label(frm, text=lab).grid(row=i, column=0, sticky="w", pady=1)
            vk = tk.StringVar(value=fmt(p[kk], 3).rstrip("0").rstrip(","))
            vl = tk.StringVar(value="" if p.get(lk) in (None, "") else fmt(float(p[lk]), 3))
            ttk.Entry(frm, textvariable=vk, width=8).grid(row=i, column=1, sticky="w", padx=4)
            ttk.Entry(frm, textvariable=vl, width=10).grid(row=i, column=2, sticky="w", padx=4)
            lb = ttk.Label(frm, text="", foreground="#0b5394", width=16)
            lb.grid(row=i, column=3, sticky="w", padx=4)
            self.v[kk], self.v[lk], self.lbl_ef[kk] = vk, vl, lb
            for var in (vk, vl):
                var.trace_add("write", lambda *a: self._efetivos())
        ttk.Label(frm, text="Preencha o comprimento em metros quando a barra foi dividida por nós que NÃO travam a "
                            "flambagem (ex.: nó só para lançar carga). Vazio = usa K × L. Lb = 0 → mesa travada "
                            "continuamente.", foreground="#666", wraplength=470, justify="left").grid(
            row=6, column=0, columnspan=4, sticky="w", pady=(4, 8))
        r = 7
        ttk.Label(frm, text="Cb ('auto' = 5.4.2.3)").grid(row=r, column=0, sticky="w")
        self.v["cb"] = tk.StringVar(value=str(p["cb"]).replace(".", ","))
        ttk.Entry(frm, textvariable=self.v["cb"], width=8).grid(row=r, column=1, sticky="w", padx=4)
        r += 1
        ttk.Separator(frm).grid(row=r, column=0, columnspan=4, sticky="we", pady=6)
        r += 1
        ttk.Label(frm, text="Flecha limite:  Lref /").grid(row=r, column=0, sticky="w")
        self.v["flecha"] = tk.StringVar(value=f"{p['flecha']:g}")
        ttk.Combobox(frm, textvariable=self.v["flecha"], width=7,
                     values=["120", "180", "250", "300", "350", "500"]).grid(row=r, column=1, sticky="w", padx=4)
        r += 1
        ttk.Label(frm, text="Vão de referência Lref (m)").grid(row=r, column=0, sticky="w")
        self.v["Lflecha"] = tk.StringVar(value="" if not p.get("Lflecha") else fmt(float(p["Lflecha"]), 3))
        ttk.Entry(frm, textvariable=self.v["Lflecha"], width=10).grid(row=r, column=1, sticky="w", padx=4)
        ttk.Label(frm, text="vazio = L da barra", foreground="#666").grid(row=r, column=2, columnspan=2, sticky="w")
        r += 1
        ttk.Label(frm, text="Deslocamento medido").grid(row=r, column=0, sticky="w")
        self.MODOS = {"Relativo à corda da barra": "corda",
                      "Absoluto (em relação à posição original)": "absoluta"}
        self.v_modo = tk.StringVar(value={v: k for k, v in self.MODOS.items()}[p.get("flecha_modo", "corda")])
        ttk.Combobox(frm, textvariable=self.v_modo, values=list(self.MODOS), state="readonly", width=38).grid(
            row=r, column=1, columnspan=3, sticky="w", padx=4)
        r += 1
        self.v_bal = tk.BooleanVar(value=p["balanco"])
        ttk.Checkbutton(frm, text="Barra em balanço (Lref = 2 × comprimento; Cb = 1,0)",
                        variable=self.v_bal).grid(row=r, column=0, columnspan=4, sticky="w", pady=(4, 0))
        r += 1
        ttk.Label(frm, text="Viga dividida por nós de carga: use 'Absoluto' e informe o vão total em Lref — o limite "
                            "passa a ser vão/x e o deslocamento é o total, não o relativo a cada trecho.  "
                            "Tabela C.1: L/350 vigas que suportam paredes; L/250 vigas de piso; L/180 terças e vigas "
                            "de cobertura.", foreground="#666", wraplength=470, justify="left").grid(
            row=r, column=0, columnspan=4, sticky="w", pady=(4, 0))
        r += 1
        ttk.Separator(frm).grid(row=r, column=0, columnspan=4, sticky="we", pady=6)
        r += 1
        for k, lab in (("ct", "Ct — coef. de redução da área líquida"), ("an", "An/Ag — área líquida / bruta")):
            ttk.Label(frm, text=lab).grid(row=r, column=0, sticky="w")
            self.v[k] = tk.StringVar(value=fmt(p[k], 3))
            ttk.Entry(frm, textvariable=self.v[k], width=8).grid(row=r, column=1, sticky="w", padx=4)
            r += 1
        bf = ttk.Frame(frm)
        bf.grid(row=r, column=0, columnspan=4, sticky="e", pady=(10, 0))
        ttk.Button(bf, text="OK", command=self._ok).pack(side="left", padx=4)
        ttk.Button(bf, text="Cancelar", command=self.destroy).pack(side="left")
        self._efetivos()
        self.grab_set()
        self.wait_window()

    def _efetivos(self):
        for kk, lk, _ in self.LINHAS:
            try:
                t = self.v[lk].get().strip()
                val = fnum(t) if t else fnum(self.v[kk].get()) * self.L
                self.lbl_ef[kk].configure(text=f"= {fmt(val, 3)} m" + (" (manual)" if t else ""))
            except ValueError:
                self.lbl_ef[kk].configure(text="?")

    def _ok(self):
        out = {"balanco": self.v_bal.get(), "flecha_modo": self.MODOS[self.v_modo.get()]}
        try:
            for kk, lk, _ in self.LINHAS:
                out[kk] = fnum(self.v[kk].get())
                t = self.v[lk].get().strip()
                out[lk] = fnum(t) if t else None
                assert out[kk] >= 0 and (out[lk] is None or out[lk] >= 0)
            assert out["kp"] > 0 or out["Lp"], "K no plano"
            assert out["kf"] > 0 or out["Lf"], "K fora do plano"
            t = self.v["cb"].get().strip()
            out["cb"] = "auto" if t.lower() in ("auto", "") else fnum(t)
            out["flecha"] = fnum(self.v["flecha"].get())
            assert out["flecha"] > 0
            t = self.v["Lflecha"].get().strip()
            out["Lflecha"] = fnum(t) if t else None
            out["ct"], out["an"] = fnum(self.v["ct"].get()), fnum(self.v["an"].get())
        except (ValueError, AssertionError):
            messagebox.showerror("Parâmetros", "Valor inválido. Comprimentos no plano e fora do plano devem ser > 0.",
                                 parent=self)
            return
        self.result = out
        self.destroy()


def abrir_html(conteudo, nome="memorial.html"):
    import os
    import tempfile
    import webbrowser
    p = os.path.join(tempfile.gettempdir(), nome)
    with open(p, "w", encoding="utf-8") as f:
        f.write(conteudo)
    webbrowser.open("file://" + p.replace(os.sep, "/"))
    return p


def janela_texto(parent, titulo, texto, arquivo_padrao="memorial.txt", html_doc=None):
    w = tk.Toplevel(parent)
    w.title(titulo)
    w.geometry("1060x760")
    bar = ttk.Frame(w, padding=4)
    bar.pack(fill="x")
    t = tk.Text(w, font=("Consolas", 9), wrap="none")

    def salvar():
        p = filedialog.asksaveasfilename(parent=w, defaultextension=".txt", initialfile=arquivo_padrao,
                                         filetypes=[("Texto", "*.txt")])
        if p:
            with open(p, "w", encoding="utf-8") as f:
                f.write(t.get("1.0", "end"))

    def copiar():
        w.clipboard_clear()
        w.clipboard_append(t.get("1.0", "end"))
    if html_doc:
        def salvar_html():
            p = filedialog.asksaveasfilename(parent=w, defaultextension=".html",
                                             initialfile=arquivo_padrao.replace(".txt", ".html"),
                                             filetypes=[("Página HTML (imprimível em PDF)", "*.html")])
            if p:
                with open(p, "w", encoding="utf-8") as f:
                    f.write(html_doc)
        tk.Button(bar, text="🖨 Abrir memorial formatado (HTML / PDF)", bg="#e8f5e9",
                  command=lambda: abrir_html(html_doc, arquivo_padrao.replace(".txt", ".html"))).pack(side="left", padx=2)
        ttk.Button(bar, text="Salvar HTML", command=salvar_html).pack(side="left", padx=2)
    ttk.Button(bar, text="Salvar .txt", command=salvar).pack(side="left", padx=2)
    ttk.Button(bar, text="Copiar", command=copiar).pack(side="left", padx=2)
    sy = ttk.Scrollbar(w, command=t.yview)
    sx = ttk.Scrollbar(w, orient="horizontal", command=t.xview)
    t.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
    sy.pack(side="right", fill="y")
    sx.pack(side="bottom", fill="x")
    t.pack(fill="both", expand=True)
    t.insert("1.0", texto)
    t.tag_configure("ok", foreground="#1b8a2f", font=("Consolas", 9, "bold"))
    t.tag_configure("falha", foreground="#8b0000", font=("Consolas", 9, "bold"))
    for padrao, tag in (("PASSA", "ok"), ("NÃO PASSA", "falha"), ("NÃO ATENDE", "falha")):
        ini = "1.0"
        while True:
            pos = t.search(padrao, ini, stopindex="end")
            if not pos:
                break
            fim = f"{pos}+{len(padrao)}c"
            t.tag_add(tag, pos, fim)
            ini = fim
    return w


class GerarCombDialog(tk.Toplevel):
    """Classificação dos casos (γ, ψ, grupos) e geração automática — NBR 8800 (Tabelas 1 e 2)."""
    NAT = {"Permanente": "permanente", "Variável": "variavel", "Não combinar": "ignorar"}

    def __init__(self, parent, modelo):
        super().__init__(parent)
        self.title("Combinações de ações — NBR 8800 (seção 4.7)")
        self.transient(parent)
        self.modelo, self.aplicou = modelo, None
        nbr.garantir_classificacao(modelo)
        frm = ttk.Frame(self, padding=10)
        frm.pack(fill="both", expand=True)
        top = ttk.Frame(frm)
        top.pack(fill="x")
        ttk.Label(top, text="Combinações últimas:").pack(side="left")
        self.v_tipo = tk.StringVar(value="Normais")
        cbt = ttk.Combobox(top, textvariable=self.v_tipo, values=["Normais", "Especiais ou de construção"],
                           state="readonly", width=26)
        cbt.pack(side="left", padx=6)
        ttk.Button(top, text="Recarregar γ e ψ das tabelas", command=self._recarregar).pack(side="left", padx=6)
        ttk.Label(frm, text="Casos com o mesmo GRUPO são mutuamente exclusivos (ex.: todos os ventos no grupo "
                            "'Vento'). Valores editáveis — confira com as Tabelas 1 e 2 da NBR 8800:2024.",
                  foreground="#555", wraplength=1000, justify="left").pack(anchor="w", pady=(6, 4))
        tab = ttk.Frame(frm)
        tab.pack(fill="x")
        cab = ["Caso", "Natureza", "Categoria (Tabela 1 / Tabela 2)", "Grupo", "γ desf.", "γ fav.", "γq", "ψ0", "ψ1", "ψ2"]
        for j, t in enumerate(cab):
            ttk.Label(tab, text=t, font=("Segoe UI", 9, "bold")).grid(row=0, column=j, padx=2, sticky="w")
        self.linhas = {}
        for i, c in enumerate(modelo.casos.values(), start=1):
            v = {}
            ttk.Label(tab, text=c.nome).grid(row=i, column=0, sticky="w", padx=2)
            nat_lab = {v_: k for k, v_ in self.NAT.items()}.get(c.natureza, "Variável")
            v["nat"] = tk.StringVar(value=nat_lab)
            cn = ttk.Combobox(tab, textvariable=v["nat"], values=list(self.NAT), state="readonly", width=12)
            cn.grid(row=i, column=1, padx=2)
            v["cat"] = tk.StringVar(value=c.categoria)
            cc = ttk.Combobox(tab, textvariable=v["cat"], state="readonly", width=62)
            cc.grid(row=i, column=2, padx=2)
            v["cb_cat"] = cc
            v["grupo"] = tk.StringVar(value=c.grupo)
            ttk.Entry(tab, textvariable=v["grupo"], width=10).grid(row=i, column=3, padx=2)
            for j, k in enumerate(("gd", "gf", "gq", "psi0", "psi1", "psi2"), start=4):
                v[k] = tk.StringVar(value=f"{getattr(c, k):g}".replace(".", ","))
                ttk.Entry(tab, textvariable=v[k], width=6).grid(row=i, column=j, padx=2)
            cn.bind("<<ComboboxSelected>>", lambda e, n=c.nome: self._natureza(n, True))
            cc.bind("<<ComboboxSelected>>", lambda e, n=c.nome: self._categoria(n))
            self.linhas[c.nome] = v
            self._natureza(c.nome, False)
        op = ttk.LabelFrame(frm, text="Gerar", padding=6)
        op.pack(fill="x", pady=8)
        self.v_elu, self.v_rara = tk.BooleanVar(value=True), tk.BooleanVar(value=True)
        self.v_freq, self.v_qp = tk.BooleanVar(value=True), tk.BooleanVar(value=True)
        self.v_fav, self.v_opc = tk.BooleanVar(value=True), tk.BooleanVar(value=True)
        self.v_subst = tk.BooleanVar(value=True)
        for j, (t, v) in enumerate((("ELU", self.v_elu), ("ELS rara", self.v_rara), ("ELS frequente", self.v_freq),
                                    ("ELS quase permanente", self.v_qp))):
            ttk.Checkbutton(op, text=t, variable=v).grid(row=0, column=j, sticky="w", padx=6)
        ttk.Checkbutton(op, text="ELU também com permanentes favoráveis (γg,fav) — ex.: sucção de vento",
                        variable=self.v_fav).grid(row=1, column=0, columnspan=4, sticky="w", padx=6)
        ttk.Checkbutton(op, text="Variáveis secundárias opcionais (não considera variável favorável)",
                        variable=self.v_opc).grid(row=2, column=0, columnspan=4, sticky="w", padx=6)
        ttk.Checkbutton(op, text="Substituir combinações geradas anteriormente (as manuais são mantidas)",
                        variable=self.v_subst).grid(row=3, column=0, columnspan=4, sticky="w", padx=6)
        self.txt = tk.Text(frm, height=12, font=("Consolas", 9), wrap="none")
        self.txt.pack(fill="both", expand=True)
        bf = ttk.Frame(frm)
        bf.pack(fill="x", pady=(6, 0))
        ttk.Button(bf, text="Pré-visualizar", command=self._previa).pack(side="left", padx=2)
        ttk.Button(bf, text="Gerar e aplicar", command=self._aplicar).pack(side="right", padx=2)
        ttk.Button(bf, text="Fechar", command=self.destroy).pack(side="right", padx=2)
        cbt.bind("<<ComboboxSelected>>", lambda e: self._recarregar())
        self._previa()
        self.grab_set()
        self.wait_window()

    def _tipo(self):
        return "normal" if self.v_tipo.get() == "Normais" else "especial"

    def _natureza(self, nome, mudou):
        v = self.linhas[nome]
        nat = self.NAT[v["nat"].get()]
        cats = list(nbr.PERMANENTES) if nat == "permanente" else (list(nbr.VARIAVEIS) if nat == "variavel" else [])
        v["cb_cat"].configure(values=cats, state="readonly" if cats else "disabled")
        if mudou:
            v["cat"].set(cats[0] if cats else "")
            self._categoria(nome)

    def _categoria(self, nome):
        v = self.linhas[nome]
        tmp = core.CasoCarga(nome)
        nbr.aplicar_categoria(tmp, self.NAT[v["nat"].get()], v["cat"].get(), self._tipo())
        for k in ("gd", "gf", "gq", "psi0", "psi1", "psi2"):
            v[k].set(f"{getattr(tmp, k):g}".replace(".", ","))

    def _recarregar(self):
        for nome in self.linhas:
            self._categoria(nome)
        self._previa()

    def _grava(self):
        for nome, v in self.linhas.items():
            c = self.modelo.casos[nome]
            c.natureza, c.categoria, c.grupo = self.NAT[v["nat"].get()], v["cat"].get(), v["grupo"].get().strip()
            for k in ("gd", "gf", "gq", "psi0", "psi1", "psi2"):
                setattr(c, k, fnum(v[k].get() or 0))

    def _gera(self):
        self._grava()
        return nbr.gerar(self.modelo, self._tipo(), self.v_elu.get(), self.v_rara.get(), self.v_freq.get(),
                         self.v_qp.get(), self.v_fav.get(), self.v_opc.get())

    def _previa(self):
        self.txt.delete("1.0", "end")
        try:
            lst = self._gera()
        except ValueError as e:
            self.txt.insert("1.0", f"Erro: {e}")
            return
        cont = {}
        for c in lst:
            cont[c.tipo] = cont.get(c.tipo, 0) + 1
        self.txt.insert("end", f"{len(lst)} combinações: " + ", ".join(f"{k}: {n}" for k, n in cont.items()) + "\n\n")
        for c in lst:
            self.txt.insert("end", f"{c.nome:<8} [{c.tipo:<6}] {c.descricao}\n")

    def _aplicar(self):
        try:
            lst = self._gera()
        except ValueError as e:
            messagebox.showerror("Combinações", str(e), parent=self)
            return
        nbr.aplicar(self.modelo, lst, self.v_subst.get())
        self.aplicou = f"{len(lst)} combinações NBR 8800 geradas."
        self.destroy()


class ImportarDXFDialog(tk.Toplevel):
    """Mapeia layers do DXF para seções/ligações e importa as barras."""
    UNID = {"mm": 0.001, "cm": 0.01, "m": 1.0}
    IGN, NOVA = "(ignorar layer)", "(criar nova seção…)"

    def __init__(self, parent, modelo, segs, caminho):
        super().__init__(parent)
        self.title(f"Importar DXF — {caminho.replace(chr(92), '/').split('/')[-1]}")
        self.transient(parent)
        self.parent, self.modelo, self.segs, self.relatorio = parent, modelo, segs, None
        frm = ttk.Frame(self, padding=10)
        frm.pack(fill="both", expand=True)
        esq = ttk.Frame(frm)
        esq.grid(row=0, column=0, sticky="nw")
        self.cv = tk.Canvas(frm, width=420, height=320, bg="white", highlightthickness=1, highlightbackground="#ccc")
        self.cv.grid(row=0, column=1, padx=(12, 0), sticky="n")
        ttk.Label(esq, text="Unidade do desenho").grid(row=0, column=0, sticky="w")
        self.v_un = tk.StringVar(value="mm")
        cu = ttk.Combobox(esq, textvariable=self.v_un, values=list(self.UNID), state="readonly", width=6)
        cu.grid(row=0, column=1, sticky="w")
        cu.bind("<<ComboboxSelected>>", lambda e: self._ext())
        ttk.Label(esq, text="Tolerância de fusão de nós (mm)").grid(row=1, column=0, sticky="w")
        self.v_tol = tk.StringVar(value="1")
        ttk.Entry(esq, textvariable=self.v_tol, width=8).grid(row=1, column=1, sticky="w")
        self.v_queb, self.v_orig = tk.BooleanVar(value=True), tk.BooleanVar(value=True)
        ttk.Checkbutton(esq, text="Quebrar barras nas interseções (banzos contínuos)", variable=self.v_queb).grid(
            row=2, column=0, columnspan=2, sticky="w")
        ttk.Checkbutton(esq, text="Transladar para a origem (menor ponto em 0,0)", variable=self.v_orig).grid(
            row=3, column=0, columnspan=2, sticky="w")
        self.lbl_ext = ttk.Label(esq, text="", foreground="#0b5394")
        self.lbl_ext.grid(row=4, column=0, columnspan=2, sticky="w", pady=(4, 8))
        tab = ttk.Frame(esq)
        tab.grid(row=5, column=0, columnspan=2, sticky="w")
        for j, t in enumerate(("Layer", "Nº", "Seção", "Ligação", "β")):
            ttk.Label(tab, text=t, font=("Segoe UI", 9, "bold")).grid(row=0, column=j, sticky="w", padx=2)
        secs = list(modelo.secoes)
        opcoes = [self.IGN, self.NOVA] + secs
        self.linhas = {}
        for i, (lay, n) in enumerate(dxf.resumo_layers(segs).items(), start=1):
            ttk.Label(tab, text=lay).grid(row=i, column=0, sticky="w", padx=2)
            ttk.Label(tab, text=str(n)).grid(row=i, column=1, sticky="e", padx=2)
            ini = lay if lay in modelo.secoes else (self.IGN if lay.upper() in ("0", "DEFPOINTS") or
                                                    any(k in lay.upper() for k in ("COTA", "TEXT", "DIM", "EIXO", "HATCH"))
                                                    else parent.v_secao.get())
            vs = tk.StringVar(value=ini)
            cb = ttk.Combobox(tab, textvariable=vs, values=opcoes, state="readonly", width=24)
            cb.grid(row=i, column=2, padx=2, pady=1)
            cb.bind("<<ComboboxSelected>>", lambda e: self._desenha())
            trel = any(k in lay.upper() for k in ("DIAG", "MONT", "TREL", "CONTRAV"))
            vl = tk.StringVar(value="Rótula – Rótula (treliça)" if trel else parent.v_lig.get())
            ttk.Combobox(tab, textvariable=vl, values=list(LIGACOES), state="readonly", width=24).grid(row=i, column=3, padx=2)
            vb = tk.StringVar(value="0")
            ttk.Combobox(tab, textvariable=vb, values=["0", "90", "180", "270"], width=5).grid(row=i, column=4, padx=2)
            self.linhas[lay] = (vs, vl, vb)
        bf = ttk.Frame(esq)
        bf.grid(row=6, column=0, columnspan=2, sticky="w", pady=(10, 0))
        ttk.Button(bf, text="Importar", command=self._importar).pack(side="left", padx=2)
        ttk.Button(bf, text="Cancelar", command=self.destroy).pack(side="left", padx=2)
        self._ext()
        self.grab_set()
        self.wait_window()

    def _ext(self):
        x0, y0, x1, y1 = dxf.extensao(self.segs)
        k = self.UNID[self.v_un.get()]
        self.lbl_ext.configure(text=f"Extensão: {fmt((x1 - x0) * k)} × {fmt((y1 - y0) * k)} m   "
                                    f"({len(self.segs)} segmentos)  — confira a unidade!")
        self._desenha()

    def _desenha(self):
        cv = self.cv
        cv.delete("all")
        x0, y0, x1, y1 = dxf.extensao(self.segs)
        W, H = int(cv["width"]), int(cv["height"])
        e = min((W - 20) / max(x1 - x0, 1e-9), (H - 20) / max(y1 - y0, 1e-9))
        cores = {}
        for lay, (vs, _, _) in self.linhas.items():
            s = self.modelo.secoes.get(vs.get())
            cores[lay] = None if vs.get() == self.IGN else (s.cor if s else "#555")
        for lay, a, b, c, d in self.segs:
            cor = cores.get(lay)
            cv.create_line(10 + (a - x0) * e, H - 10 - (b - y0) * e, 10 + (c - x0) * e, H - 10 - (d - y0) * e,
                           fill=cor or "#dddddd", width=2 if cor else 1, dash=() if cor else (2, 2))

    def _importar(self):
        mapa = {}
        try:
            tol = fnum(self.v_tol.get()) / 1000
            assert tol > 0
        except (ValueError, AssertionError):
            messagebox.showerror("DXF", "Tolerância inválida.", parent=self)
            return
        for lay, (vs, vl, vb) in self.linhas.items():
            sec = vs.get()
            if sec == self.IGN:
                mapa[lay] = (None, False, False, 0.0)
                continue
            if sec == self.NOVA:
                cor = PALETA[len(self.modelo.secoes) % len(PALETA)]
                d = SecaoDialog(self, self.modelo, cor=cor, nome=lay if lay not in self.modelo.secoes else None)
                if not d.result:
                    return
                self.modelo.add_secao(d.result)
                sec = d.result.nome
                vs.set(sec)
            try:
                beta = fnum(vb.get() or 0) % 360
            except ValueError:
                beta = 0.0
            mapa[lay] = (sec, *LIGACOES[vl.get()], beta)
        rel = dxf.importar(self.modelo, self.segs, self.UNID[self.v_un.get()], mapa, tol,
                           self.v_queb.get(), self.v_orig.get())
        self.relatorio = rel
        messagebox.showinfo("DXF importado", (
            f"Segmentos lidos: {rel['lidos']}\n"
            f"Ignorados (layer): {rel['ignorados_layer']}\n"
            f"Nulos (< tolerância): {rel['nulos']}\n"
            f"Duplicados (já existiam): {rel['duplicados']}\n\n"
            f"Nós criados: {rel['nos_novos']}\nBarras criadas: {rel['barras_novas']}"), parent=self)
        self.destroy()


# =============================================================================
# APLICAÇÃO
# =============================================================================
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Pórtico 2D — Treliças, vigas e pórticos planos (OpenSees)")
        self.geometry("1400x860")
        self.modelo = self._modelo_inicial()
        self.resultados = None
        self.arquivo = None
        self.sel = []                 # [('no', id) | ('barra', id)]
        self.pendente = None          # nó inicial de barra em lançamento
        self.cursor = None
        self.v_modo = tk.StringVar(value="barra")
        self.v_secao = tk.StringVar(value=next(iter(self.modelo.secoes)))
        self.v_lig = tk.StringVar(value=list(LIGACOES)[0])
        self.v_caso = tk.StringVar(value=next(iter(self.modelo.casos)))
        self.v_res = tk.StringVar()
        self.v_diag = tk.StringVar(value="Estrutura")
        self.v_escala = tk.DoubleVar(value=1.0)
        self.v_valores = tk.BooleanVar(value=True)
        self.v_num_nos = tk.BooleanVar(value=True)
        self.v_num_bar = tk.BooleanVar(value=False)
        self.v_cargas = tk.BooleanVar(value=True)
        self.v_snap = tk.BooleanVar(value=True)
        self.v_eixos = tk.BooleanVar(value=False)
        self.v_extr = tk.BooleanVar(value=False)
        self.v_beta_nova = tk.StringVar(value="0")
        self.arraste = None
        self.verif = {}
        self.visivel = {}
        self._build_menu()
        self._build_ui()
        self._atualiza_listas()
        self.ax.set_xlim(-1, 13)
        self.ax.set_ylim(-2, 7)
        self.redesenha()

    # ------------------------------------------------------------ modelo inicial
    def _modelo_inicial(self):
        m = core.Modelo()
        m.add_secao(core.Secao("W200x15,0", "W", "ASTM A572 Gr50", {"catalogo": "W200x15,0"}, cor=PALETA[0]))
        m.add_secao(core.Secao("Ue 150x60x20x2,65", "Ue", "ASTM A36",
                               {"bw": 150, "bf": 60, "D": 20, "t": 2.65}, cor=PALETA[1]))
        m.add_secao(core.Secao('L 2"x3/16"', "L", "ASTM A572 Gr50", {"b": 50.8, "t": 4.76}, cor=PALETA[2]))
        m.add_caso("G — permanente", peso_proprio=True)
        m.add_caso("Q — sobrecarga")
        m.add_caso("V — vento")
        return m

    # ------------------------------------------------------------ menu / UI
    def _build_menu(self):
        mb = tk.Menu(self)
        ma = tk.Menu(mb, tearoff=0)
        ma.add_command(label="Novo", command=self.novo, accelerator="Ctrl+N")
        ma.add_command(label="Abrir…", command=self.abrir, accelerator="Ctrl+O")
        ma.add_command(label="Salvar", command=self.salvar, accelerator="Ctrl+S")
        ma.add_command(label="Salvar como…", command=lambda: self.salvar(True))
        ma.add_separator()
        ma.add_command(label="Importar barras de DXF…", command=self.importar_dxf)
        ma.add_separator()
        ma.add_command(label="Exportar memorial (.txt)…", command=self.exportar_memorial)
        ma.add_command(label="Exportar figura (.png/.pdf)…", command=self.exportar_figura)
        ma.add_separator()
        ma.add_command(label="Sair", command=self.destroy)
        mb.add_cascade(label="Arquivo", menu=ma)
        me = tk.Menu(mb, tearoff=0)
        me.add_command(label="Nós por coordenadas…", command=self.nos_por_coordenadas)
        me.add_command(label="Barra por nós…", command=self.barra_por_nos)
        me.add_command(label="Mover seleção (dx, dy)…", command=self.mover_selecao)
        me.add_command(label="Dividir barras selecionadas…", command=self.dividir_barras)
        me.add_command(label="Conectar nós soltos às barras", command=lambda: (
            messagebox.showinfo("Conectar", f"{self._conectar_soltos()} divisão(ões) de barra."), self._modificou()))
        me.add_command(label="Propriedades das barras selecionadas…", command=self.propriedades_barras)
        me.add_command(label="Girar perfil das selecionadas +90°", command=self.girar_barras, accelerator="R")
        me.add_command(label="Grade…", command=self.config_grade)
        me.add_separator()
        me.add_command(label="Materiais…", command=self.editar_materiais)
        me.add_command(label="Título do modelo…", command=self.editar_titulo)
        mb.add_cascade(label="Editar", menu=me)
        mv = tk.Menu(mb, tearoff=0)
        mv.add_command(label="Zoom extensão", command=self.zoom_extensao, accelerator="F")
        mv.add_checkbutton(label="Numeração dos nós", variable=self.v_num_nos, command=self.redesenha)
        mv.add_checkbutton(label="Numeração das barras", variable=self.v_num_bar, command=self.redesenha)
        mv.add_checkbutton(label="Mostrar cargas do caso ativo", variable=self.v_cargas, command=self.redesenha)
        mv.add_checkbutton(label="Eixos locais das barras", variable=self.v_eixos, command=self.redesenha)
        mv.add_checkbutton(label="Vista extrudada (seção real)", variable=self.v_extr, command=self.redesenha)
        mb.add_cascade(label="Exibir", menu=mv)
        mh = tk.Menu(mb, tearoff=0)
        mh.add_command(label="Convenções e atalhos", command=self.ajuda)
        mb.add_cascade(label="Ajuda", menu=mh)
        self.config(menu=mb)
        self.bind("<Control-n>", lambda e: self.novo())
        self.bind("<Control-o>", lambda e: self.abrir())
        self.bind("<Control-s>", lambda e: self.salvar())
        self.bind("<Escape>", lambda e: self._cancelar())
        self.bind("<Delete>", lambda e: self.excluir_selecao())
        self.bind("<F5>", lambda e: self.calcular())
        self.bind("f", lambda e: self.zoom_extensao() if not isinstance(e.widget, (tk.Entry, ttk.Entry)) else None)
        self.bind("r", lambda e: self.girar_barras() if not isinstance(e.widget, (tk.Entry, ttk.Entry, tk.Text)) else None)

    def _build_ui(self):
        top = ttk.Frame(self, padding=(4, 4))
        top.pack(side="top", fill="x")
        style = ttk.Style(self)
        style.configure("Modo.Toolbutton", padding=(8, 4))
        for k, lab in MODOS:
            ttk.Radiobutton(top, text=lab, value=k, variable=self.v_modo, style="Modo.Toolbutton",
                            command=self._troca_modo).pack(side="left", padx=1)
        ttk.Separator(top, orient="vertical").pack(side="left", fill="y", padx=6)
        ttk.Button(top, text="Nós por coordenadas", command=self.nos_por_coordenadas).pack(side="left", padx=1)
        ttk.Button(top, text="Barra por nós", command=self.barra_por_nos).pack(side="left", padx=1)
        ttk.Button(top, text="Grade…", command=self.config_grade).pack(side="left", padx=1)
        ttk.Checkbutton(top, text="Snap", variable=self.v_snap).pack(side="left", padx=4)
        ttk.Checkbutton(top, text="Eixos locais", variable=self.v_eixos, command=self.redesenha).pack(side="left", padx=4)
        ttk.Checkbutton(top, text="Extrudado", variable=self.v_extr, command=self.redesenha).pack(side="left", padx=4)
        ttk.Button(top, text="Zoom ext.", command=self.zoom_extensao).pack(side="left", padx=1)
        tk.Button(top, text="▶ CALCULAR (F5)", bg="#0b5394", fg="white", font=("Segoe UI", 10, "bold"),
                  command=self.calcular).pack(side="right", padx=4)

        pw = ttk.PanedWindow(self, orient="horizontal")
        pw.pack(fill="both", expand=True)
        left = ttk.Frame(pw, width=360)
        pw.add(left, weight=0)
        right = ttk.Frame(pw)
        pw.add(right, weight=1)

        self.nb = ttk.Notebook(left)
        self.nb.pack(fill="both", expand=True)
        self._tab_secoes()
        self._tab_camadas()
        self._tab_cargas()
        self._tab_resultados()
        self._tab_dimensionamento()

        self.fig = Figure(figsize=(9, 7), dpi=100)
        self.fig.subplots_adjust(left=0.05, right=0.99, top=0.97, bottom=0.05)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.toolbar = NavigationToolbar2Tk(self.canvas, right, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self.canvas.mpl_connect("motion_notify_event", self._on_move)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("axes_leave_event", lambda e: (self._limpa_hover(), self.canvas.draw_idle()))
        self.canvas.mpl_connect("button_release_event", self._on_release)
        self.canvas.mpl_connect("resize_event", lambda e: self.after_idle(self.redesenha))

        self.status = tk.StringVar(value="")
        ttk.Label(self, textvariable=self.status, anchor="w", relief="sunken", padding=(6, 2)).pack(side="bottom", fill="x")
        self._troca_modo()

    # ---------- aba seções
    def _tab_secoes(self):
        f = ttk.Frame(self.nb, padding=6)
        self.nb.add(f, text="Seções")
        ttk.Label(f, text="Seção ativa (novas barras)").pack(anchor="w")
        self.cb_secao = ttk.Combobox(f, textvariable=self.v_secao, state="readonly")
        self.cb_secao.pack(fill="x", pady=(0, 6))
        ttk.Label(f, text="Ligação das novas barras").pack(anchor="w")
        ttk.Combobox(f, textvariable=self.v_lig, values=list(LIGACOES), state="readonly").pack(fill="x", pady=(0, 4))
        fb = ttk.Frame(f)
        fb.pack(fill="x", pady=(0, 8))
        ttk.Label(fb, text="Giro β das novas barras (°)").pack(side="left")
        ttk.Combobox(fb, textvariable=self.v_beta_nova, values=["0", "90", "180", "270"], width=6).pack(side="left", padx=4)
        self.tv_sec = ttk.Treeview(f, columns=("fam", "mat", "A", "Ix", "Iy"), show="tree headings", height=10)
        for c, t, w in (("#0", "Seção", 125), ("fam", "Fam.", 36), ("mat", "Material", 70), ("A", "A cm²", 48),
                        ("Ix", "Ix cm⁴", 52), ("Iy", "Iy cm⁴", 50)):
            self.tv_sec.heading(c, text=t)
            self.tv_sec.column(c, width=w, anchor="w" if c in ("#0", "mat") else "e")
        self.tv_sec.pack(fill="both", expand=True)
        self.tv_sec.bind("<Double-1>", lambda e: self.editar_secao())
        bf = ttk.Frame(f)
        bf.pack(fill="x", pady=4)
        ttk.Button(bf, text="Nova", command=self.nova_secao).pack(side="left", padx=1)
        ttk.Button(bf, text="Editar", command=self.editar_secao).pack(side="left", padx=1)
        ttk.Button(bf, text="Excluir", command=self.excluir_secao).pack(side="left", padx=1)
        ttk.Button(bf, text="Tornar ativa", command=self._secao_tv_ativa).pack(side="left", padx=1)
        ttk.Separator(f).pack(fill="x", pady=6)
        ttk.Label(f, text="Barras selecionadas (Shift+clique soma):").pack(anchor="w")
        bf2 = ttk.Frame(f)
        bf2.pack(fill="x")
        ttk.Button(bf2, text="Atribuir seção ativa", command=self.atribuir_secao).pack(side="left", padx=1)
        ttk.Button(bf2, text="Atribuir ligação", command=self.atribuir_ligacao).pack(side="left", padx=1)
        ttk.Button(bf2, text="Inverter i↔j", command=self.inverter_barras).pack(side="left", padx=1)
        bf3 = ttk.Frame(f)
        bf3.pack(fill="x", pady=(2, 0))
        ttk.Button(bf3, text="Propriedades…", command=self.propriedades_barras).pack(side="left", padx=1)
        ttk.Button(bf3, text="↻ Girar +90° (R)", command=self.girar_barras).pack(side="left", padx=1)
        ttk.Button(bf3, text="Peso próprio on/off", command=self.alternar_pp).pack(side="left", padx=1)
        self.lbl_info = ttk.Label(f, text="", foreground="#333", wraplength=330, justify="left")
        self.lbl_info.pack(anchor="w", pady=6)

    # ---------- aba camadas
    def _tab_camadas(self):
        f = ttk.Frame(self.nb, padding=6)
        self.nb.add(f, text="Camadas")
        ttk.Label(f, text="Uma camada por perfil. Desmarque para ocultar\n(no modelo e nos resultados). Clique na cor para trocá-la.",
                  foreground="#444").pack(anchor="w")
        bf = ttk.Frame(f)
        bf.pack(fill="x", pady=4)
        ttk.Button(bf, text="Todas", command=lambda: self._camadas_todas(True)).pack(side="left")
        ttk.Button(bf, text="Nenhuma", command=lambda: self._camadas_todas(False)).pack(side="left", padx=4)
        self.frm_camadas = ttk.Frame(f)
        self.frm_camadas.pack(fill="both", expand=True)

    # ---------- aba cargas
    def _tab_cargas(self):
        f = ttk.Frame(self.nb, padding=6)
        self.nb.add(f, text="Carregamentos")
        ttk.Label(f, text="Casos de carga").pack(anchor="w")
        self.lb_casos = tk.Listbox(f, height=6, exportselection=False)
        self.lb_casos.pack(fill="x")
        self.lb_casos.bind("<<ListboxSelect>>", self._caso_sel)
        bf = ttk.Frame(f)
        bf.pack(fill="x", pady=3)
        ttk.Button(bf, text="Novo", command=self.novo_caso).pack(side="left", padx=1)
        ttk.Button(bf, text="Editar", command=self.editar_caso).pack(side="left", padx=1)
        ttk.Button(bf, text="Excluir", command=self.excluir_caso).pack(side="left", padx=1)
        self.lbl_caso = ttk.Label(f, text="", foreground="#0b5394")
        self.lbl_caso.pack(anchor="w", pady=(2, 4))
        ttk.Label(f, text="Cargas do caso ativo").pack(anchor="w")
        self.tv_cargas = ttk.Treeview(f, columns=("alvo", "desc"), show="headings", height=8)
        self.tv_cargas.heading("alvo", text="Onde")
        self.tv_cargas.heading("desc", text="Carga")
        self.tv_cargas.column("alvo", width=70)
        self.tv_cargas.column("desc", width=250)
        self.tv_cargas.pack(fill="both", expand=True)
        ttk.Button(f, text="Excluir carga selecionada", command=self.excluir_carga).pack(anchor="w", pady=3)
        ttk.Separator(f).pack(fill="x", pady=6)
        ttk.Label(f, text="Combinações").pack(anchor="w")
        self.lb_comb = tk.Listbox(f, height=5, exportselection=False)
        self.lb_comb.pack(fill="x")
        self.lb_comb.bind("<Double-1>", lambda e: self.editar_comb())
        bf = ttk.Frame(f)
        bf.pack(fill="x", pady=3)
        ttk.Button(bf, text="Nova", command=self.nova_comb).pack(side="left", padx=1)
        ttk.Button(bf, text="Editar", command=self.editar_comb).pack(side="left", padx=1)
        ttk.Button(bf, text="Excluir", command=self.excluir_comb).pack(side="left", padx=1)
        ttk.Button(bf, text="Excluir geradas", command=self.excluir_comb_auto).pack(side="left", padx=1)
        tk.Button(f, text="⚙ Gerar combinações NBR 8800…", bg="#e8f0fe", command=self.gerar_nbr8800).pack(
            fill="x", pady=(2, 0))

    # ---------- aba resultados
    def _tab_resultados(self):
        f = ttk.Frame(self.nb, padding=6)
        self.nb.add(f, text="Resultados")
        ttk.Label(f, text="Caso / combinação").pack(anchor="w")
        self.cb_res = ttk.Combobox(f, textvariable=self.v_res, state="readonly")
        self.cb_res.pack(fill="x", pady=(0, 6))
        self.cb_res.bind("<<ComboboxSelected>>", lambda e: (self.redesenha(), self._info_resultado()))
        lf = ttk.LabelFrame(f, text="Exibir", padding=4)
        lf.pack(fill="x")
        for d in DIAGRAMAS:
            ttk.Radiobutton(lf, text=d, value=d, variable=self.v_diag, command=self.redesenha).pack(anchor="w")
        ttk.Label(f, text="Escala do diagrama").pack(anchor="w", pady=(6, 0))
        ttk.Scale(f, from_=0.2, to=4.0, variable=self.v_escala, command=lambda v: self.redesenha()).pack(fill="x")
        ttk.Checkbutton(f, text="Mostrar valores", variable=self.v_valores, command=self.redesenha).pack(anchor="w")
        ttk.Label(f, text="Clique numa barra/nó para detalhes:").pack(anchor="w", pady=(8, 0))
        self.txt_res = tk.Text(f, height=16, width=44, font=("Consolas", 9), wrap="none")
        self.txt_res.pack(fill="both", expand=True)
        bf = ttk.Frame(f)
        bf.pack(fill="x", pady=3)
        ttk.Button(bf, text="Memorial completo", command=self.ver_memorial).pack(side="left", padx=1)
        ttk.Button(bf, text="Exportar .txt", command=self.exportar_memorial).pack(side="left", padx=1)

    # ---------- aba dimensionamento
    def _tab_dimensionamento(self):
        f = ttk.Frame(self.nb, padding=6)
        self.nb.add(f, text="Dimensionamento")
        lf = ttk.LabelFrame(f, text="NBR 8800 — perfis I/W (laminados e soldados)", padding=4)
        lf.pack(fill="x")
        c = vf.cfg(self.modelo)
        self.v_ga1 = tk.StringVar(value=fmt(c["ga1"]))
        self.v_ga2 = tk.StringVar(value=fmt(c["ga2"]))
        self.v_halma = tk.StringVar()
        self.v_sigqa = tk.StringVar()
        self.v_celu = tk.StringVar()
        self.v_cels = tk.StringVar()
        r = 0
        for lab, var in (("γa1", self.v_ga1), ("γa2", self.v_ga2)):
            ttk.Label(lf, text=lab).grid(row=r, column=0, sticky="w")
            ttk.Entry(lf, textvariable=var, width=7).grid(row=r, column=1, sticky="w")
            r += 1
        self._OPC = {
            "h_alma": (self.v_halma, {"NBR: d − 2tf − 2r (laminados)": "nbr", "Simplificado: d − 2tf": "simplificado"}),
            "sigma_qa": (self.v_sigqa, {"σ = χ·fy (Q = 1) — F.3.2": "chi", "σ = fy (conservador)": "fy"}),
            "comb_elu": (self.v_celu, {"Combinações ELU": "ELU", "Todas as combinações (exceto ELS)": "todas"}),
            "comb_els": (self.v_cels, {"Todas as ELS": "ELS", "ELS raras": "ELS-R", "ELS frequentes": "ELS-F",
                                       "ELS quase permanentes": "ELS-QP", "Todas as combinações (exceto ELU)": "todas"}),
        }
        for chave, lab in (("h_alma", "Altura h da alma"), ("sigma_qa", "Tensão σ p/ Qa"),
                           ("comb_elu", "Esforços (ELU)"), ("comb_els", "Flecha (ELS)")):
            var, ops = self._OPC[chave]
            ttk.Label(lf, text=lab).grid(row=r, column=0, sticky="w")
            ttk.Combobox(lf, textvariable=var, values=list(ops), state="readonly", width=30).grid(
                row=r, column=1, sticky="w", pady=1)
            r += 1
        self._carrega_cfg_proj()
        for var in (self.v_ga1, self.v_ga2, self.v_halma, self.v_sigqa, self.v_celu, self.v_cels):
            var.trace_add("write", lambda *a: self._grava_cfg_proj())
        bf = ttk.Frame(f)
        bf.pack(fill="x", pady=4)
        ttk.Button(bf, text="Parâmetros das selecionadas…", command=self.parametros_projeto).pack(side="left", padx=1)
        bf2 = ttk.Frame(f)
        bf2.pack(fill="x")
        tk.Button(bf2, text="✓ Verificar selecionadas", bg="#e8f5e9", command=lambda: self.verificar(False)).pack(
            side="left", padx=1)
        tk.Button(bf2, text="✓✓ Verificar todas", bg="#c8e6c9", command=lambda: self.verificar(True)).pack(
            side="left", padx=1)
        ttk.Button(bf2, text="Limpar", command=self.limpar_verificacao).pack(side="left", padx=1)
        self.v_cor_verif = tk.BooleanVar(value=True)
        ttk.Checkbutton(f, text="Colorir barras: verde = passa, vermelho escuro = não passa",
                        variable=self.v_cor_verif, command=self.redesenha).pack(anchor="w", pady=(4, 2))
        self.tv_verif = ttk.Treeview(f, columns=("perfil", "eta", "gov", "st"), show="tree headings", height=12)
        for col, t, w in (("#0", "Barra", 50), ("perfil", "Perfil", 95), ("eta", "η máx %", 60),
                          ("gov", "Governa", 105), ("st", "Estado", 60)):
            self.tv_verif.heading(col, text=t)
            self.tv_verif.column(col, width=w, anchor="w" if col in ("#0", "perfil", "gov") else "e")
        self.tv_verif.tag_configure("OK", foreground="#1b8a2f")
        self.tv_verif.tag_configure("FALHA", foreground="#8b0000")
        self.tv_verif.tag_configure("NA", foreground="#888")
        self.tv_verif.pack(fill="both", expand=True)
        self.tv_verif.bind("<Double-1>", lambda e: self._memorial_tv())
        self.tv_verif.bind("<<TreeviewSelect>>", self._sel_tv_verif)
        bf3 = ttk.Frame(f)
        bf3.pack(fill="x", pady=3)
        ttk.Button(bf3, text="Memorial da barra", command=self._memorial_tv).pack(side="left", padx=1)
        ttk.Button(bf3, text="Exportar todos (HTML/txt)", command=self.exportar_memoriais).pack(side="left", padx=1)
        ttk.Label(f, text="Dica: botão direito sobre uma barra no desenho → verificar / memorial / parâmetros.",
                  foreground="#666", wraplength=330).pack(anchor="w")

    def _carrega_cfg_proj(self):
        c = vf.cfg(self.modelo)
        self._carregando = True
        self.v_ga1.set(fmt(c["ga1"]))
        self.v_ga2.set(fmt(c["ga2"]))
        for chave, (var, ops) in self._OPC.items():
            inv = {v: k for k, v in ops.items()}
            var.set(inv.get(c[chave], list(ops)[0]))
        self._carregando = False

    def _grava_cfg_proj(self):
        if getattr(self, "_carregando", False):
            return
        try:
            self.modelo.projeto["ga1"] = fnum(self.v_ga1.get())
            self.modelo.projeto["ga2"] = fnum(self.v_ga2.get())
        except ValueError:
            return
        for chave, (var, ops) in self._OPC.items():
            if var.get() in ops:
                self.modelo.projeto[chave] = ops[var.get()]
        self.limpar_verificacao(silencioso=True)

    def parametros_projeto(self, ids=None):
        ids = ids or self._barras_sel()
        if not ids:
            messagebox.showinfo("Parâmetros", "Selecione uma ou mais barras (Shift+clique soma).")
            return
        d = ProjetoBarraDialog(self, self.modelo, [self.modelo.barras[i] for i in ids])
        if d.result:
            for i in ids:
                self.modelo.barras[i].proj = dict(d.result)
                self.verif.pop(i, None)
            self._lista_verif()
            self.redesenha()

    def _garantir_analise(self):
        if self.resultados:
            return True
        erros, _ = self.modelo.validar()
        if erros:
            messagebox.showerror("Modelo incompleto", "\n".join(erros))
            return False
        try:
            self.resultados = core.analisar(self.modelo)
        except core.ErroAnalise as e:
            messagebox.showerror("Erro na análise", str(e))
            return False
        nomes = list(self.resultados)
        self.cb_res.configure(values=nomes)
        if self.v_res.get() not in nomes:
            self.v_res.set(nomes[0])
        return True

    def verificar(self, todas, ids=None):
        if ids is None:
            ids = list(self.modelo.barras) if todas else self._barras_sel()
        if not ids:
            messagebox.showinfo("Verificação", "Selecione barras (ou use 'Verificar todas').")
            return
        if not self._garantir_analise():
            return
        self.config(cursor="watch")
        self.update()
        try:
            novos = vf.verificar(self.modelo, self.resultados, ids)
        except Exception as e:
            messagebox.showerror("Verificação", f"{type(e).__name__}: {e}")
            return
        finally:
            self.config(cursor="")
        self.verif.update(novos)
        self._lista_verif()
        self.v_cor_verif.set(True)
        self.v_diag.set("Estrutura")
        self.redesenha()
        ok = sum(1 for v in novos.values() if v.status == "OK")
        fa = sum(1 for v in novos.values() if v.status == "FALHA")
        na = sum(1 for v in novos.values() if v.status == "NA")
        self.status.set(f"Verificadas {len(novos)} barra(s): {ok} passam, {fa} não passam"
                        + (f", {na} não verificáveis (família não implementada)" if na else "") + ".")

    def limpar_verificacao(self, silencioso=False):
        self.verif = {}
        if hasattr(self, "tv_verif"):
            self._lista_verif()
        if not silencioso:
            self.redesenha()

    def _lista_verif(self):
        tv = self.tv_verif
        tv.delete(*tv.get_children())
        for bid in sorted(self.verif):
            v = self.verif[bid]
            if bid not in self.modelo.barras:
                continue
            est = {"OK": "PASSA", "FALHA": "NÃO PASSA", "NA": "N/A"}[v.status]
            tv.insert("", "end", iid=str(bid), text=str(bid), tags=(v.status,),
                      values=(self.modelo.barras[bid].secao, fmt(v.eta_max * 100, 1) if v.status != "NA" else "—",
                              v.governa if v.status != "NA" else v.motivo[:30], est))

    def _sel_tv_verif(self, e=None):
        s = self.tv_verif.selection()
        if s:
            self.sel = [("barra", int(s[0]))]
            self._info_selecao()
            self.redesenha()

    def _memorial_tv(self):
        s = self.tv_verif.selection()
        if s:
            self.memorial_barra(int(s[0]))

    def memorial_barra(self, bid):
        v = self.verif.get(bid)
        if v is None:
            if not messagebox.askyesno("Memorial", f"A barra {bid} ainda não foi verificada. Verificar agora?"):
                return
            self.verificar(False, [bid])
            v = self.verif.get(bid)
            if v is None:
                return
        janela_texto(self, f"Memorial de cálculo — barra {bid}", v.memorial, f"memorial_barra_{bid}.txt",
                     v.html() if v.mem else None)

    def exportar_memoriais(self):
        if not self.verif:
            messagebox.showinfo("Exportar", "Nenhuma barra verificada.")
            return
        p = filedialog.asksaveasfilename(defaultextension=".html", initialfile="memorial_verificacao.html",
                                         filetypes=[("Página HTML (imprimível em PDF)", "*.html"), ("Texto", "*.txt")])
        if not p:
            return
        if p.lower().endswith((".html", ".htm")):
            with open(p, "w", encoding="utf-8") as f:
                f.write(vf.memoriais_html(self.modelo, self.verif))
            self.status.set(f"Memoriais exportados: {p}")
            return
        partes = [f"MEMORIAL DE VERIFICAÇÃO — {self.modelo.titulo}", "",
                  f"{'Barra':>6}  {'Perfil':<22}{'η máx (%)':>10}  {'Governa':<22}Estado"]
        for bid in sorted(self.verif):
            v = self.verif[bid]
            partes.append(f"{bid:>6}  {self.modelo.barras[bid].secao:<22}{fmt(v.eta_max * 100, 1):>10}  "
                          f"{v.governa:<22}{ {'OK': 'PASSA', 'FALHA': 'NÃO PASSA', 'NA': 'N/A'}[v.status]}")
        partes.append("")
        for bid in sorted(self.verif):
            partes += [self.verif[bid].memorial, "", ""]
        with open(p, "w", encoding="utf-8") as f:
            f.write("\n".join(partes))
        self.status.set(f"Memoriais exportados: {p}")

    def _menu_barra(self, ev, b):
        m = tk.Menu(self, tearoff=0)
        m.add_command(label=f"Barra {b.id} — {b.secao}", state="disabled")
        m.add_separator()
        m.add_command(label="✓ Verificar esta barra (NBR 8800)", command=lambda: self.verificar(False, [b.id]))
        m.add_command(label="📄 Memorial de cálculo", command=lambda: self.memorial_barra(b.id))
        m.add_command(label="Parâmetros de dimensionamento…", command=lambda: self.parametros_projeto([b.id]))
        m.add_separator()
        m.add_command(label="Propriedades da barra…", command=lambda: self._aplica_dialogo_barras([b]))
        x, y = self.canvas.get_tk_widget().winfo_rootx(), self.canvas.get_tk_widget().winfo_rooty()
        h = self.canvas.get_tk_widget().winfo_height()
        try:
            m.tk_popup(int(x + ev.x), int(y + h - ev.y))
        finally:
            m.grab_release()

    # ============================================================ listas
    def _atualiza_listas(self):
        m = self.modelo
        nomes = list(m.secoes)
        self.cb_secao.configure(values=nomes)
        if self.v_secao.get() not in m.secoes:
            self.v_secao.set(nomes[0] if nomes else "")
        self.tv_sec.delete(*self.tv_sec.get_children())
        for s in m.secoes.values():
            p = s.propriedades()
            tag = "t_" + s.nome
            self.tv_sec.tag_configure(tag, foreground=s.cor)
            self.tv_sec.insert("", "end", iid=s.nome, text="■ " + s.nome, tags=(tag,),
                               values=(s.familia, s.material.replace("ASTM ", ""), f"{p['A']:.2f}", f"{p['Ix']:.0f}", f"{p['Iy']:.0f}"))
        # camadas
        for w in self.frm_camadas.winfo_children():
            w.destroy()
        for s in m.secoes.values():
            v = self.visivel.setdefault(s.nome, tk.BooleanVar(value=True))
            row = ttk.Frame(self.frm_camadas)
            row.pack(fill="x", pady=1)
            tk.Button(row, text="   ", bg=s.cor, relief="groove",
                      command=lambda n=s.nome: self._cor_camada(n)).pack(side="left")
            nb = sum(1 for b in m.barras.values() if b.secao == s.nome)
            ttk.Checkbutton(row, text=f"{s.nome}  ({nb} barras)", variable=v, command=self.redesenha).pack(side="left", padx=4)
        # casos
        self.lb_casos.delete(0, "end")
        for c in m.casos.values():
            self.lb_casos.insert("end", c.nome + ("   [+PP]" if c.peso_proprio else ""))
        if self.v_caso.get() not in m.casos and m.casos:
            self.v_caso.set(next(iter(m.casos)))
        if self.v_caso.get() in m.casos:
            i = list(m.casos).index(self.v_caso.get())
            self.lb_casos.selection_clear(0, "end")
            self.lb_casos.selection_set(i)
        self.lbl_caso.configure(text=f"Caso ativo para lançar cargas: {self.v_caso.get()}")
        self.lb_comb.delete(0, "end")
        for c in m.combinacoes.values():
            self.lb_comb.insert("end", (f"[{c.tipo}] " if c.tipo else "") + f"{c.nome}: " +
                                " + ".join(f"{f:g}·{k.split(' ')[0]}" for k, f in c.fatores.items() if f))
        self._lista_cargas()

    def _lista_cargas(self):
        self.tv_cargas.delete(*self.tv_cargas.get_children())
        caso = self.v_caso.get()
        for i, c in enumerate(self.modelo.cargas_no):
            if c.caso == caso:
                self.tv_cargas.insert("", "end", iid=f"n{i}", values=(f"Nó {c.no}", f"Fx={c.fx:g} Fy={c.fy:g} kN  Mz={c.mz:g} kN·m"))
        for i, c in enumerate(self.modelo.cargas_barra):
            if c.caso == caso:
                if c.tipo == "distribuida":
                    d = f"q={c.valor:g} kN/m  {core.DIRECOES[c.direcao]}" + (" proj." if c.projetada else "")
                else:
                    d = f"P={c.valor:g} kN  {core.DIRECOES[c.direcao]}  a={c.pos:g} m"
                self.tv_cargas.insert("", "end", iid=f"b{i}", values=(f"Barra {c.barra}", d))

    def _cor_camada(self, nome):
        s = self.modelo.secoes[nome]
        c = colorchooser.askcolor(s.cor, parent=self)[1]
        if c:
            s.cor = c
            self._atualiza_listas()
            self.redesenha()

    def _camadas_todas(self, v):
        for var in self.visivel.values():
            var.set(v)
        self.redesenha()

    def barra_visivel(self, b):
        v = self.visivel.get(b.secao)
        return v.get() if v else True

    # ============================================================ geometria da tela
    def _span(self):
        x0, x1 = self.ax.get_xlim()
        y0, y1 = self.ax.get_ylim()
        return max(x1 - x0, y1 - y0, 1e-6)

    def _snap(self, x, y):
        if not self.v_snap.get():
            return x, y
        dx, dy = self.modelo.grade["dx"], self.modelo.grade["dy"]
        return round(x / dx) * dx, round(y / dy) * dy

    def _pix(self, x, y):
        return self.ax.transData.transform((x, y))

    def _no_perto(self, ev, tol=10):
        best, dmin = None, tol
        for n in self.modelo.nos.values():
            px, py = self._pix(n.x, n.y)
            d = math.hypot(px - ev.x, py - ev.y)
            if d < dmin:
                best, dmin = n, d
        return best

    def _barra_perto(self, ev, tol=8):
        best, dmin = None, tol
        for b in self.modelo.barras.values():
            if not self.barra_visivel(b):
                continue
            n1, n2 = self.modelo.nos[b.ni], self.modelo.nos[b.nj]
            x1, y1 = self._pix(n1.x, n1.y)
            x2, y2 = self._pix(n2.x, n2.y)
            vx, vy = x2 - x1, y2 - y1
            L2 = vx * vx + vy * vy
            if L2 == 0:
                continue
            t = max(0, min(1, ((ev.x - x1) * vx + (ev.y - y1) * vy) / L2))
            d = math.hypot(ev.x - (x1 + t * vx), ev.y - (y1 + t * vy))
            if d < dmin:
                best, dmin = b, d
        return best

    # ============================================================ eventos
    def _troca_modo(self):
        self.pendente = None
        dicas = {"sel": "Clique: seleciona | Shift+clique: soma à seleção | Duplo clique: editar | Del: excluir",
                 "no": "Clique no grid para criar nó (snap na grade). Coordenadas exatas: 'Nós por coordenadas'.",
                 "barra": "Clique no nó inicial e depois no final (cria nós no grid se preciso). Lançamento contínuo; Esc encerra.",
                 "apoio": "Clique num nó para definir/editar o apoio.",
                 "cno": "Clique num nó para lançar carga no caso ativo.",
                 "cbarra": "Clique numa barra para lançar carga distribuída ou pontual no caso ativo.",
                 "del": "Clique num nó (remove barras ligadas) ou numa barra para excluir.",
                 "mover": "Arraste um nó ou uma barra (se estiverem selecionados, move toda a seleção). "
                          "Deslocamento exato: Editar › Mover seleção."}
        self.dica = dicas.get(self.v_modo.get(), "")
        self.status.set(self.dica)
        if self.v_modo.get() != "sel" and self.v_diag.get() != "Estrutura":
            self.v_diag.set("Estrutura")
        self.redesenha()

    def _cancelar(self):
        self.pendente = None
        self.sel = []
        self.redesenha()

    def _on_move(self, ev):
        if ev.inaxes != self.ax or ev.xdata is None:
            return
        x, y = self._snap(ev.xdata, ev.ydata)
        if self.arraste:
            a = self.arraste
            dx, dy = x - a["p0"][0], y - a["p0"][1]
            nos = self.modelo.nos
            antes = {nid: (nos[nid].x, nos[nid].y) for nid in a["orig"]}
            for nid, (ox, oy) in a["orig"].items():
                nos[nid].x, nos[nid].y = round(ox + dx, 9), round(oy + dy, 9)
            if self._conflito_nos(set(a["orig"])):
                for nid, (px, py) in antes.items():
                    nos[nid].x, nos[nid].y = px, py
                self.status.set("Posição ocupada por outro nó — escolha outro ponto.")
                return
            self.status.set(f"Movendo: dx = {fmt(dx, 3)} m   dy = {fmt(dy, 3)} m")
            self.redesenha()
            return
        self.status.set(f"x = {fmt(x, 3)} m   y = {fmt(y, 3)} m      |  {self.dica}")
        self._hover_resultado(ev)
        if self.pendente and self.v_modo.get() == "barra":
            n = self._no_perto(ev)
            self.cursor = (n.x, n.y) if n else (x, y)
            if hasattr(self, "_linha_prev") and self._linha_prev is not None:
                p = self.modelo.nos[self.pendente]
                self._linha_prev.set_data([p.x, self.cursor[0]], [p.y, self.cursor[1]])
                self.canvas.draw_idle()

    def _limpa_hover(self):
        for a in getattr(self, "_hover_art", []):
            try:
                a.remove()
            except (ValueError, NotImplementedError):
                pass
        self._hover_art = []

    def _hover_resultado(self, ev):
        """Leitura contínua do diagrama/deformada no ponto da barra mais próximo do cursor."""
        if not getattr(self, "_hover", None) or not self._hover_meta or self.arraste:
            if getattr(self, "_hover_art", None):
                self._limpa_hover()
                self.canvas.draw_idle()
            return
        self._limpa_hover()
        cur = np.array([ev.x, ev.y])
        melhor = None
        for bid, (xs, base, curva) in self._hover.items():
            for pts in (curva, base):
                P = self.ax.transData.transform(np.asarray(pts))
                A, B = P[:-1], P[1:]
                AB = B - A
                L2 = (AB ** 2).sum(axis=1)
                L2[L2 == 0] = 1e-12
                t = np.clip(((cur - A) * AB).sum(axis=1) / L2, 0, 1)
                D = np.hypot(*(A + AB * t[:, None] - cur).T)
                k = int(D.argmin())
                if melhor is None or D[k] < melhor[0]:
                    melhor = (D[k], bid, xs[k] + t[k] * (xs[k + 1] - xs[k]))
        if melhor is None or melhor[0] > 14:
            self.canvas.draw_idle()
            return
        _, bid, xl = melhor
        m, r, meta = self.modelo, self._res_atual(), self._hover_meta
        b = m.barras[bid]
        n1 = m.nos[b.ni]
        L, c, s = m.geometria(b)
        bx, by = n1.x + c * xl, n1.y + s * xl
        if meta["tipo"] == "desloc":
            _, DX, DY = core.desloc_global_barra(m, r, bid, [xl])
            dx, dy = DX[0], DY[0]
            px, py = bx + dx * meta["esc"], by + dy * meta["esc"]
            _, U, V = core.deformada_barra(m, r, bid, xs=[xl])
            corda = r.ext[bid]
            v1 = -corda[0][0] * s + corda[0][1] * c
            v2 = -corda[1][0] * s + corda[1][1] * c
            flecha = V[0] - (v1 + (v2 - v1) * xl / L)
            txt = (f"Barra {bid}   x = {fmt(xl, 3)} m\n"
                   f"δx = {fmt(dx * 1e3)} mm   δy = {fmt(dy * 1e3)} mm\n"
                   f"|δ| = {fmt(math.hypot(dx, dy) * 1e3)} mm\n"
                   f"flecha rel. à corda = {fmt(flecha * 1e3)} mm")
            cor = "#b71c1c"
        else:
            xs, N, V, M = core.esforcos_barra(m, r, bid, xs=[xl])
            val = (N, V, M)[meta["idx"] - 1][0]
            off = val * meta["esc"] * meta["sinal"]
            px, py = bx - s * off, by + c * off
            txt = (f"Barra {bid}   x = {fmt(xl, 3)} m\n"
                   f"{meta['nome']} = {fmt(val)} {meta['un']}")
            if meta["idx"] == 1:
                txt += "  (tração)" if val > 1e-9 else ("  (compressão)" if val < -1e-9 else "")
            cor = "#222"
        ax = self.ax
        self._hover_art += ax.plot([bx, px], [by, py], color="#ff9800", lw=1.4, zorder=20)
        self._hover_art += ax.plot([px], [py], "o", ms=7, mfc="#ff9800", mec="k", mew=0.8, zorder=21)
        self._hover_art += ax.plot([bx], [by], "o", ms=4, color="#ff9800", zorder=21)
        # caixa de texto deslocada do cursor, para dentro da área visível
        W = ax.get_window_extent()
        dxp = -15 if ev.x > W.x0 + 0.7 * W.width else 15
        dyp = -15 if ev.y > W.y0 + 0.7 * W.height else 15
        self._hover_art.append(ax.annotate(
            txt, xy=(px, py), xytext=(dxp, dyp), textcoords="offset points",
            ha="right" if dxp < 0 else "left", va="top" if dyp < 0 else "bottom",
            fontsize=8.5, color=cor, zorder=22, family="monospace",
            bbox=dict(boxstyle="round,pad=0.35", fc="#fffbe6", ec="#ff9800", lw=0.8)))
        self.canvas.draw_idle()

    def _on_scroll(self, ev):
        if ev.inaxes != self.ax:
            return
        f = 0.8 if ev.button == "up" else 1.25
        x0, x1 = self.ax.get_xlim()
        y0, y1 = self.ax.get_ylim()
        cx, cy = ev.xdata, ev.ydata
        self.ax.set_xlim(cx + (x0 - cx) * f, cx + (x1 - cx) * f)
        self.ax.set_ylim(cy + (y0 - cy) * f, cy + (y1 - cy) * f)
        self.redesenha()

    def _on_release(self, ev):
        if self.arraste:
            a, self.arraste = self.arraste, None
            m = self.modelo
            conflito = self._conflito_nos(set(a["orig"]))
            if conflito:
                for nid, (ox, oy) in a["orig"].items():
                    m.nos[nid].x, m.nos[nid].y = ox, oy
                messagebox.showwarning("Mover", conflito)
            moveu = any((m.nos[n].x, m.nos[n].y) != o for n, o in a["orig"].items())
            if moveu:
                self._conectar_soltos()
                self._modificou()
            else:
                self.redesenha()
            return
        if self.toolbar.mode:
            self.redesenha()

    def _conflito_nos(self, movidos):
        """Mensagem de erro se algum nó movido coincide com outro nó ou gera barra nula."""
        m = self.modelo
        for nid in movidos:
            n = m.nos[nid]
            for o in m.nos.values():
                if o.id != nid and abs(o.x - n.x) < 1e-6 and abs(o.y - n.y) < 1e-6:
                    return f"O nó {nid} ficaria sobre o nó {o.id}. Movimento desfeito."
        return None

    def _nos_da_selecao(self):
        ids = set()
        for t, i in self.sel:
            if t == "no" and i in self.modelo.nos:
                ids.add(i)
            elif t == "barra" and i in self.modelo.barras:
                ids |= {self.modelo.barras[i].ni, self.modelo.barras[i].nj}
        return ids

    def _conectar_soltos(self, avisar=True):
        """Nó sobre o trecho interno de uma barra → divide a barra no nó (como no Ftool)."""
        k = self.modelo.inserir_nos_nas_barras()
        if k and avisar:
            self.status.set(f"{k} barra(s) dividida(s) em nós sobre elas.")
        return k

    def dividir_barras(self):
        ids = self._barras_sel()
        if not ids:
            messagebox.showinfo("Dividir", "Selecione uma ou mais barras.")
            return
        d = FormDialog(self, "Dividir barras", [dict(key="n", label="Número de trechos iguais", type="int", default=2)],
                       texto="Cria nós intermediários igualmente espaçados. Cargas distribuídas são copiadas "
                             "e pontuais redistribuídas pela posição.")
        if d.result and d.result["n"] >= 2:
            for i in ids:
                self.modelo.dividir_em_partes(i, d.result["n"])
            self.sel = []
            self._modificou()

    def mover_selecao(self):
        ids = self._nos_da_selecao()
        if not ids:
            messagebox.showinfo("Mover", "Selecione nós e/ou barras (modo Selecionar, Shift+clique soma).")
            return
        d = FormDialog(self, "Mover seleção", [dict(key="dx", label="dx (m)", default=0.0),
                                               dict(key="dy", label="dy (m)", default=0.0)],
                       texto=f"{len(ids)} nó(s) serão deslocados. As barras ligadas acompanham.")
        if not d.result:
            return
        m = self.modelo
        orig = {i: (m.nos[i].x, m.nos[i].y) for i in ids}
        for i in ids:
            m.nos[i].x = round(orig[i][0] + d.result["dx"], 9)
            m.nos[i].y = round(orig[i][1] + d.result["dy"], 9)
        c = self._conflito_nos(ids)
        if c:
            for i, (ox, oy) in orig.items():
                m.nos[i].x, m.nos[i].y = ox, oy
            messagebox.showwarning("Mover", c)
        self._modificou()

    def _on_click(self, ev):
        if self.toolbar.mode or ev.inaxes != self.ax or ev.xdata is None:
            return
        if ev.button == 3:
            b = None if self.pendente else self._barra_perto(ev)
            if b and not self._no_perto(ev):
                self.sel = [("barra", b.id)]
                self._info_selecao()
                self.redesenha()
                self._menu_barra(ev, b)
                return
            self._cancelar()
            return
        modo = self.v_modo.get()
        m = self.modelo
        x, y = self._snap(ev.xdata, ev.ydata)
        if modo == "no":
            m.add_no(x, y)
            self._conectar_soltos()
            self._modificou()
        elif modo == "barra":
            n = self._no_perto(ev) or m.add_no(x, y)
            if self.pendente is None:
                self.pendente = n.id
            elif n.id != self.pendente:
                ri, rj = LIGACOES[self.v_lig.get()]
                try:
                    nb = m.add_barra(self.pendente, n.id, self.v_secao.get(), ri, rj)
                    try:
                        nb.beta = fnum(self.v_beta_nova.get() or 0) % 360
                    except ValueError:
                        pass
                except ValueError as e:
                    self.status.set(str(e))
                self.pendente = n.id
            self._conectar_soltos()
            self._modificou()
        elif modo == "mover":
            n = self._no_perto(ev)
            item = ("no", n.id) if n else None
            if not item:
                b = self._barra_perto(ev)
                item = ("barra", b.id) if b else None
            if not item:
                return
            if item not in self.sel:
                self.sel = [item]
            ids = self._nos_da_selecao()
            p0 = (n.x, n.y) if n else (x, y)
            self.arraste = {"orig": {i: (m.nos[i].x, m.nos[i].y) for i in ids}, "p0": p0}
        elif modo == "apoio":
            n = self._no_perto(ev)
            if n:
                d = ApoioDialog(self, n.id, m.apoios.get(n.id))
                if d.remover:
                    m.apoios.pop(n.id, None)
                elif d.result:
                    m.apoios[n.id] = d.result
                self._modificou()
        elif modo == "cno":
            n = self._no_perto(ev)
            if n:
                self.dialogo_carga_no(n.id)
        elif modo == "cbarra":
            b = self._barra_perto(ev)
            if b:
                self.dialogo_carga_barra(b.id, ev)
        elif modo == "del":
            n = self._no_perto(ev)
            if n:
                m.remove_no(n.id)
            else:
                b = self._barra_perto(ev)
                if b:
                    m.remove_barra(b.id)
            self._modificou()
        elif modo == "sel":
            n = self._no_perto(ev)
            item = ("no", n.id) if n else None
            if not item:
                b = self._barra_perto(ev)
                item = ("barra", b.id) if b else None
            if ev.dblclick and item:
                self.sel = [item]
                self.editar_item(item)
                return
            if ev.key == "shift" and item:
                if item in self.sel:
                    self.sel.remove(item)
                else:
                    self.sel.append(item)
            else:
                self.sel = [item] if item else []
            self._info_selecao()
            self.redesenha()

    # ============================================================ edição
    def _modificou(self):
        self.verif = {}
        if hasattr(self, "tv_verif"):
            self._lista_verif()
        self.resultados = None
        self.cb_res.configure(values=[])
        self.v_res.set("")
        self._atualiza_listas()
        self.redesenha()

    def _info_selecao(self):
        m = self.modelo
        lin = []
        for tipo, i in self.sel[:6]:
            if tipo == "no" and i in m.nos:
                n = m.nos[i]
                lin.append(f"Nó {i}: ({fmt(n.x, 3)}; {fmt(n.y, 3)}) m" + (f" — apoio {m.apoios[i].descricao()}" if i in m.apoios else ""))
            elif tipo == "barra" and i in m.barras:
                b = m.barras[i]
                L = m.geometria(b)[0]
                lin.append(f"Barra {i}: nós {b.ni}→{b.nj}, L={fmt(L, 3)} m, β={b.beta:g}°"
                           f"{'' if b.pp else ', SEM peso próprio'}\n  {b.secao} | {b.tipo}")
        if len(self.sel) > 6:
            lin.append(f"... +{len(self.sel) - 6}")
        self.lbl_info.configure(text="\n".join(lin))
        if self.resultados:
            self._info_resultado()

    def editar_item(self, item):
        m = self.modelo
        tipo, i = item
        if tipo == "no":
            n = m.nos[i]
            d = FormDialog(self, f"Nó {i}", [dict(key="x", label="X (m)", default=n.x),
                                              dict(key="y", label="Y (m)", default=n.y)])
            if d.result:
                o = m.no_em(d.result["x"], d.result["y"])
                if o and o.id != i:
                    messagebox.showerror("Nó", f"Já existe o nó {o.id} nessa posição.")
                    return
                n.x, n.y = d.result["x"], d.result["y"]
        else:
            self._aplica_dialogo_barras([m.barras[i]])
            return
        self._modificou()

    def _aplica_dialogo_barras(self, barras):
        d = BarraDialog(self, self.modelo, barras)
        if d.result:
            for b in barras:
                b.secao = d.result["secao"]
                b.rot_i, b.rot_j = d.result["lig"]
                b.beta = d.result["beta"]
                b.pp = d.result["pp"]
            self._modificou()
            self._info_selecao()

    def propriedades_barras(self):
        ids = self._barras_sel()
        if not ids:
            messagebox.showinfo("Barras", "Selecione uma ou mais barras.")
            return
        self._aplica_dialogo_barras([self.modelo.barras[i] for i in ids])

    def girar_barras(self):
        ids = self._barras_sel()
        for i in ids:
            b = self.modelo.barras[i]
            b.beta = (b.beta + 90.0) % 360
        if ids:
            self._modificou()
            self._info_selecao()
            self.status.set(f"Perfil girado +90° em {len(ids)} barra(s).")

    def alternar_pp(self):
        ids = self._barras_sel()
        if not ids:
            return
        novo = not all(self.modelo.barras[i].pp for i in ids)
        for i in ids:
            self.modelo.barras[i].pp = novo
        self._modificou()
        self._info_selecao()
        self.status.set(f"Peso próprio {'LIGADO' if novo else 'DESLIGADO'} em {len(ids)} barra(s).")

    def _barras_sel(self):
        return [i for t, i in self.sel if t == "barra" and i in self.modelo.barras]

    def atribuir_secao(self):
        for i in self._barras_sel():
            self.modelo.barras[i].secao = self.v_secao.get()
        self._modificou()

    def atribuir_ligacao(self):
        for i in self._barras_sel():
            self.modelo.barras[i].rot_i, self.modelo.barras[i].rot_j = LIGACOES[self.v_lig.get()]
        self._modificou()

    def inverter_barras(self):
        for i in self._barras_sel():
            b = self.modelo.barras[i]
            b.ni, b.nj, b.rot_i, b.rot_j = b.nj, b.ni, b.rot_j, b.rot_i
            for c in self.modelo.cargas_barra:
                if c.barra == i and c.tipo == "pontual":
                    c.pos = self.modelo.geometria(b)[0] - c.pos
        self._modificou()

    def excluir_selecao(self):
        if isinstance(self.focus_get(), (tk.Entry, ttk.Entry, tk.Text)):
            return
        for t, i in self.sel:
            if t == "barra":
                self.modelo.remove_barra(i)
        for t, i in self.sel:
            if t == "no":
                self.modelo.remove_no(i)
        self.sel = []
        self._modificou()

    def nos_por_coordenadas(self):
        w = tk.Toplevel(self)
        w.title("Nós por coordenadas")
        w.transient(self)
        ttk.Label(w, text="Um nó por linha:  X  Y   (m, separados por espaço ou ';')\n"
                          "Ex.:  0 0\n       6 0\n       3;1,5", padding=8, justify="left").pack(anchor="w")
        t = tk.Text(w, width=36, height=12, font=("Consolas", 10))
        t.pack(padx=8)
        t.focus_set()

        def ok():
            criados = 0
            for ln in t.get("1.0", "end").splitlines():
                ln = ln.strip().replace(";", " ").replace("\t", " ")
                if not ln:
                    continue
                p = ln.split()
                try:
                    x, y = fnum(p[0]), fnum(p[1])
                except (ValueError, IndexError):
                    messagebox.showerror("Nós", f"Linha inválida: '{ln}'", parent=w)
                    return
                n0 = len(self.modelo.nos)
                self.modelo.add_no(x, y)
                criados += len(self.modelo.nos) - n0
            w.destroy()
            self._conectar_soltos()
            self._modificou()
            self.zoom_extensao()
            self.status.set(f"{criados} nó(s) criado(s).")
        ttk.Button(w, text="Criar nós", command=ok).pack(pady=8)
        w.grab_set()

    def barra_por_nos(self):
        d = FormDialog(self, "Barra por nós", [
            dict(key="ni", label="Nó inicial (i)", type="int", default=""),
            dict(key="nj", label="Nó final (j)", type="int", default=""),
            dict(key="secao", label="Seção", type="combo", options=list(self.modelo.secoes), default=self.v_secao.get()),
            dict(key="lig", label="Ligação", type="combo", options=list(LIGACOES), default=self.v_lig.get())])
        if d.result:
            r = d.result
            if r["ni"] not in self.modelo.nos or r["nj"] not in self.modelo.nos:
                messagebox.showerror("Barra", "Nó inexistente.")
                return
            try:
                self.modelo.add_barra(r["ni"], r["nj"], r["secao"], *LIGACOES[r["lig"]])
            except ValueError as e:
                messagebox.showerror("Barra", str(e))
            self._modificou()

    def config_grade(self):
        g = self.modelo.grade
        d = FormDialog(self, "Grade", [dict(key="dx", label="Espaçamento X (m)", default=g["dx"]),
                                       dict(key="dy", label="Espaçamento Y (m)", default=g["dy"])])
        if d.result and d.result["dx"] > 0 and d.result["dy"] > 0:
            g.update(d.result)
            self.redesenha()

    def editar_materiais(self):
        m = self.modelo
        nome = FormDialog(self, "Materiais", [dict(key="mat", label="Material", type="combo",
                                                   options=list(m.materiais), default="ASTM A36")])
        if not nome.result:
            return
        mt = m.materiais[nome.result["mat"]]
        d = FormDialog(self, mt.nome, [dict(key="E", label="E (MPa)", default=mt.E), dict(key="G", label="G (MPa)", default=mt.G),
                                       dict(key="fy", label="fy (MPa)", default=mt.fy), dict(key="fu", label="fu (MPa)", default=mt.fu),
                                       dict(key="gamma", label="γ (kN/m³)", default=mt.gamma)])
        if d.result:
            for k, v in d.result.items():
                setattr(mt, k, v)
            self._modificou()

    def editar_titulo(self):
        d = FormDialog(self, "Título", [dict(key="t", label="Título", type="str", default=self.modelo.titulo)])
        if d.result:
            self.modelo.titulo = d.result["t"]

    # ---------- seções
    def nova_secao(self):
        cor = PALETA[len(self.modelo.secoes) % len(PALETA)]
        d = SecaoDialog(self, self.modelo, cor=cor)
        if d.result:
            self.modelo.add_secao(d.result)
            self.v_secao.set(d.result.nome)
            self._modificou()

    def _tv_secao(self):
        s = self.tv_sec.selection()
        return s[0] if s else None

    def editar_secao(self):
        nome = self._tv_secao()
        if not nome:
            return
        d = SecaoDialog(self, self.modelo, self.modelo.secoes[nome])
        if d.result:
            s = d.result
            self.modelo.secoes.pop(nome)
            self.modelo.secoes[s.nome] = s
            if s.nome != nome:
                for b in self.modelo.barras.values():
                    if b.secao == nome:
                        b.secao = s.nome
                if nome in self.visivel:
                    self.visivel[s.nome] = self.visivel.pop(nome)
                if self.v_secao.get() == nome:
                    self.v_secao.set(s.nome)
            self._modificou()

    def excluir_secao(self):
        nome = self._tv_secao()
        if nome:
            try:
                self.modelo.remove_secao(nome)
            except ValueError as e:
                messagebox.showerror("Seção", str(e))
            self._modificou()

    def _secao_tv_ativa(self):
        nome = self._tv_secao()
        if nome:
            self.v_secao.set(nome)

    # ---------- casos / cargas / combinações
    def _caso_sel(self, e=None):
        s = self.lb_casos.curselection()
        if s:
            self.v_caso.set(list(self.modelo.casos)[s[0]])
            self.lbl_caso.configure(text=f"Caso ativo para lançar cargas: {self.v_caso.get()}")
            self._lista_cargas()
            self.redesenha()

    def novo_caso(self):
        d = FormDialog(self, "Novo caso de carga", [
            dict(key="nome", label="Nome", type="str", default=f"Caso {len(self.modelo.casos) + 1}"),
            dict(key="desc", label="Descrição", type="str", default=""),
            dict(key="pp", label="Incluir peso próprio automático das barras", type="check", default=False)])
        if d.result and d.result["nome"].strip():
            if d.result["nome"] in self.modelo.casos:
                messagebox.showerror("Caso", "Nome já existe.")
                return
            self.modelo.add_caso(d.result["nome"].strip(), d.result["pp"], d.result["desc"])
            self.v_caso.set(d.result["nome"].strip())
            self._modificou()

    def editar_caso(self):
        c = self.modelo.casos.get(self.v_caso.get())
        if not c:
            return
        d = FormDialog(self, "Editar caso", [
            dict(key="nome", label="Nome", type="str", default=c.nome),
            dict(key="desc", label="Descrição", type="str", default=c.descricao),
            dict(key="pp", label="Incluir peso próprio automático das barras", type="check", default=c.peso_proprio)])
        if d.result:
            novo = d.result["nome"].strip()
            if novo and novo != c.nome:
                try:
                    self.modelo.renomear_caso(c.nome, novo)
                except ValueError as e:
                    messagebox.showerror("Caso", str(e))
                    return
                self.v_caso.set(novo)
            c.peso_proprio, c.descricao = d.result["pp"], d.result["desc"]
            self._modificou()

    def excluir_caso(self):
        c = self.v_caso.get()
        if c and messagebox.askyesno("Excluir caso", f"Excluir o caso '{c}' e todas as suas cargas?"):
            self.modelo.remove_caso(c)
            self.v_caso.set("")
            self._modificou()

    def _exige_caso(self):
        if self.v_caso.get() not in self.modelo.casos:
            messagebox.showwarning("Carga", "Crie/selecione um caso de carga na aba Carregamentos.")
            return False
        return True

    def dialogo_carga_no(self, nid):
        if not self._exige_caso():
            return
        d = FormDialog(self, f"Carga no nó {nid}", [
            dict(key="fx", label="Fx (kN)  → +", default=0.0),
            dict(key="fy", label="Fy (kN)  ↑ +", default=0.0),
            dict(key="mz", label="Mz (kN·m)  ↺ +", default=0.0)],
            texto=f"Caso: {self.v_caso.get()}\nEixos globais. Carga para baixo → Fy negativo.")
        if d.result and any(d.result.values()):
            self.modelo.cargas_no.append(core.CargaNo(self.v_caso.get(), nid, d.result["fx"], d.result["fy"], d.result["mz"]))
            self._modificou()

    def dialogo_carga_barra(self, bid, ev=None):
        if not self._exige_caso():
            return
        b = self.modelo.barras[bid]
        L = self.modelo.geometria(b)[0]
        pos = L / 2
        if ev is not None and ev.xdata is not None:
            n1 = self.modelo.nos[b.ni]
            _, c, s = self.modelo.geometria(b)
            pos = max(0.0, min(L, (ev.xdata - n1.x) * c + (ev.ydata - n1.y) * s))
            dx = self.modelo.grade["dx"]
            pos = round(pos / (dx / 2)) * (dx / 2) if self.v_snap.get() else pos
        dirs = {v: k for k, v in core.DIRECOES.items()}
        d = FormDialog(self, f"Carga na barra {bid}", [
            dict(key="tipo", label="Tipo", type="combo", options=["Distribuída (kN/m)", "Pontual (kN)"], default="Distribuída (kN/m)"),
            dict(key="dir", label="Direção", type="combo", options=list(dirs), default="Global Y"),
            dict(key="valor", label="Valor (q em kN/m ou P em kN)", default=-1.0),
            dict(key="pos", label="Posição da pontual a partir do nó i (m)", default=round(pos, 4)),
            dict(key="proj", label="Distribuída por comprimento projetado (GX/GY)", type="check", default=False)],
            texto=f"Caso: {self.v_caso.get()}\nBarra {bid}: nós {b.ni}→{b.nj},  L = {fmt(L, 3)} m.\n"
                  "Sinal pelo eixo: Global Y para baixo → negativo. Local y = x local girado 90° anti-horário.")
        if d.result:
            r = d.result
            tipo = "distribuida" if r["tipo"].startswith("Dist") else "pontual"
            if tipo == "pontual" and not (0 <= r["pos"] <= L + 1e-9):
                messagebox.showerror("Carga", f"Posição fora da barra (0 a {L:.3f} m).")
                return
            self.modelo.cargas_barra.append(core.CargaBarra(self.v_caso.get(), bid, tipo, dirs[r["dir"]],
                                                            r["valor"], r["pos"], r["proj"]))
            self._modificou()

    def excluir_carga(self):
        for iid in self.tv_cargas.selection():
            idx = int(iid[1:])
            (self.modelo.cargas_no if iid[0] == "n" else self.modelo.cargas_barra)[idx] = None
        self.modelo.cargas_no = [c for c in self.modelo.cargas_no if c]
        self.modelo.cargas_barra = [c for c in self.modelo.cargas_barra if c]
        self._modificou()

    def nova_comb(self):
        if not self.modelo.casos:
            return
        d = CombDialog(self, self.modelo)
        if d.result:
            self.modelo.combinacoes[d.result.nome] = d.result
            self._modificou()

    def _comb_sel(self):
        s = self.lb_comb.curselection()
        return list(self.modelo.combinacoes)[s[0]] if s else None

    def editar_comb(self):
        nome = self._comb_sel()
        if nome:
            d = CombDialog(self, self.modelo, self.modelo.combinacoes[nome])
            if d.result:
                self.modelo.combinacoes.pop(nome)
                self.modelo.combinacoes[d.result.nome] = d.result
                self._modificou()

    def excluir_comb_auto(self):
        n = sum(1 for c in self.modelo.combinacoes.values() if c.auto)
        if n and messagebox.askyesno("Combinações", f"Excluir as {n} combinações geradas automaticamente?"):
            for k in [k for k, c in self.modelo.combinacoes.items() if c.auto]:
                self.modelo.combinacoes.pop(k)
            self._modificou()

    def gerar_nbr8800(self):
        if not self.modelo.casos:
            messagebox.showinfo("Combinações", "Crie os casos de carga primeiro.")
            return
        d = GerarCombDialog(self, self.modelo)
        if d.aplicou:
            self._modificou()
            self.status.set(d.aplicou)

    def importar_dxf(self):
        p = filedialog.askopenfilename(title="Importar barras de DXF", filetypes=[("AutoCAD DXF", "*.dxf")])
        if not p:
            return
        try:
            segs = dxf.ler_segmentos(p)
        except Exception as e:
            messagebox.showerror("DXF", f"Não foi possível ler o arquivo:\n{e}")
            return
        if not segs:
            messagebox.showwarning("DXF", "Nenhuma LINE/POLYLINE encontrada no espaço do modelo.")
            return
        d = ImportarDXFDialog(self, self.modelo, segs, p)
        if d.relatorio:
            self._modificou()
            self.zoom_extensao()

    def excluir_comb(self):
        nome = self._comb_sel()
        if nome:
            self.modelo.combinacoes.pop(nome)
            self._modificou()

    # ============================================================ cálculo
    def calcular(self):
        soltos = [nid for nid in self.modelo.nos if self.modelo.barra_sob_no(nid)]
        if soltos and messagebox.askyesno(
                "Nós soltos sobre barras",
                f"Os nós {soltos} estão sobre barras mas não estão ligados a elas "
                "(cargas e apoios nesses nós seriam ignorados).\n\nDividir as barras nesses nós?"):
            self._conectar_soltos()
            self._modificou()
        erros, avisos = self.modelo.validar()
        if erros:
            messagebox.showerror("Modelo incompleto", "\n".join(erros))
            return
        try:
            self.config(cursor="watch")
            self.update()
            self.resultados = core.analisar(self.modelo)
        except core.ErroAnalise as e:
            messagebox.showerror("Erro na análise", str(e))
            self.resultados = None
            return
        except Exception as e:
            messagebox.showerror("Erro inesperado", f"{type(e).__name__}: {e}")
            self.resultados = None
            return
        finally:
            self.config(cursor="")
        nomes = list(self.resultados)
        self.cb_res.configure(values=nomes)
        if self.v_res.get() not in nomes:
            self.v_res.set(nomes[0])
        self.nb.select(3)
        self.v_modo.set("sel")
        self.pendente = None
        if self.v_diag.get() == "Estrutura":
            self.v_diag.set("Momento (M)")
        self.redesenha()
        self._info_resultado()
        self.status.set("Análise concluída. " + (" | ".join(avisos) if avisos else ""))

    def _res_atual(self):
        if not self.resultados:
            return None
        return self.resultados.get(self.v_res.get())

    def _info_resultado(self):
        r = self._res_atual()
        t = self.txt_res
        t.delete("1.0", "end")
        if not r:
            return
        m = self.modelo
        out = [f"[{self.v_res.get()}]"]
        itens = self.sel or []
        for tipo, i in itens[:4]:
            if tipo == "barra" and i in m.barras:
                xs, N, V, M = core.esforcos_barra(m, r, i)
                (xM, Mx), (xm, Mn) = core.extremos(xs, M)
                (xV, Vx), (xv, Vn) = core.extremos(xs, V)
                xf, fl = core.deslocamento_max_barra(m, r, i)
                xa, dxa, dya = core.desloc_max_absoluto(m, r, i)
                L = xs[-1]
                out += [f"Barra {i} ({m.barras[i].secao}), L={L:.3f} m",
                        f"          nó i      nó j",
                        f"  N  {N[0]:>9.2f} {N[-1]:>9.2f}  kN",
                        f"  V  {V[0]:>9.2f} {V[-1]:>9.2f}  kN",
                        f"  M  {M[0]:>9.2f} {M[-1]:>9.2f}  kN·m",
                        f"  Mmax={Mx:.2f} @x={xM:.3f}  Mmin={Mn:.2f} @x={xm:.3f}",
                        f"  Vmax={Vx:.2f}  Vmin={Vn:.2f}",
                        f"  Flecha rel. à corda: {fl:.2f} mm @x={xf:.2f} (L/{abs(L * 1e3 / fl):.0f})" if abs(fl) > 1e-9 else "  Flecha ≈ 0",
                        f"  δ máx. absoluto: {math.hypot(dxa, dya) * 1e3:.2f} mm @x={xa:.2f} "
                        f"(δx={dxa * 1e3:.2f}, δy={dya * 1e3:.2f})",
                        ""]
            elif tipo == "no" and i in r.desloc:
                d = r.desloc[i]
                out.append(f"Nó {i}: ux={d[0]*1e3:.3f} mm  uy={d[1]*1e3:.3f} mm  rz={d[2]:.3e} rad")
                if i in r.reacoes:
                    R = r.reacoes[i]
                    out.append(f"   Rx={R[0]:.2f} kN  Ry={R[1]:.2f} kN  Mz={R[2]:.2f} kN·m")
                out.append("")
        out.append("REAÇÕES DE APOIO")
        out.append(f"{'Nó':>4}{'Rx kN':>10}{'Ry kN':>10}{'Mz kN·m':>10}")
        sx = sy = 0
        for nid, R in sorted(r.reacoes.items()):
            out.append(f"{nid:>4}{R[0]:>10.2f}{R[1]:>10.2f}{R[2]:>10.2f}")
            sx += R[0]
            sy += R[1]
        out.append(f"{'Σ':>4}{sx:>10.2f}{sy:>10.2f}")
        cx, cy = r.soma_cargas
        out.append(f"{'Σcarg':>6}{cx:>8.2f}{cy:>10.2f}   (cargas aplicadas)")
        erro = math.hypot(sx + cx, sy + cy)
        if erro > 1e-3 * max(1.0, math.hypot(cx, cy)):
            out.append(f"⚠ EQUILÍBRIO NÃO FECHA: diferença = {erro:.3f} kN")
        else:
            out.append("   ✔ equilíbrio global: ΣR + Σcargas = 0")
        dmax = max(r.desloc.items(), key=lambda kv: math.hypot(kv[1][0], kv[1][1]))
        out.append(f"\nMaior desloc. nodal: nó {dmax[0]} = {math.hypot(*dmax[1][:2]) * 1e3:.2f} mm")
        melhor = None
        for bid in m.barras:
            xm, dxm, dym = core.desloc_max_absoluto(m, r, bid)
            if melhor is None or math.hypot(dxm, dym) > melhor[0]:
                melhor = (math.hypot(dxm, dym), bid, xm, dxm, dym)
        if melhor:
            out.append(f"Maior desloc. nas barras: barra {melhor[1]} @x={melhor[2]:.2f} m = {melhor[0] * 1e3:.2f} mm")
            out.append(f"   (δx = {melhor[3] * 1e3:.2f} mm, δy = {melhor[4] * 1e3:.2f} mm)")
        t.insert("1.0", "\n".join(out))

    # ============================================================ desenho
    def redesenha(self):
        ax = self.ax
        xl, yl = ax.get_xlim(), ax.get_ylim()
        ax.clear()
        xl, yl = self._limites_iguais(xl, yl)
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        self._grade()
        self._linha_prev = None
        self._hover = {}          # bid -> (xs, pontos do eixo, pontos da curva)
        self._hover_meta = None
        self._hover_art = []
        sz = 0.025 * self._span()
        r = self._res_atual()
        diag = self.v_diag.get() if r else "Estrutura"
        m = self.modelo
        cinza = diag not in ("Estrutura",)
        # barras
        for b in m.barras.values():
            if not self.barra_visivel(b):
                continue
            n1, n2 = m.nos[b.ni], m.nos[b.nj]
            s = m.secoes.get(b.secao)
            cor = "#9a9a9a" if cinza else (s.cor if s else "k")
            vv = self.verif.get(b.id) if (not cinza and self.v_cor_verif.get()) else None
            if vv:
                cor = {"OK": "#1b8a2f", "FALHA": "#8b0000"}.get(vv.status, "#9a9a9a")
            selec = ("barra", b.id) in self.sel
            L, c, sn = m.geometria(b)
            if self.v_extr.get() and not cinza and s:
                self._desenha_extrudado(b, s, n1, n2, L, c, sn, selec, cor)
            else:
                ax.plot([n1.x, n2.x], [n1.y, n2.y], color="#ffbf00" if selec else cor,
                        lw=5 if selec else (1.6 if cinza else 2.6), zorder=3, solid_capstyle="round")
                if selec:
                    ax.plot([n1.x, n2.x], [n1.y, n2.y], color=cor, lw=2.2, zorder=3)
            if self.v_eixos.get():
                self._desenha_eixos_locais(b, n1, L, c, sn, sz)
            cc = m.casos.get(self.v_caso.get())
            if not b.pp and not cinza and self.v_cargas.get() and cc and cc.peso_proprio:
                ax.text((n1.x + n2.x) / 2 + sn * 0.9 * sz, (n1.y + n2.y) / 2 - c * 0.9 * sz, "sem PP",
                        fontsize=7, color="#888", style="italic", ha="center", va="center", zorder=7)
            for rot, n, sg in ((b.rot_i, n1, 1), (b.rot_j, n2, -1)):
                if rot:
                    off = min(0.7 * sz, 0.2 * L)
                    ax.plot(n.x + sg * c * off, n.y + sg * sn * off, "o", ms=5.5, mfc="white",
                            mec="#9a9a9a" if cinza else cor, mew=1.3, zorder=6)
            if vv and vv.status != "NA":
                ax.text((n1.x + n2.x) / 2 + sn * 0.9 * sz, (n1.y + n2.y) / 2 - c * 0.9 * sz,
                        f"{vv.eta_max * 100:.0f}%", fontsize=7, color=cor, fontweight="bold",
                        ha="center", va="center", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))
            if self.v_num_bar.get():
                ax.text((n1.x + n2.x) / 2 - sn * 0.6 * sz, (n1.y + n2.y) / 2 + c * 0.6 * sz, f"{b.id}",
                        fontsize=8, color="#555", ha="center", va="center", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.15", fc="#f5f5f5", ec="#bbb", lw=0.5))
        # nós
        visiveis = {b.ni for b in m.barras.values() if self.barra_visivel(b)} | \
                   {b.nj for b in m.barras.values() if self.barra_visivel(b)} | \
                   {n for n in m.nos if not any(n in (b.ni, b.nj) for b in m.barras.values())}
        for n in m.nos.values():
            if n.id not in visiveis:
                continue
            selec = ("no", n.id) in self.sel or n.id == self.pendente
            ax.plot(n.x, n.y, "o", ms=8 if selec else 4.5, color="#ff7f00" if selec else "k", zorder=8)
            if self.v_num_nos.get():
                ax.text(n.x + 0.35 * sz, n.y + 0.35 * sz, str(n.id), fontsize=8, color="#b00", zorder=9)
        # apoios
        for a in m.apoios.values():
            if a.no in m.nos and a.no in visiveis:
                self._desenha_apoio(a, sz)
        # cargas (modo modelo)
        if diag == "Estrutura" and self.v_cargas.get():
            self._desenha_cargas(sz)
        if r and diag != "Estrutura":
            if diag == "Deformada":
                self._desenha_deformada(r)
            elif diag == "Reações":
                self._desenha_reacoes(r, sz)
            else:
                self._desenha_diagrama(r, diag, sz)
        if self.pendente and self.v_modo.get() == "barra":
            p = m.nos[self.pendente]
            c = self.cursor or (p.x, p.y)
            self._linha_prev, = ax.plot([p.x, c[0]], [p.y, c[1]], "--", color="#ff7f00", lw=1.2, zorder=10)
        titulo = f"{m.titulo}"
        if r and diag != "Estrutura":
            titulo += f"   —   {diag}   —   {self.v_res.get()}"
        elif self.v_cargas.get() and self.v_caso.get():
            titulo += f"   —   cargas: {self.v_caso.get()}"
        ax.set_title(titulo, fontsize=10, loc="left")
        self.canvas.draw_idle()

    def _limites_iguais(self, xl, yl):
        """Escala 1:1 — expande o eixo necessário conforme a proporção em pixels do gráfico."""
        bb = self.ax.get_window_extent()
        if bb.width < 10 or bb.height < 10:
            return xl, yl
        w, h = xl[1] - xl[0], yl[1] - yl[0]
        alvo = bb.width / bb.height
        cx, cy = sum(xl) / 2, sum(yl) / 2
        if w / h < alvo:
            w = h * alvo
        else:
            h = w / alvo
        return (cx - w / 2, cx + w / 2), (cy - h / 2, cy + h / 2)

    def _desenha_extrudado(self, b, s, n1, n2, L, c, sn, selec, cor=None):
        """Vista em elevação da barra com a altura real da seção (girada de β)."""
        ax = self.ax
        ymin, ymax, ints = s.linhas_elevacao(b.beta)
        nx, ny = -sn, c
        def faixa(y):
            y /= 1000.0
            return ([n1.x + nx * y, n2.x + nx * y], [n1.y + ny * y, n2.y + ny * y])
        a0, a1 = ymin / 1000.0, ymax / 1000.0
        pts = [(n1.x + nx * a0, n1.y + ny * a0), (n2.x + nx * a0, n2.y + ny * a0),
               (n2.x + nx * a1, n2.y + ny * a1), (n1.x + nx * a1, n1.y + ny * a1)]
        cor = cor or s.cor
        ax.add_patch(Polygon(pts, closed=True, fc=cor, alpha=0.30, ec="none", zorder=2))
        ec = "#ff9800" if selec else cor
        ax.add_patch(Polygon(pts, closed=True, fill=False, ec=ec, lw=2.2 if selec else 1.1, zorder=3))
        for y in ints:
            ax.plot(*faixa(y), color=cor, lw=0.6, zorder=3)
        ax.plot([n1.x, n2.x], [n1.y, n2.y], color="#555", lw=0.6, ls=(0, (8, 3, 2, 3)), zorder=3)

    def _desenha_eixos_locais(self, b, n1, L, c, sn, sz):
        ax = self.ax
        la = min(2.0 * sz, 0.3 * L)
        ox, oy = n1.x + c * 0.3 * L, n1.y + sn * 0.3 * L
        ax.annotate("", xy=(ox + c * la, oy + sn * la), xytext=(ox, oy),
                    arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=1.3, mutation_scale=9), zorder=13)
        ax.annotate("", xy=(ox - sn * la, oy + c * la), xytext=(ox, oy),
                    arrowprops=dict(arrowstyle="-|>", color="#2e7d32", lw=1.3, mutation_scale=9), zorder=13)
        ax.text(ox + c * la * 1.25, oy + sn * la * 1.25, "x", color="#d62728", fontsize=8,
                fontweight="bold", ha="center", va="center", zorder=13)
        ax.text(ox - sn * la * 1.25, oy + c * la * 1.25, "y", color="#2e7d32", fontsize=8,
                fontweight="bold", ha="center", va="center", zorder=13)
        if b.beta:
            ax.text(ox + (c - sn) * la * 0.6, oy + (sn + c) * la * 0.6, f"β={b.beta:g}°", color="#6a1b9a",
                    fontsize=7, ha="left", va="bottom", zorder=13)

    def _grade(self):
        ax = self.ax
        g = self.modelo.grade
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        fx = max(1, math.ceil((x1 - x0) / g["dx"] / 60))
        fy = max(1, math.ceil((y1 - y0) / g["dy"] / 60))
        ax.xaxis.set_major_locator(MultipleLocator(g["dx"] * fx))
        ax.yaxis.set_major_locator(MultipleLocator(g["dy"] * fy))
        ax.grid(True, color="#dddddd", lw=0.6, zorder=0)
        ax.tick_params(labelsize=7, colors="#777")
        for sp in ax.spines.values():
            sp.set_color("#cccccc")

    def _desenha_apoio(self, a, sz):
        ax = self.ax
        n = self.modelo.nos[a.no]
        x, y = n.x, n.y
        t = (a.ux, a.uy, a.rz)
        cor = "#2d6a2d"
        if t == ("fixo", "fixo", "fixo"):
            ax.add_patch(Rectangle((x - 1.2 * sz, y - 0.5 * sz), 2.4 * sz, 0.5 * sz, fc="none", ec=cor, hatch="////", lw=1, zorder=4))
            ax.plot([x - 1.2 * sz, x + 1.2 * sz], [y, y], color=cor, lw=2.2, zorder=4)
            return
        if a.uy in ("fixo",):
            tri = [(x, y), (x - 0.8 * sz, y - 1.2 * sz), (x + 0.8 * sz, y - 1.2 * sz)]
            ax.add_patch(Polygon(tri, closed=True, fc="#dff0df", ec=cor, lw=1.3, zorder=4))
            base = y - 1.2 * sz
            if a.ux == "fixo":
                ax.add_patch(Rectangle((x - 1.1 * sz, base - 0.35 * sz), 2.2 * sz, 0.35 * sz, fc="none", ec=cor, hatch="////", lw=0.8, zorder=4))
            else:
                ax.plot([x - 1.1 * sz, x + 1.1 * sz], [base - 0.3 * sz] * 2, color=cor, lw=1.3, zorder=4)
        elif a.ux == "fixo":
            tri = [(x, y), (x - 1.2 * sz, y - 0.8 * sz), (x - 1.2 * sz, y + 0.8 * sz)]
            ax.add_patch(Polygon(tri, closed=True, fc="#dff0df", ec=cor, lw=1.3, zorder=4))
            ax.plot([x - 1.5 * sz] * 2, [y - 1.1 * sz, y + 1.1 * sz], color=cor, lw=1.3, zorder=4)
        if a.rz == "fixo":
            ax.add_patch(Rectangle((x - 0.35 * sz, y - 0.35 * sz), 0.7 * sz, 0.7 * sz, fc=cor, ec=cor, zorder=5))
        # molas
        def zig(x0, y0, dx, dy):
            k = 6
            px, py = -dy, dx
            pts = [(x0, y0)]
            for i in range(1, k + 1):
                f = i / (k + 1)
                s = 0.35 * sz * (1 if i % 2 else -1)
                pts.append((x0 + dx * f + px * s, y0 + dy * f + py * s))
            pts.append((x0 + dx, y0 + dy))
            ax.plot(*zip(*pts), color="#8b4513", lw=1.2, zorder=4)
            ax.plot([x0 + dx - 0.5 * sz * abs(dy) / (abs(dx) + abs(dy)), x0 + dx + 0.5 * sz * abs(dy) / (abs(dx) + abs(dy))],
                    [y0 + dy - 0.5 * sz * abs(dx) / (abs(dx) + abs(dy)), y0 + dy + 0.5 * sz * abs(dx) / (abs(dx) + abs(dy))],
                    color="#8b4513", lw=2, zorder=4)
        if a.uy == "mola":
            zig(x, y, 0, -2.2 * sz)
            ax.text(x + 0.5 * sz, y - 1.6 * sz, f"ky={a.ky:g}", fontsize=7, color="#8b4513")
        if a.ux == "mola":
            zig(x, y, -2.2 * sz, 0)
            ax.text(x - 2.4 * sz, y + 0.4 * sz, f"kx={a.kx:g}", fontsize=7, color="#8b4513")
        if a.rz == "mola":
            ax.annotate("", xy=(x + 0.9 * sz, y), xytext=(x, y + 0.9 * sz),
                        arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.6", color="#8b4513", lw=1.5))
            ax.text(x + 0.9 * sz, y + 0.6 * sz, f"kθ={a.kr:g}", fontsize=7, color="#8b4513")

    def _seta(self, x, y, ux, uy, comp, cor, txt=None, lw=1.4):
        """Seta com a ponta em (x, y), apontando na direção (ux, uy)."""
        ax = self.ax
        ax.annotate("", xy=(x, y), xytext=(x - ux * comp, y - uy * comp),
                    arrowprops=dict(arrowstyle="-|>", color=cor, lw=lw, mutation_scale=10), zorder=11)
        if txt:
            ax.text(x - ux * comp * 1.1, y - uy * comp * 1.1, txt, fontsize=8, color=cor,
                    ha="center", va="center", zorder=12,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))

    def _desenha_cargas(self, sz):
        m = self.modelo
        caso = self.v_caso.get()
        cor = "#c00000"
        comp = 2.5 * sz
        for cg in m.cargas_no:
            if cg.caso != caso or cg.no not in m.nos:
                continue
            n = m.nos[cg.no]
            if cg.fx:
                self._seta(n.x, n.y, math.copysign(1, cg.fx), 0, comp, cor, f"{fmt(abs(cg.fx))} kN")
            if cg.fy:
                self._seta(n.x, n.y, 0, math.copysign(1, cg.fy), comp, cor, f"{fmt(abs(cg.fy))} kN")
            if cg.mz:
                r = 1.2 * sz
                ccw = cg.mz > 0
                self.ax.annotate("", xy=(n.x - r, n.y) if ccw else (n.x + r, n.y), xytext=(n.x + r, n.y) if ccw else (n.x - r, n.y),
                                 arrowprops=dict(arrowstyle="-|>", color=cor, lw=1.4,
                                                 connectionstyle="arc3,rad=" + ("-0.9" if ccw else "0.9")), zorder=11)
                self.ax.text(n.x, n.y + 1.7 * r, f"{fmt(abs(cg.mz))} kN·m", fontsize=8, color=cor, ha="center")
        for cg in m.cargas_barra:
            if cg.caso != caso or cg.barra not in m.barras:
                continue
            b = m.barras[cg.barra]
            if not self.barra_visivel(b):
                continue
            n1 = m.nos[b.ni]
            L, c, s = m.geometria(b)
            if cg.direcao == "GX":
                ux, uy = 1.0, 0.0
            elif cg.direcao == "GY":
                ux, uy = 0.0, 1.0
            elif cg.direcao == "LX":
                ux, uy = c, s
            else:
                ux, uy = -s, c
            sg = 1 if cg.valor >= 0 else -1
            ux, uy = ux * sg, uy * sg
            if cg.tipo == "pontual":
                px, py = n1.x + c * cg.pos, n1.y + s * cg.pos
                self._seta(px, py, ux, uy, comp, cor, f"{fmt(abs(cg.valor))} kN")
            else:
                k = max(3, min(12, int(L / (2 * sz)) + 1))
                cc = "#b35900"
                if cg.direcao in ("LX",):
                    for i in range(k):
                        f = (i + 0.5) / k
                        self._seta(n1.x + c * L * f, n1.y + s * L * f, ux, uy, 0.8 * comp / k * 3, cc, lw=1.0)
                    self.ax.text(n1.x + c * L / 2 - s * sz, n1.y + s * L / 2 + c * sz,
                                 f"{fmt(abs(cg.valor))} kN/m", fontsize=8, color=cc, ha="center")
                    continue
                h = 0.7 * comp
                tails = []
                for i in range(k + 1):
                    f = i / k
                    px, py = n1.x + c * L * f, n1.y + s * L * f
                    self._seta(px, py, ux, uy, h, cc, lw=1.0)
                    tails.append((px - ux * h, py - uy * h))
                self.ax.plot(*zip(*tails), color=cc, lw=1.0, zorder=11)
                mx, my = tails[len(tails) // 2]
                self.ax.text(mx - ux * 0.5 * sz, my - uy * 0.5 * sz,
                             f"{fmt(abs(cg.valor))} kN/m" + (" (proj.)" if cg.projetada else ""),
                             fontsize=8, color=cc, ha="center", va="center", zorder=12,
                             bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
        # peso próprio
        cc = m.casos.get(caso)
        if cc and cc.peso_proprio and m.barras:
            self.ax.text(0.01, 0.01, "+ peso próprio automático neste caso", transform=self.ax.transAxes,
                         fontsize=8, color="#b35900")

    def _desenha_diagrama(self, r, diag, sz):
        ax, m = self.ax, self.modelo
        idx = {"Normal (N)": 1, "Cortante (V)": 2, "Momento (M)": 3}[diag]
        un = "kN·m" if idx == 3 else "kN"
        dados = {}
        vmax = 1e-12
        for b in m.barras.values():
            if not self.barra_visivel(b):
                continue
            d = core.esforcos_barra(m, r, b.id, 40)
            dados[b.id] = (d[0], d[idx])
            vmax = max(vmax, max(abs(v) for v in d[idx]))
        esc = self.v_escala.get() * 0.10 * self._span() / vmax
        sinal = -1 if idx == 3 else 1          # M desenhado do lado tracionado (−y local p/ M>0)
        self._hover_meta = {"tipo": "esforco", "idx": idx, "esc": esc, "sinal": sinal, "un": un,
                            "nome": {1: "N", 2: "V", 3: "M"}[idx]}
        cpos, cneg = {1: ("#1f5fbf", "#c62828"), 2: ("#2e7d32", "#ef6c00"), 3: ("#6a1b9a", "#6a1b9a")}[idx]
        for bid, (xs, vals) in dados.items():
            b = m.barras[bid]
            n1 = m.nos[b.ni]
            L, c, s = m.geometria(b)
            nx, ny = -s, c
            pts_base = [(n1.x + c * x, n1.y + s * x) for x in xs]
            pts = [(bx + nx * v * esc * sinal, by + ny * v * esc * sinal) for (bx, by), v in zip(pts_base, vals)]
            for parte, cor in ((lambda v: max(v, 0.0), cpos), (lambda v: min(v, 0.0), cneg)):
                pp = [(bx + nx * parte(v) * esc * sinal, by + ny * parte(v) * esc * sinal) for (bx, by), v in zip(pts_base, vals)]
                if any(abs(parte(v)) > 1e-9 for v in vals):
                    ax.add_patch(Polygon(pts_base + pp[::-1], closed=True, fc=cor, alpha=0.25, ec="none", zorder=2))
            ax.plot(*zip(*pts), color=cpos if idx == 3 else "#333", lw=1.1, zorder=4)
            self._hover[bid] = (xs, pts_base, pts)
            if self.v_valores.get():
                (xM, vM), (xm, vm) = core.extremos(xs, vals)
                marcar = {0: vals[0], len(xs) - 1: vals[-1]}
                for xe, ve in ((xM, vM), (xm, vm)):
                    k = xs.index(xe)
                    if 0 < k < len(xs) - 1 and abs(ve - vals[0]) > 1e-3 * vmax and abs(ve - vals[-1]) > 1e-3 * vmax:
                        marcar[k] = ve
                if idx == 1 and abs(vals[0] - vals[-1]) < 1e-6 * vmax:
                    marcar = {len(xs) // 2: vals[len(xs) // 2]}
                for k, v in marcar.items():
                    if abs(v) < 1e-3 * vmax and len(marcar) > 1:
                        continue
                    bx, by = pts_base[k]
                    off = v * esc * sinal
                    tx, ty = bx + nx * (off + math.copysign(0.6 * sz, off if off else 1)), by + ny * (off + math.copysign(0.6 * sz, off if off else 1))
                    ax.text(tx, ty, fmt(v), fontsize=7.5, ha="center", va="center", zorder=12,
                            color=cpos if v >= 0 else cneg,
                            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75))
        leg = {1: "N: azul = tração (+), vermelho = compressão (−)",
               2: "V: verde = positivo, laranja = negativo",
               3: "M: desenhado do lado tracionado"}[idx]
        ax.text(0.01, 0.01, f"{leg}   |   máx |{diag[0]}| = {fmt(vmax)} {un}", transform=ax.transAxes, fontsize=8, color="#333")

    def _desenha_deformada(self, r):
        ax, m = self.ax, self.modelo
        dmax = 1e-12
        curvas = {}
        for b in m.barras.values():
            if not self.barra_visivel(b):
                continue
            L = m.geometria(b)[0]
            xs, DX, DY = core.desloc_global_barra(m, r, b.id, [L * k / 40 for k in range(41)])
            curvas[b.id] = (xs, DX, DY)
            dmax = max(dmax, max(math.hypot(u, v) for u, v in zip(DX, DY)))
        esc = self.v_escala.get() * 0.08 * self._span() / dmax
        self._hover_meta = {"tipo": "desloc", "esc": esc}
        sz = 0.025 * self._span()
        for bid, (xs, DX, DY) in curvas.items():
            b = m.barras[bid]
            n1 = m.nos[b.ni]
            L, c, s = m.geometria(b)
            base = [(n1.x + c * x, n1.y + s * x) for x in xs]
            curva = [(bx + dx * esc, by + dy * esc) for (bx, by), dx, dy in zip(base, DX, DY)]
            ax.plot(*zip(*curva), color="#d62728", lw=1.8, zorder=5)
            self._hover[bid] = (xs, base, curva)
            if self.v_valores.get():
                xm, dxm, dym = core.desloc_max_absoluto(m, r, bid)
                dm = math.hypot(dxm, dym)
                d_ext = max(math.hypot(DX[0], DY[0]), math.hypot(DX[-1], DY[-1]))
                if dm > 1e-2 * dmax and dm > 1.001 * d_ext and 1e-6 * L < xm < L * (1 - 1e-6):
                    px, py = n1.x + c * xm + dxm * esc, n1.y + s * xm + dym * esc
                    ax.plot(px, py, "o", ms=5, color="#d62728", zorder=12)
                    ax.plot([n1.x + c * xm, px], [n1.y + s * xm, py], ":", color="#d62728", lw=0.8, zorder=4)
                    ux, uy = (dxm / dm, dym / dm)
                    ax.text(px + ux * 0.8 * sz, py + uy * 0.8 * sz,
                            f"δmáx = {fmt(dm * 1e3)} mm\n(x = {fmt(xm)} m)", fontsize=7.5, color="#b71c1c",
                            ha="center", va="center", zorder=13,
                            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#e57373", lw=0.6, alpha=0.9))
        if self.v_valores.get():
            for nid, d in r.desloc.items():
                n = m.nos[nid]
                dd = math.hypot(d[0], d[1]) * 1e3
                if dd > 1e-3 * dmax * 1e3:
                    ax.text(n.x + d[0] * esc, n.y + d[1] * esc, f" {fmt(dd)} mm", fontsize=7, color="#d62728", zorder=12)
        ax.text(0.01, 0.01, f"Deformada ampliada {esc:.0f}×   |   máx = {fmt(dmax * 1e3)} mm   |   "
                            "passe o mouse sobre a barra para ler o deslocamento",
                transform=ax.transAxes, fontsize=8, color="#333")

    def _desenha_reacoes(self, r, sz):
        m = self.modelo
        cor = "#0b5394"
        for nid, R in r.reacoes.items():
            n = m.nos[nid]
            comp = 3 * sz
            if abs(R[0]) > 1e-6:
                self._seta(n.x, n.y, math.copysign(1, R[0]), 0, comp, cor,
                           f"{fmt(abs(R[0]))} kN")
            if abs(R[1]) > 1e-6:
                self._seta(n.x, n.y - 1.4 * sz, 0, math.copysign(1, R[1]), comp, cor, f"{fmt(abs(R[1]))} kN")
            if abs(R[2]) > 1e-6:
                self.ax.text(n.x + 1.2 * sz, n.y - 2.2 * sz, f"M = {fmt(R[2])} kN·m {'↺' if R[2] > 0 else '↻'}",
                             fontsize=8, color=cor, zorder=12)
        self.ax.text(0.01, 0.01, "Reações nos eixos globais (setas no sentido real)", transform=self.ax.transAxes,
                     fontsize=8, color="#333")

    def zoom_extensao(self):
        m = self.modelo
        if not m.nos:
            return
        xs = [n.x for n in m.nos.values()]
        ys = [n.y for n in m.nos.values()]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        mg = 0.18 * max(w, h, 1.0)
        self.ax.set_xlim(min(xs) - mg, max(xs) + mg)
        self.ax.set_ylim(min(ys) - mg, max(ys) + mg)
        self.redesenha()

    # ============================================================ arquivos
    def novo(self):
        if messagebox.askyesno("Novo", "Descartar o modelo atual?"):
            self.modelo = self._modelo_inicial()
            self._carrega_cfg_proj()
            self.arquivo = None
            self.sel = []
            self.visivel = {}
            self._modificou()

    def abrir(self):
        p = filedialog.askopenfilename(filetypes=[("Modelo Pórtico 2D", "*.p2d.json *.json")])
        if not p:
            return
        try:
            self.modelo = core.Modelo.abrir(p)
        except Exception as e:
            messagebox.showerror("Abrir", str(e))
            return
        self.arquivo = p
        self.sel, self.visivel = [], {}
        self.v_secao.set(next(iter(self.modelo.secoes), ""))
        self.v_caso.set(next(iter(self.modelo.casos), ""))
        self._carrega_cfg_proj()
        self._modificou()
        self.zoom_extensao()

    def salvar(self, como=False):
        if como or not self.arquivo:
            p = filedialog.asksaveasfilename(defaultextension=".p2d.json",
                                             filetypes=[("Modelo Pórtico 2D", "*.p2d.json")])
            if not p:
                return
            self.arquivo = p
        self.modelo.salvar(self.arquivo)
        self.status.set(f"Salvo em {self.arquivo}")

    def ver_memorial(self):
        txt = core.memorial_txt(self.modelo, self.resultados)
        w = tk.Toplevel(self)
        w.title("Memorial")
        w.geometry("980x700")
        t = tk.Text(w, font=("Consolas", 9), wrap="none")
        sb = ttk.Scrollbar(w, command=t.yview)
        t.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        t.pack(fill="both", expand=True)
        t.insert("1.0", txt)

    def exportar_memorial(self):
        p = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Texto", "*.txt")])
        if p:
            with open(p, "w", encoding="utf-8") as f:
                f.write(core.memorial_txt(self.modelo, self.resultados))
            self.status.set(f"Memorial exportado: {p}")

    def exportar_figura(self):
        p = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png"), ("PDF", "*.pdf")])
        if p:
            self.fig.savefig(p, dpi=200)

    def ajuda(self):
        messagebox.showinfo("Convenções", (
            "UNIDADES: kN, m, kN·m (seções em mm).\n\n"
            "SINAIS\n• Cargas nodais: eixos globais (X →, Y ↑, Mz ↺ +). Gravidade = Fy negativo.\n"
            "• Barra: x local do nó i → j; y local = x girado 90° anti-horário.\n"
            "• N > 0 tração. M desenhado do lado tracionado.\n\n"
            "LIGAÇÕES: 'Rótula–Rótula' = barra de treliça (M = 0 nas pontas; carga no vão gera M no vão).\n\n"
            "ATALHOS: Esc cancela | Del exclui seleção | F5 calcula | F zoom extensão | roda do mouse = zoom | "
            "botão direito cancela a barra em lançamento.\n\n"
            "Análise elástica linear de 1ª ordem (OpenSees elasticBeamColumn). Combinações por superposição."))


if __name__ == "__main__":
    App().mainloop()
