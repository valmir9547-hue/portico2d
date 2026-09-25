# -*- coding: utf-8 -*-
"""
verificacao_nbr8800.py — Dimensionamento de barras de aço (perfis I/H laminados W e
soldados I, dois eixos de simetria) pela ABNT NBR 8800 — SEM interface gráfica.

Formulação implementada (artigos da NBR 8800:2008, mesma numeração citada no memorial
de referência; CONFIRA eventuais alterações da revisão 2024):
  • 5.2   Tração — escoamento da seção bruta e ruptura da seção líquida efetiva
  • 5.3   Compressão — χ (5.3.3), Q (Anexo F: Qs mesas grupos 4/5, Qa alma grupo 2),
          Ne (Anexo E: flexão x, flexão y, torção — seção duplamente simétrica)
  • 5.3.4 / 5.2.8.1  Limitação do índice de esbeltez (200 comprimidas / 300 tracionadas)
  • 5.4.2 + Anexo G  Momento fletor — FLT, FLM, FLA (eixo x) e FLM, FLA (eixo y);
          limite 1,50·W·fy/γa1; Cb pela 5.4.2.3 (Rm = 1,0)
  • 5.4.3 Força cortante — alma (eixo de maior inércia, kv = 5,0) e mesas (menor, kv = 1,2)
  • 5.5.1.2 Interação força axial + momento fletor
  • Anexo C  Deslocamento máximo (flecha) em combinações de serviço
Modelo 2D: esforços fora do plano (My ou Mx, V perpendicular, torção) são nulos.

Unidades internas: kN e cm (tensões em kN/cm²).  Relatório: kN, kN·m, cm, MPa.
O memorial é montado como uma lista de blocos estruturados (classe Memorial) e
renderizado em texto (conferência rápida) ou HTML (impressão / auditoria).
"""
from __future__ import annotations

import html
import math

import portico2d_core as core

CONFIG_PADRAO = {
    "ga1": 1.10, "ga2": 1.35,
    "h_alma": "nbr",       # "nbr": h = d − 2tf − 2r (laminados) | "simplificado": h = d − 2tf
    "sigma_qa": "chi",     # "chi": σ = χ(Q=1)·fy (F.3.2) | "fy": σ = fy (conservador)
    "comb_elu": "ELU",     # "ELU" (combinações geradas tipo ELU) | "todas"
    "comb_els": "ELS",     # "ELS" | "ELS-R" | "ELS-F" | "ELS-QP" | "todas"
}
PROJ_PADRAO = {
    "kp": 1.0,        # K flambagem por flexão NO PLANO   (× L da barra)
    "kf": 1.0,        # K flambagem por flexão FORA DO PLANO (× L)
    "kz": 1.0,        # K flambagem por torção (× L)
    "klb": 1.0,       # Lb / L — distância entre contenções laterais (FLT); 0 = travada continuamente
    # comprimentos informados manualmente (m) — quando preenchidos, substituem K·L
    "Lp": None,       # KL no plano
    "Lf": None,       # KL fora do plano
    "Lz": None,       # KzLz
    "Lb": None,       # Lb (FLT)
    "cb": "auto",     # "auto" (5.4.2.3) ou valor
    "flecha": 250.0,  # limite L / valor (Anexo C, Tabela C.1)
    "flecha_modo": "corda",   # "corda": relativa à corda da barra | "absoluta": deslocamento total ⊥ barra
    "Lflecha": None,  # vão de referência para o limite (m) — None = L da barra (ou 2L no balanço)
    "balanco": False, # barra em balanço: flecha com L = 2·comprimento, Cb = 1,0
    "ct": 1.0,        # coeficiente de redução da área líquida (5.2.5)
    "an": 1.0,        # An / Ag (furos)
}
FAMILIAS_SUPORTADAS = {"W", "I"}


def cfg(modelo):
    c = dict(CONFIG_PADRAO)
    c.update(getattr(modelo, "projeto", {}) or {})
    return c


def pproj(barra):
    p = dict(PROJ_PADRAO)
    p.update(barra.proj or {})
    for k in ("Lp", "Lf", "Lz", "Lb", "Lflecha"):
        if p.get(k) in ("", 0, 0.0) and k != "Lb":
            p[k] = None
    return p


def n(v, nd=2):
    """Número com vírgula decimal e ponto de milhar."""
    if v is None:
        return "—"
    if isinstance(v, float) and math.isinf(v):
        return "∞"
    if abs(v) < 0.5 * 10 ** -nd:
        v = 0.0
    return f"{v:,.{nd}f}".replace(",", "X").replace(".", ",").replace("X", ".")


# =============================================================================
# MEMORIAL ESTRUTURADO
# =============================================================================
class Memorial:
    """Lista de blocos: ("h1", título, artigo) ("h2", título) ("p", texto) ("eq", ...)
    ("cond", texto, ok) ("ver", ...) ("tab", cabeçalho, linhas, alinhamentos) ("aviso", texto)."""

    def __init__(self):
        self.b = []

    def h1(self, titulo, artigo=""):
        self.b.append(("h1", titulo, artigo))

    def h2(self, titulo):
        self.b.append(("h2", titulo))

    def p(self, texto):
        self.b.append(("p", texto))

    def eq(self, simb, expr, subst=None, valor=None, un="", nota=""):
        """simb = expr = subst = valor un   (ex.: Nt,Rd = Ag·fy/γa1 = 19,40×34,50/1,10 = 608,45 kN)"""
        self.b.append(("eq", simb, expr, subst, valor, un, nota))

    def cond(self, texto, ok=None):
        self.b.append(("cond", texto, ok))

    def ver(self, expr, subst, eta, ok=None):
        self.b.append(("ver", expr, subst, eta, eta <= 1.0 + 1e-9 if ok is None else ok))

    def tab(self, cab, linhas, alinh=None):
        self.b.append(("tab", cab, linhas, alinh or ["l"] * len(cab)))

    def aviso(self, texto):
        self.b.append(("aviso", texto))

    def extend(self, outro):
        self.b.extend(outro.b if isinstance(outro, Memorial) else outro)

    # ------------------------------------------------------------------ texto
    def texto(self):
        L = []
        for bl in self.b:
            t = bl[0]
            if t == "h1":
                L += ["", "─" * 96, f"{bl[1]}" + (f"   [{bl[2]}]" if bl[2] else ""), "─" * 96]
            elif t == "h2":
                L += ["", f"  ▸ {bl[1]}"]
            elif t == "p":
                L.append(f"    {bl[1]}")
            elif t == "eq":
                _, s, e, sub, v, un, nota = bl
                pad = " " * (6 + len(s))
                linha = f"      {s} = {e}" if e else f"      {s}"
                if sub is None and v is not None:
                    linha += f" = {v} {un}".rstrip()
                    L.append(linha)
                else:
                    L.append(linha)
                    if sub:
                        L.append(f"{pad} = {sub}")
                    if v is not None:
                        L.append(f"{pad} = {v} {un}".rstrip())
                if nota:
                    L.append(f"{pad}   ({nota})")
            elif t == "cond":
                marca = "" if bl[2] is None else ("  ✔" if bl[2] else "  ✘")
                L.append(f"      → {bl[1]}{marca}")
            elif t == "ver":
                _, e, sub, eta, ok = bl
                L += [f"      VERIFICAÇÃO:  η = {e} = {sub} = {n(eta, 3)}  ({n(eta * 100, 1)} %)  "
                      f"{'✔ ATENDE' if ok else '✘ NÃO ATENDE'}"]
            elif t == "tab":
                _, cab, linhas, al = bl
                larg = [max(len(str(c)), *(len(str(l[i])) for l in linhas)) if linhas else len(str(c))
                        for i, c in enumerate(cab)]

                def fmt_l(vals):
                    return "    " + "  ".join(str(v).rjust(larg[i]) if al[i] == "r" else str(v).ljust(larg[i])
                                              for i, v in enumerate(vals))
                L.append(fmt_l(cab))
                L.append("    " + "  ".join("-" * w for w in larg))
                L += [fmt_l(l) for l in linhas]
            elif t == "aviso":
                L.append(f"    ⚠ {bl[1]}")
        return "\n".join(L)

    # ------------------------------------------------------------------ HTML
    def html(self):
        e = html.escape
        H = []
        for bl in self.b:
            t = bl[0]
            if t == "h1":
                H.append(f'<h3><span>{e(bl[1])}</span>' + (f' <small>({e(bl[2])})</small>' if bl[2] else "") + "</h3>")
            elif t == "h2":
                H.append(f"<h4>{e(bl[1])}</h4>")
            elif t == "p":
                H.append(f"<p>{e(bl[1])}</p>")
            elif t == "eq":
                _, s, ex, sub, v, un, nota = bl
                partes = [f"<b>{e(s)}</b> = {e(ex)}" if ex else f"<b>{e(s)}</b>"]
                if sub:
                    partes.append(f"= {e(sub)}")
                val = f"<b>{e(str(v))}</b> {e(un)}" if v is not None else ""
                H.append(f'<div class="eq"><span class="f">{" ".join(partes)}</span>'
                         f'<span class="v">{val}</span></div>' + (f'<div class="nota">{e(nota)}</div>' if nota else ""))
            elif t == "cond":
                cls = "" if bl[2] is None else (" ok" if bl[2] else " nok")
                marca = "" if bl[2] is None else (" ✔" if bl[2] else " ✘")
                H.append(f'<div class="cond{cls}">→ {e(bl[1])}{marca}</div>')
            elif t == "ver":
                _, ex, sub, eta, ok = bl
                H.append(f'<div class="ver {"ok" if ok else "nok"}"><span>η = {e(ex)} = {e(sub)}</span>'
                         f'<span><b>η = {n(eta, 3)} ({n(eta * 100, 1)} %)</b> {"✔ ATENDE" if ok else "✘ NÃO ATENDE"}'
                         f'</span></div>')
            elif t == "tab":
                _, cab, linhas, al = bl
                H.append("<table><tr>" + "".join(f"<th>{e(str(c))}</th>" for c in cab) + "</tr>")
                for l in linhas:
                    cls = ""
                    if any(str(x) in ("✘", "NÃO PASSA") for x in l):
                        cls = ' class="nok"'
                    H.append(f"<tr{cls}>" + "".join(
                        f'<td class="{"r" if al[i] == "r" else ""}">{e(str(x))}</td>' for i, x in enumerate(l)) + "</tr>")
                H.append("</table>")
            elif t == "aviso":
                H.append(f'<div class="aviso">⚠ {e(bl[1])}</div>')
        return "\n".join(H)


CSS = """
body{font-family:'Segoe UI',Arial,sans-serif;font-size:12.5px;color:#222;max-width:900px;margin:18px auto;line-height:1.45}
h1{font-size:17px;border-bottom:3px solid #1b5e20;padding-bottom:4px;margin-bottom:4px}
h2{font-size:15px;background:#1b5e20;color:#fff;padding:5px 8px;margin-top:26px}
h3{font-size:13px;background:#dcedc8;border-left:5px solid #558b2f;padding:4px 8px;margin:20px 0 6px}
h3 small{font-weight:normal;color:#33691e}
h4{font-size:12.5px;margin:10px 0 3px 8px;color:#1b5e20}
p{margin:3px 0 3px 18px}
.eq{display:flex;justify-content:space-between;margin:2px 0 2px 26px;font-family:Consolas,monospace;font-size:12px}
.eq .v{white-space:nowrap;padding-left:18px}
.nota{margin:0 0 3px 46px;color:#666;font-size:11px}
.cond{margin:2px 0 2px 26px;font-style:italic;color:#333}
.cond.ok{color:#1b5e20}.cond.nok{color:#8b0000}
.ver{display:flex;justify-content:space-between;margin:6px 0 10px 18px;padding:5px 8px;border:1px solid #9ccc65;background:#f1f8e9;font-family:Consolas,monospace}
.ver.nok{border-color:#e57373;background:#ffebee;color:#8b0000}
table{border-collapse:collapse;margin:6px 0 8px 18px;font-size:11.5px}
th{background:#eceff1;border:1px solid #b0bec5;padding:3px 6px;text-align:left}
td{border:1px solid #cfd8dc;padding:2px 6px}td.r{text-align:right}
tr.nok td{color:#8b0000;font-weight:bold}
.aviso{margin:4px 0 4px 18px;padding:4px 8px;background:#fff8e1;border-left:4px solid #ffb300}
.estado{font-size:15px;font-weight:bold;padding:6px 10px;display:inline-block;margin:6px 0}
.estado.ok{background:#1b8a2f;color:#fff}.estado.nok{background:#8b0000;color:#fff}
.barra{page-break-before:always}.barra:first-of-type{page-break-before:auto}
@media print{body{max-width:none;margin:10mm}h2{-webkit-print-color-adjust:exact;print-color-adjust:exact}}
"""


def documento_html(titulo, corpo):
    return (f"<!DOCTYPE html><html lang='pt-BR'><head><meta charset='utf-8'><title>{html.escape(titulo)}</title>"
            f"<style>{CSS}</style></head><body>{corpo}</body></html>")


# =============================================================================
# PROPRIEDADES GEOMÉTRICAS DO PERFIL I
# =============================================================================
def props_I(secao: core.Secao, c):
    """Propriedades em cm (cm², cm³, cm⁴, cm⁶) + fontes. Laminado (W) ou soldado (I)."""
    d_ = secao.dims
    fontes = {}
    if secao.familia == "W":
        d, bf, tw, tf, A, Ix, Iy = core.CATALOGO_W[d_.get("catalogo", "W200x15,0")]
        laminado = True
        fontes["A, Ix, Iy"] = "catálogo"
        A_geo = (2 * bf * tf + (d - 2 * tf) * tw) / 100
        if "r" in d_:
            r = float(d_["r"])
            fontes["r"] = "informado"
        else:
            r = math.sqrt(max(A - A_geo, 0.0) * 100 / (4 * (1 - math.pi / 4)))
            fontes["r"] = "derivado da área do catálogo (A − A sem raios)"
    else:
        d, bf, tw, tf = d_["d"], d_["bf"], d_["tw"], d_["tf"]
        A, _, _, Ix, Iy = core.props_I(d, bf, tw, tf)
        A, Ix, Iy = A / 100, Ix / 1e4, Iy / 1e4
        laminado = False
        r = 0.0
        fontes["A, Ix, Iy"] = "calculado (chapas)"
    af = (1 - math.pi / 4) * r * r            # mm² por raio
    ef = r * (10 - 3 * math.pi) / (12 - 3 * math.pi)
    h_geo = d - 2 * tf
    if laminado and c["h_alma"] == "nbr":
        h = h_geo - 2 * r
        fontes["h"] = "d − 2tf − 2r (NBR 8800, perfis laminados)"
    else:
        h = h_geo
        fontes["h"] = "d − 2tf (faces internas das mesas)"
    Zx = (bf * tf * (d - tf) + tw * h_geo ** 2 / 4 + 4 * af * (d / 2 - tf - ef)) / 1e3
    Zy = (tf * bf ** 2 / 2 + h_geo * tw ** 2 / 4 + 4 * af * (tw / 2 + ef)) / 1e3
    J = (2 * bf * tf ** 3 + (d - tf) * tw ** 3) / 3 / 1e4
    Cw = Iy * ((d - tf) / 10) ** 2 / 4
    fontes.update({"Zx": "calculado", "Zy": "calculado", "J": "calculado (sem raios — conservador)",
                   "Cw": "Iy·(d − tf)²/4"})
    for k in ("Zx", "Zy", "J", "Cw"):
        if d_.get(k) not in (None, "", 0, 0.0):
            fontes[k] = "informado"
    Zx = float(d_.get("Zx") or Zx)
    Zy = float(d_.get("Zy") or Zy)
    J = float(d_.get("J") or J)
    Cw = float(d_.get("Cw") or Cw)
    return {"d": d / 10, "bf": bf / 10, "tw": tw / 10, "tf": tf / 10, "r": r / 10, "h": h / 10,
            "A": A, "Ix": Ix, "Iy": Iy, "Wx": Ix / (d / 20), "Wy": Iy / (bf / 20),
            "Zx": Zx, "Zy": Zy, "rx": math.sqrt(Ix / A), "ry": math.sqrt(Iy / A),
            "J": J, "Cw": Cw, "laminado": laminado, "fontes": fontes}

# =============================================================================
# RESISTÊNCIAS (cada função devolve o valor e o passo a passo em um Memorial)
# =============================================================================
def _kNm(v_kNcm, nd=2):
    return f"{n(v_kNcm, 1)} kN·cm = {n(v_kNcm / 100, nd)} kN·m"


def _regime(m, nome, lam, lp, lr, Mpl, Mr, Mcr, Cb=None):
    """Classifica λ e escreve a expressão do momento resistente nominal (antes de γa1)."""
    if lam <= lp:
        m.cond(f"λ = {n(lam)} ≤ λp = {n(lp)} → seção compacta para {nome}")
        m.eq(f"M,{nome}", "Mpl", None, _kNm(Mpl))
        return Mpl
    if lam <= lr:
        m.cond(f"λp = {n(lp)} < λ = {n(lam)} ≤ λr = {n(lr)} → comportamento inelástico")
        base = Mpl - (Mpl - Mr) * (lam - lp) / (lr - lp)
        if Cb is not None:
            M = Cb * base
            m.eq(f"M,{nome}", "Cb·[Mpl − (Mpl − Mr)·(λ − λp)/(λr − λp)] ≤ Mpl",
                 f"{n(Cb, 3)} × [{n(Mpl, 1)} − ({n(Mpl, 1)} − {n(Mr, 1)})·({n(lam)} − {n(lp)})/({n(lr)} − {n(lp)})]",
                 _kNm(min(M, Mpl)), nota="limitado a Mpl" if M > Mpl else "")
            return min(M, Mpl)
        m.eq(f"M,{nome}", "Mpl − (Mpl − Mr)·(λ − λp)/(λr − λp)",
             f"{n(Mpl, 1)} − ({n(Mpl, 1)} − {n(Mr, 1)})·({n(lam)} − {n(lp)})/({n(lr)} − {n(lp)})", _kNm(base))
        return base
    m.cond(f"λ = {n(lam)} > λr = {n(lr)} → regime elástico: M = Mcr ≤ Mpl")
    M = min(Mcr, Mpl)
    m.eq(f"M,{nome}", "Mcr ≤ Mpl", None, _kNm(M), nota="limitado a Mpl" if Mcr > Mpl else "")
    return M


def tracao(p, mat, c, pj):
    m = Memorial()
    fy, fu, A = mat["fy"], mat["fu"], p["A"]
    m.h2("(a) Escoamento da seção bruta")
    N1 = A * fy / c["ga1"]
    m.eq("Nt,Rd,1", "Ag·fy / γa1", f"{n(A)} × {n(fy)} / {n(c['ga1'])}", n(N1), "kN")
    m.h2("(b) Ruptura da seção líquida efetiva (5.2.3 a 5.2.5)")
    An = pj["an"] * A
    Ae = pj["ct"] * An
    N2 = Ae * fu / c["ga2"]
    m.eq("An", "(An/Ag)·Ag", f"{n(pj['an'], 3)} × {n(A)}", n(An), "cm²")
    m.eq("Ae", "Ct·An", f"{n(pj['ct'], 3)} × {n(An)}", n(Ae), "cm²")
    m.eq("Nt,Rd,2", "Ae·fu / γa2", f"{n(Ae)} × {n(fu)} / {n(c['ga2'])}", n(N2), "kN")
    NtRd = min(N1, N2)
    m.eq("Nt,Rd", "mín(Nt,Rd,1 ; Nt,Rd,2)", f"mín({n(N1)} ; {n(N2)})", n(NtRd), "kN")
    return NtRd, m, {"N1": N1, "N2": N2}


def compressao(p, mat, c, KLx, KLy, KLz):
    E, G, fy, A = mat["E"], mat["G"], mat["fy"], p["A"]
    m = Memorial()
    pi2 = math.pi ** 2
    m.h2("Força axial de flambagem elástica — Anexo E (seção com dois eixos de simetria)")
    if KLx > 0:
        Nex = pi2 * E * p["Ix"] / KLx ** 2
        m.eq("Nex", "π²·E·Ix / (KxLx)²", f"π² × {n(E, 0)} × {n(p['Ix'], 1)} / {n(KLx, 1)}²", n(Nex), "kN")
    else:
        Nex = math.inf
        m.cond("KxLx = 0 → Nex = ∞")
    if KLy > 0:
        Ney = pi2 * E * p["Iy"] / KLy ** 2
        m.eq("Ney", "π²·E·Iy / (KyLy)²", f"π² × {n(E, 0)} × {n(p['Iy'], 1)} / {n(KLy, 1)}²", n(Ney), "kN")
    else:
        Ney = math.inf
        m.cond("KyLy = 0 → Ney = ∞")
    r0 = math.sqrt(p["rx"] ** 2 + p["ry"] ** 2)
    m.eq("r0", "√(rx² + ry² + x0² + y0²)", f"√({n(p['rx'])}² + {n(p['ry'])}² + 0 + 0)", n(r0), "cm",
         nota="centro de cisalhamento coincide com o CG (x0 = y0 = 0)")
    if KLz > 0:
        Nez = (pi2 * E * p["Cw"] / KLz ** 2 + G * p["J"]) / r0 ** 2
        m.eq("Nez", "[π²·E·Cw / (KzLz)² + G·J] / r0²",
             f"[π² × {n(E, 0)} × {n(p['Cw'], 0)} / {n(KLz, 1)}² + {n(G, 0)} × {n(p['J'], 3)}] / {n(r0)}²", n(Nez), "kN")
    else:
        Nez = math.inf
        m.cond("KzLz = 0 → Nez = ∞")
    Ne = min(Nex, Ney, Nez)
    m.eq("Ne", "mín(Nex ; Ney ; Nez)", f"mín({n(Nex)} ; {n(Ney)} ; {n(Nez)})", n(Ne), "kN")

    m.h2("Flambagem local — Anexo F")
    rq = math.sqrt(E / fy)
    bt = p["bf"] / 2 / p["tf"]
    if p["laminado"]:
        m.p("Mesas: elementos AL — Grupo 4 da Tabela F.1 (mesas de perfis I laminados)")
        m.eq("b/t", "(bf/2) / tf", f"({n(p['bf'] * 10, 1)}/2) / {n(p['tf'] * 10)}", n(bt))
        l1, l2 = 0.56 * rq, 1.03 * rq
        m.eq("(b/t)lim", "0,56·√(E/fy)", f"0,56 × √({n(E, 0)}/{n(fy)})", n(l1))
        if bt <= l1:
            Qs = 1.0
            m.cond(f"b/t = {n(bt)} ≤ {n(l1)} → Qs = 1,000")
        elif bt <= l2:
            Qs = 1.415 - 0.74 * bt * math.sqrt(fy / E)
            m.cond(f"0,56√(E/fy) = {n(l1)} < b/t ≤ 1,03√(E/fy) = {n(l2)}")
            m.eq("Qs", "1,415 − 0,74·(b/t)·√(fy/E)", f"1,415 − 0,74 × {n(bt)} × √({n(fy)}/{n(E, 0)})", n(Qs, 3))
        else:
            Qs = 0.69 * E / (fy * bt ** 2)
            m.cond(f"b/t = {n(bt)} > 1,03√(E/fy) = {n(l2)}")
            m.eq("Qs", "0,69·E / [fy·(b/t)²]", f"0,69 × {n(E, 0)} / [{n(fy)} × {n(bt)}²]", n(Qs, 3))
    else:
        m.p("Mesas: elementos AL — Grupo 5 da Tabela F.1 (mesas de perfis I soldados)")
        kc0 = 4 / math.sqrt(p["h"] / p["tw"])
        kc = min(max(kc0, 0.35), 0.76)
        m.eq("kc", "4 / √(h/tw)  (0,35 ≤ kc ≤ 0,76)", f"4 / √({n(p['h'] / p['tw'])})", n(kc, 3))
        m.eq("b/t", "(bf/2) / tf", f"({n(p['bf'] * 10, 1)}/2) / {n(p['tf'] * 10)}", n(bt))
        l1, l2 = 0.64 * math.sqrt(E * kc / fy), 1.17 * math.sqrt(E * kc / fy)
        m.eq("(b/t)lim", "0,64·√[E/(fy/kc)]", None, n(l1))
        if bt <= l1:
            Qs = 1.0
            m.cond(f"b/t = {n(bt)} ≤ {n(l1)} → Qs = 1,000")
        elif bt <= l2:
            Qs = 1.415 - 0.65 * bt * math.sqrt(fy / (kc * E))
            m.cond(f"{n(l1)} < b/t ≤ 1,17√[E/(fy/kc)] = {n(l2)}")
            m.eq("Qs", "1,415 − 0,65·(b/t)·√[fy/(kc·E)]", None, n(Qs, 3))
        else:
            Qs = 0.90 * E * kc / (fy * bt ** 2)
            m.cond(f"b/t = {n(bt)} > {n(l2)}")
            m.eq("Qs", "0,90·E·kc / [fy·(b/t)²]", None, n(Qs, 3))
    m.p("Alma: elemento AA — Grupo 2 da Tabela F.1")
    btw = p["h"] / p["tw"]
    lim = 1.49 * rq
    m.eq("b/t", "h / tw", f"{n(p['h'] * 10, 1)} / {n(p['tw'] * 10)}", n(btw))
    m.eq("(b/t)lim", "1,49·√(E/fy)", f"1,49 × √({n(E, 0)}/{n(fy)})", n(lim))
    if btw <= lim:
        Qa = 1.0
        m.cond(f"b/t = {n(btw)} ≤ {n(lim)} → alma não esbelta: Qa = 1,000")
    else:
        m.cond(f"b/t = {n(btw)} > {n(lim)} → alma esbelta: largura efetiva (F.3.2)")
        if c["sigma_qa"] == "fy":
            sig = fy
            m.eq("σ", "fy", None, n(sig), "kN/cm²", nota="opção conservadora")
        else:
            l0 = math.sqrt(A * fy / Ne)
            chi1 = 0.658 ** (l0 ** 2) if l0 <= 1.5 else 0.877 / l0 ** 2
            sig = chi1 * fy
            m.eq("λ0 (Q = 1)", "√(Ag·fy / Ne)", f"√({n(A)} × {n(fy)} / {n(Ne)})", n(l0, 3))
            m.eq("χ (Q = 1)", "0,658^(λ0²)" if l0 <= 1.5 else "0,877 / λ0²", None, n(chi1, 3))
            m.eq("σ", "χ·fy", f"{n(chi1, 3)} × {n(fy)}", n(sig, 3), "kN/cm²")
        raiz = math.sqrt(E / sig)
        bef0 = 1.92 * p["tw"] * raiz * (1 - 0.34 / btw * raiz)
        bef = min(bef0, p["h"])
        m.eq("bef", "1,92·t·√(E/σ)·[1 − (ca/(b/t))·√(E/σ)] ≤ b",
             f"1,92 × {n(p['tw'], 3)} × √({n(E, 0)}/{n(sig, 3)}) × [1 − (0,34/{n(btw)}) × √({n(E, 0)}/{n(sig, 3)})]",
             n(bef * 10, 1), "mm", nota=f"ca = 0,34; b = h = {n(p['h'] * 10, 1)} mm" + ("; bef limitado a b" if bef0 > p["h"] else ""))
        Aef = A - (p["h"] - bef) * p["tw"]
        Qa = Aef / A
        m.eq("Aef", "Ag − (b − bef)·t", f"{n(A)} − ({n(p['h'], 3)} − {n(bef, 3)}) × {n(p['tw'], 3)}", n(Aef), "cm²")
        m.eq("Qa", "Aef / Ag", f"{n(Aef)} / {n(A)}", n(Qa, 3))
    Q = Qs * Qa
    m.eq("Q", "Qs·Qa", f"{n(Qs, 3)} × {n(Qa, 3)}", n(Q, 3))

    m.h2("Fator de redução χ — 5.3.3")
    l0 = math.sqrt(Q * A * fy / Ne)
    m.eq("λ0", "√(Q·Ag·fy / Ne)", f"√({n(Q, 3)} × {n(A)} × {n(fy)} / {n(Ne)})", n(l0, 3))
    if l0 <= 1.5:
        chi = 0.658 ** (l0 ** 2)
        m.cond(f"λ0 = {n(l0, 3)} ≤ 1,5")
        m.eq("χ", "0,658^(λ0²)", f"0,658^({n(l0, 3)}²)", n(chi, 3))
    else:
        chi = 0.877 / l0 ** 2
        m.cond(f"λ0 = {n(l0, 3)} > 1,5")
        m.eq("χ", "0,877 / λ0²", f"0,877 / {n(l0, 3)}²", n(chi, 3))
    m.h2("Força axial de compressão resistente de cálculo — 5.3.2")
    NcRd = chi * Q * A * fy / c["ga1"]
    m.eq("Nc,Rd", "χ·Q·Ag·fy / γa1", f"{n(chi, 3)} × {n(Q, 3)} × {n(A)} × {n(fy)} / {n(c['ga1'])}", n(NcRd), "kN")
    return NcRd, m, {"Nex": Nex, "Ney": Ney, "Nez": Nez, "Ne": Ne, "Q": Q, "Qs": Qs, "Qa": Qa, "chi": chi, "l0": l0}


def flexao_x(p, mat, c, Lb, Cb):
    """Flexão em torno do eixo de maior inércia (x). Retorna (MRd kN·cm, Memorial, info)."""
    E, fy = mat["E"], mat["fy"]
    m = Memorial()
    sr = 0.3 * fy
    W, Z = p["Wx"], p["Zx"]
    Mpl = Z * fy
    Mlim = 1.5 * W * fy
    m.h2("Grandezas comuns")
    m.eq("Mpl", "Zx·fy", f"{n(Z)} × {n(fy)}", _kNm(Mpl))
    m.eq("σr", "0,30·fy", f"0,30 × {n(fy)}", n(sr, 3), "kN/cm²", nota="tensão residual")
    m.eq("Mr", "(fy − σr)·Wx", f"({n(fy)} − {n(sr, 3)}) × {n(W)}", _kNm((fy - sr) * W), nota="FLT e FLM")
    cand = []
    # ---- FLT
    m.h2("(a) FLT — flambagem lateral com torção (Anexo G, Tabela G.1)")
    if Lb > 0:
        lam = Lb / p["ry"]
        lp = 1.76 * math.sqrt(E / fy)
        b1 = (fy - sr) * W / (E * p["J"])
        lr = 1.38 * math.sqrt(p["Iy"] * p["J"]) / (p["ry"] * p["J"] * b1) * \
            math.sqrt(1 + math.sqrt(1 + 27 * p["Cw"] * b1 ** 2 / p["Iy"]))
        Mr = (fy - sr) * W
        Mcr = Cb * math.pi ** 2 * E * p["Iy"] / Lb ** 2 * math.sqrt(p["Cw"] / p["Iy"] * (1 + 0.039 * p["J"] * Lb ** 2 / p["Cw"]))
        m.eq("λ", "Lb / ry", f"{n(Lb, 1)} / {n(p['ry'])}", n(lam))
        m.eq("λp", "1,76·√(E/fy)", f"1,76 × √({n(E, 0)}/{n(fy)})", n(lp))
        m.eq("β1", "(fy − σr)·Wx / (E·J)", f"({n(fy)} − {n(sr, 3)}) × {n(W)} / ({n(E, 0)} × {n(p['J'], 3)})",
             n(b1, 5), "cm⁻¹")
        m.eq("λr", "1,38·√(Iy·J) / (ry·J·β1) · √[1 + √(1 + 27·Cw·β1²/Iy)]",
             f"1,38 × √({n(p['Iy'], 1)} × {n(p['J'], 3)}) / ({n(p['ry'])} × {n(p['J'], 3)} × {n(b1, 5)}) × "
             f"√[1 + √(1 + 27 × {n(p['Cw'], 0)} × {n(b1, 5)}² / {n(p['Iy'], 1)})]", n(lr))
        m.eq("Mcr", "Cb·π²·E·Iy / Lb² · √[Cw/Iy · (1 + 0,039·J·Lb²/Cw)]",
             f"{n(Cb, 3)} × π² × {n(E, 0)} × {n(p['Iy'], 1)} / {n(Lb, 1)}² × "
             f"√[{n(p['Cw'], 0)}/{n(p['Iy'], 1)} × (1 + 0,039 × {n(p['J'], 3)} × {n(Lb, 1)}²/{n(p['Cw'], 0)})]",
             _kNm(Mcr))
        cand.append(("FLT", _regime(m, "FLT", lam, lp, lr, Mpl, Mr, Mcr, Cb)))
    else:
        m.cond("Lb = 0 (mesa comprimida contida continuamente) → FLT não se aplica")
    # ---- FLM
    m.h2("(b) FLM — flambagem local da mesa comprimida (Anexo G, Tabela G.1)")
    lam = p["bf"] / 2 / p["tf"]
    lp = 0.38 * math.sqrt(E / fy)
    m.eq("λ", "bf / (2·tf)", f"{n(p['bf'] * 10, 1)} / (2 × {n(p['tf'] * 10)})", n(lam))
    m.eq("λp", "0,38·√(E/fy)", f"0,38 × √({n(E, 0)}/{n(fy)})", n(lp))
    if p["laminado"]:
        lr = 0.83 * math.sqrt(E / (fy - sr))
        Mcr = 0.69 * E * W / lam ** 2
        m.eq("λr", "0,83·√[E/(fy − σr)]", f"0,83 × √[{n(E, 0)}/({n(fy)} − {n(sr, 3)})]", n(lr))
        m.eq("Mcr", "0,69·E·Wx / λ²", f"0,69 × {n(E, 0)} × {n(W)} / {n(lam)}²", _kNm(Mcr))
    else:
        kc = min(max(4 / math.sqrt(p["h"] / p["tw"]), 0.35), 0.76)
        lr = 0.95 * math.sqrt(E / ((fy - sr) / kc))
        Mcr = 0.90 * E * kc * W / lam ** 2
        m.eq("kc", "4/√(h/tw)  (0,35 ≤ kc ≤ 0,76)", None, n(kc, 3))
        m.eq("λr", "0,95·√[E/((fy − σr)/kc)]", None, n(lr))
        m.eq("Mcr", "0,90·E·kc·Wx / λ²", None, _kNm(Mcr))
    cand.append(("FLM", _regime(m, "FLM", lam, lp, lr, Mpl, (fy - sr) * W, Mcr)))
    # ---- FLA
    m.h2("(c) FLA — flambagem local da alma (Anexo G, Tabela G.1)")
    lam = p["h"] / p["tw"]
    lp = 3.76 * math.sqrt(E / fy)
    lr = 5.70 * math.sqrt(E / fy)
    m.eq("λ", "h / tw", f"{n(p['h'] * 10, 1)} / {n(p['tw'] * 10)}", n(lam))
    m.eq("λp", "3,76·√(E/fy)", f"3,76 × √({n(E, 0)}/{n(fy)})", n(lp))
    m.eq("λr", "5,70·√(E/fy)", f"5,70 × √({n(E, 0)}/{n(fy)})", n(lr))
    esbelta = lam > lr
    if esbelta:
        m.cond(f"λ = {n(lam)} > λr = {n(lr)} → viga de alma esbelta (Anexo H)", False)
        m.aviso("Anexo H não implementado — barra considerada NÃO APROVADA.")
    else:
        m.eq("Mr", "fy·Wx", f"{n(fy)} × {n(W)}", _kNm(fy * W))
        cand.append(("FLA", _regime(m, "FLA", lam, lp, lr, Mpl, fy * W, math.inf)))
    m.h2("Momento fletor resistente de cálculo — 5.4.2.2")
    m.eq("M,máx", "1,50·Wx·fy", f"1,50 × {n(W)} × {n(fy)}", _kNm(Mlim), nota="limite para validade da análise elástica")
    cand.append(("1,50·W·fy", Mlim))
    gov = min(cand, key=lambda t: t[1])
    MRd = gov[1] / c["ga1"]
    termos = " ; ".join(f"{n(v / 100)}" for _, v in cand)
    m.eq("MRd", "mín(" + " ; ".join(f"M,{k}" if not k.startswith("1,5") else k for k, _ in cand) + ") / γa1",
         f"mín({termos}) / {n(c['ga1'])}", n(MRd / 100), "kN·m", nota=f"governa: {gov[0]}")
    return MRd, m, {"esbelta": esbelta, "gov": gov[0]}


def flexao_y(p, mat, c):
    E, fy = mat["E"], mat["fy"]
    m = Memorial()
    sr = 0.3 * fy
    W, Z = p["Wy"], p["Zy"]
    Mpl = Z * fy
    Mlim = 1.5 * W * fy
    m.h2("Grandezas comuns (FLT não se aplica à flexão em torno do eixo de menor inércia)")
    m.eq("Mpl", "Zy·fy", f"{n(Z)} × {n(fy)}", _kNm(Mpl))
    m.eq("σr", "0,30·fy", f"0,30 × {n(fy)}", n(sr, 3), "kN/cm²")
    cand = []
    m.h2("(a) FLM — flambagem local das mesas (Anexo G, Tabela G.1)")
    lam = p["bf"] / 2 / p["tf"]
    lp = 0.38 * math.sqrt(E / fy)
    m.eq("λ", "bf / (2·tf)", f"{n(p['bf'] * 10, 1)} / (2 × {n(p['tf'] * 10)})", n(lam))
    m.eq("λp", "0,38·√(E/fy)", None, n(lp))
    if p["laminado"]:
        lr = 0.83 * math.sqrt(E / (fy - sr))
        Mcr = 0.69 * E * W / lam ** 2
        m.eq("λr", "0,83·√[E/(fy − σr)]", None, n(lr))
        m.eq("Mcr", "0,69·E·Wy / λ²", f"0,69 × {n(E, 0)} × {n(W)} / {n(lam)}²", _kNm(Mcr))
    else:
        kc = min(max(4 / math.sqrt(p["h"] / p["tw"]), 0.35), 0.76)
        lr = 0.95 * math.sqrt(E / ((fy - sr) / kc))
        Mcr = 0.90 * E * kc * W / lam ** 2
        m.eq("λr", "0,95·√[E/((fy − σr)/kc)]", None, n(lr))
        m.eq("Mcr", "0,90·E·kc·Wy / λ²", None, _kNm(Mcr))
    m.eq("Mr", "(fy − σr)·Wy", f"({n(fy)} − {n(sr, 3)}) × {n(W)}", _kNm((fy - sr) * W))
    cand.append(("FLM", _regime(m, "FLM", lam, lp, lr, Mpl, (fy - sr) * W, Mcr)))
    m.h2("(b) FLA — flambagem local da alma (Anexo G, Tabela G.1)")
    lam = p["h"] / p["tw"]
    lp = 1.12 * math.sqrt(E / fy)
    lr = 1.40 * math.sqrt(E / fy)
    Wef = W
    m.eq("λ", "h / tw", f"{n(p['h'] * 10, 1)} / {n(p['tw'] * 10)}", n(lam))
    m.eq("λp", "1,12·√(E/fy)", None, n(lp))
    m.eq("λr", "1,40·√(E/fy)", None, n(lr))
    m.eq("Mr", "fy·Wef", f"{n(fy)} × {n(Wef)}", _kNm(fy * Wef), nota="Wef ≈ Wy (alma sobre o eixo neutro)")
    m.eq("Mcr", "Wef²/Wy · fy", f"{n(Wef)}² / {n(W)} × {n(fy)}", _kNm(Wef ** 2 / W * fy))
    cand.append(("FLA", _regime(m, "FLA", lam, lp, lr, Mpl, fy * Wef, Wef ** 2 / W * fy)))
    m.h2("Momento fletor resistente de cálculo — 5.4.2.2")
    m.eq("M,máx", "1,50·Wy·fy", f"1,50 × {n(W)} × {n(fy)}", _kNm(Mlim))
    cand.append(("1,50·W·fy", Mlim))
    gov = min(cand, key=lambda t: t[1])
    MRd = gov[1] / c["ga1"]
    m.eq("MRd", "mín(M,FLM ; M,FLA ; 1,50·W·fy) / γa1",
         f"mín({' ; '.join(n(v / 100) for _, v in cand)}) / {n(c['ga1'])}", n(MRd / 100), "kN·m",
         nota=f"governa: {gov[0]}")
    return MRd, m, {"esbelta": False, "gov": gov[0]}


def cortante(p, mat, c, eixo_forte):
    E, fy = mat["E"], mat["fy"]
    m = Memorial()
    if eixo_forte:
        m.p("Força cortante paralela à alma (flexão em torno de x) — 5.4.3.1, alma sem enrijecedores transversais")
        Aw, lam, kv = p["d"] * p["tw"], p["h"] / p["tw"], 5.0
        m.eq("Aw", "d·tw", f"{n(p['d'], 2)} × {n(p['tw'], 3)}", n(Aw), "cm²")
        m.eq("λ", "h / tw", f"{n(p['h'] * 10, 1)} / {n(p['tw'] * 10)}", n(lam))
        m.eq("kv", "valor para alma sem enrijecedores (ou a/h > 3)", None, "5,0")
    else:
        m.p("Força cortante paralela às mesas (flexão em torno de y) — 5.4.3.5")
        Aw, lam, kv = 2 * p["bf"] * p["tf"], p["bf"] / 2 / p["tf"], 1.2
        m.eq("Aw", "2·bf·tf", f"2 × {n(p['bf'], 2)} × {n(p['tf'], 3)}", n(Aw), "cm²")
        m.eq("λ", "bf / (2·tf)", None, n(lam))
        m.eq("kv", "valor para as mesas (5.4.3.5)", None, "1,2")
    lp = 1.10 * math.sqrt(kv * E / fy)
    lr = 1.37 * math.sqrt(kv * E / fy)
    Vpl = 0.60 * Aw * fy
    m.eq("λp", "1,10·√(kv·E/fy)", f"1,10 × √({n(kv, 1)} × {n(E, 0)}/{n(fy)})", n(lp))
    m.eq("λr", "1,37·√(kv·E/fy)", f"1,37 × √({n(kv, 1)} × {n(E, 0)}/{n(fy)})", n(lr))
    m.eq("Vpl", "0,60·Aw·fy", f"0,60 × {n(Aw)} × {n(fy)}", n(Vpl), "kN")
    if lam <= lp:
        V = Vpl
        m.cond(f"λ = {n(lam)} ≤ λp = {n(lp)}")
        m.eq("VRd", "Vpl / γa1", f"{n(Vpl)} / {n(c['ga1'])}", n(V / c["ga1"]), "kN")
    elif lam <= lr:
        V = lp / lam * Vpl
        m.cond(f"λp = {n(lp)} < λ = {n(lam)} ≤ λr = {n(lr)}")
        m.eq("VRd", "(λp/λ)·Vpl / γa1", f"({n(lp)}/{n(lam)}) × {n(Vpl)} / {n(c['ga1'])}", n(V / c["ga1"]), "kN")
    else:
        V = 1.24 * (lp / lam) ** 2 * Vpl
        m.cond(f"λ = {n(lam)} > λr = {n(lr)}")
        m.eq("VRd", "1,24·(λp/λ)²·Vpl / γa1", f"1,24 × ({n(lp)}/{n(lam)})² × {n(Vpl)} / {n(c['ga1'])}",
             n(V / c["ga1"]), "kN")
    return V / c["ga1"], m


def calcular_cb(xs, M, L, Lb, balanco):
    """Cb (5.4.2.3, Rm = 1): menor valor entre trechos de comprimento Lb (conservador).
    Retorna (Cb, info)."""
    if balanco:
        return 1.0, {"motivo": "barra em balanço → Cb = 1,0 (5.4.2.3)"}
    if Lb > L * (1 + 1e-6):
        return 1.0, {"motivo": "Lb maior que a barra: o diagrama de momentos fora da barra não é conhecido → "
                               "Cb = 1,0 (conservador). Informe Cb manualmente se desejar."}
    nseg = max(1, math.ceil(L / Lb - 1e-9))

    def m_at(x):
        for k in range(1, len(xs)):
            if xs[k] >= x:
                x0, x1 = xs[k - 1], xs[k]
                t = 0 if x1 == x0 else (x - x0) / (x1 - x0)
                return M[k - 1] + (M[k] - M[k - 1]) * t
        return M[-1]
    melhor = None
    for s in range(nseg):
        a, b = L * s / nseg, L * (s + 1) / nseg
        pts = [abs(mm) for x, mm in zip(xs, M) if a - 1e-9 <= x <= b + 1e-9]
        Mmax = max(pts + [abs(m_at(a)), abs(m_at(b))])
        MA, MB, MC = (abs(m_at(a + (b - a) * f)) for f in (0.25, 0.5, 0.75))
        cb = 1.0 if Mmax < 1e-9 else min(12.5 * Mmax / (2.5 * Mmax + 3 * MA + 4 * MB + 3 * MC), 3.0)
        info = {"Mmax": Mmax, "MA": MA, "MB": MB, "MC": MC, "a": a, "b": b, "nseg": nseg, "motivo": ""}
        if melhor is None or cb < melhor[0]:
            melhor = (cb, info)
    return melhor


# =============================================================================
# VERIFICAÇÃO DE UMA BARRA
# =============================================================================
class Verificacao:
    def __init__(self, bid):
        self.bid = bid
        self.status = "NA"        # "OK" | "FALHA" | "NA"
        self.eta_max = 0.0
        self.governa = ""
        self.itens = []
        self.memorial = ""        # texto
        self.mem = None           # Memorial estruturado
        self.titulo = ""
        self.motivo = ""

    def html_corpo(self):
        e = html.escape
        cls = {"OK": "ok", "FALHA": "nok"}.get(self.status, "")
        est = {"OK": "PASSA", "FALHA": "NÃO PASSA", "NA": "NÃO VERIFICADA"}[self.status]
        topo = (f"<h1>{e(self.titulo)}</h1><div class='estado {cls}'>{est} — η máx = {n(self.eta_max * 100, 1)} %"
                f" ({e(self.governa)})</div>")
        return topo + (self.mem.html() if self.mem else f"<p>{e(self.memorial)}</p>")

    def html(self):
        return documento_html(self.titulo, self.html_corpo())


def memoriais_html(modelo, verifs):
    e = html.escape
    linhas = "".join(
        f"<tr class='{'nok' if v.status == 'FALHA' else ''}'><td>{bid}</td><td>{e(modelo.barras[bid].secao)}</td>"
        f"<td class='r'>{n(v.eta_max * 100, 1) if v.status != 'NA' else '—'}</td><td>{e(v.governa or v.motivo)}</td>"
        f"<td>{ {'OK': 'PASSA', 'FALHA': 'NÃO PASSA', 'NA': 'N/A'}[v.status]}</td></tr>"
        for bid, v in sorted(verifs.items()) if bid in modelo.barras)
    corpo = (f"<h1>Memorial de verificação — {e(modelo.titulo)}</h1>"
             "<p>ABNT NBR 8800 — perfis I/W. Resumo das barras verificadas:</p>"
             "<table><tr><th>Barra</th><th>Perfil</th><th>η máx (%)</th><th>Governa</th><th>Estado</th></tr>"
             f"{linhas}</table>")
    for bid, v in sorted(verifs.items()):
        corpo += f"<div class='barra'>{v.html_corpo()}</div>"
    return documento_html(f"Memorial — {modelo.titulo}", corpo)


def _selecionar(modelo, resultados, tipo, filtro):
    out = []
    for k, r in resultados.items():
        if not k.startswith("COMB: "):
            continue
        cb = modelo.combinacoes.get(k[6:])
        if not cb:
            continue
        if filtro == "todas":
            if tipo == "ELU" and cb.tipo.startswith("ELS"):
                continue
            if tipo == "ELS" and cb.tipo == "ELU":
                continue
            out.append((cb.nome, r, cb))
        elif tipo == "ELU" and cb.tipo == "ELU":
            out.append((cb.nome, r, cb))
        elif tipo == "ELS" and (cb.tipo == filtro or (filtro == "ELS" and cb.tipo.startswith("ELS"))):
            out.append((cb.nome, r, cb))
    aviso = ""
    if not out:
        aviso = (f"Nenhuma combinação {tipo} ({filtro}) encontrada — usados os casos de carga isolados "
                 "(valores característicos). Gere as combinações NBR 8800.")
        out = [(k, r, None) for k, r in resultados.items() if not k.startswith("COMB: ")]
    return out, aviso


def comprimentos(modelo, b, pj=None):
    """Comprimentos efetivos (cm) e sua origem: dict chave -> (valor_cm, K ou None, texto_origem)."""
    pj = pj or pproj(b)
    L = modelo.geometria(b)[0] * 100
    out = {}
    for chave, kk, lk in (("p", "kp", "Lp"), ("f", "kf", "Lf"), ("z", "kz", "Lz")):
        if pj.get(lk) not in (None, ""):
            out[chave] = (float(pj[lk]) * 100, None, "informado manualmente")
        else:
            out[chave] = (pj[kk] * L, pj[kk], f"K × L = {n(pj[kk], 3)} × {n(L / 100, 3)} m")
    if pj.get("Lb") not in (None, ""):
        out["b"] = (float(pj["Lb"]) * 100, None, "informado manualmente")
    else:
        out["b"] = (pj["klb"] * L, pj["klb"], f"(Lb/L) × L = {n(pj['klb'], 3)} × {n(L / 100, 3)} m")
    return out, L


def verificar_barra(modelo, resultados, bid):
    v = Verificacao(bid)
    b = modelo.barras[bid]
    s = modelo.secoes[b.secao]
    c = cfg(modelo)
    pj = pproj(b)
    v.titulo = f"Barra {bid} (N{b.ni} → N{b.nj}) — {s.nome}"
    if s.familia not in FAMILIAS_SUPORTADAS:
        v.motivo = (f"Família '{core.FAMILIAS.get(s.familia, s.familia)}' ainda não implementada "
                    "(esta etapa cobre perfis I/W laminados e I soldados).")
        v.memorial = v.motivo
        return v
    beta = b.beta % 360
    if min(abs(beta - a) for a in (0, 90, 180, 270, 360)) > 1e-6:
        v.motivo = "Giro β diferente de 0°/90°/180°/270° — verificação em eixos não principais não implementada."
        v.memorial = v.motivo
        return v
    forte = abs(beta % 180) < 1e-6
    mt = modelo.materiais.get(s.material, core.MATERIAIS_PADRAO["ASTM A36"])
    mat = {"E": mt.E / 10, "G": mt.G / 10, "fy": mt.fy / 10, "fu": mt.fu / 10}
    p = props_I(s, c)
    comp, L_cm = comprimentos(modelo, b, pj)
    KLp, KLf, KLz, Lb = comp["p"][0], comp["f"][0], comp["z"][0], comp["b"][0]
    KLx, KLy = (KLp, KLf) if forte else (KLf, KLp)

    elu, aviso_elu = _selecionar(modelo, resultados, "ELU", c["comb_elu"])
    els, aviso_els = _selecionar(modelo, resultados, "ELS", c["comb_els"])

    NtRd, m_t, _ = tracao(p, mat, c, pj)
    NcRd, m_c, dc = compressao(p, mat, c, KLx, KLy, KLz)
    VRd, m_v = cortante(p, mat, c, forte)
    cb_fixo = None if str(pj["cb"]).lower() == "auto" else float(pj["cb"])

    ext = {k: {"eta": -1.0} for k in ("Nt", "Nc", "M", "V", "NM")}
    comprimida = tracionada = False
    cache = {}
    for nome, r, _ in elu:
        xs, N, V, M = core.esforcos_barra(modelo, r, bid, 40)
        if forte:
            if Lb > 0:
                if cb_fixo:
                    Cb, cbinfo = cb_fixo, {"motivo": "Cb informado manualmente"}
                else:
                    Cb, cbinfo = calcular_cb(xs, M, L_cm / 100, Lb / 100, pj["balanco"])
            else:
                Cb, cbinfo = 1.0, {"motivo": "Lb = 0"}
            chave = round(Cb, 4)
            if chave not in cache:
                cache[chave] = flexao_x(p, mat, c, Lb, Cb)
        else:
            Cb, cbinfo, chave = None, None, "y"
            if chave not in cache:
                cache[chave] = flexao_y(p, mat, c)
        MRd = cache[chave][0]
        for x, Ni, Vi, Mi in zip(xs, N, V, M):
            Mi_cm = abs(Mi) * 100
            if Ni < -1e-6:
                comprimida = True
            elif Ni > 1e-6:
                tracionada = True
            for k, sd, rd in (("Nt", max(Ni, 0.0), NtRd), ("Nc", max(-Ni, 0.0), NcRd),
                              ("V", abs(Vi), VRd), ("M", Mi_cm, MRd)):
                eta = sd / rd if rd else math.inf
                if eta > ext[k]["eta"]:
                    ext[k] = {"eta": eta, "combo": nome, "x": x, "sd": sd, "rd": rd, "Cb": Cb,
                              "cbinfo": cbinfo, "chave": chave}
            NRd = NtRd if Ni >= 0 else NcRd
            rN = abs(Ni) / NRd
            eta = rN + 8 / 9 * Mi_cm / MRd if rN >= 0.2 else rN / 2 + Mi_cm / MRd
            if eta > ext["NM"]["eta"]:
                ext["NM"] = {"eta": eta, "combo": nome, "x": x, "N": Ni, "M": Mi_cm, "NRd": NRd, "MRd": MRd,
                             "rN": rN, "Cb": Cb, "sd": 0.0, "rd": 0.0}
    for k in ext:
        if ext[k]["eta"] < 0:
            ext[k] = {"eta": 0.0, "combo": "—", "x": None, "sd": 0.0, "rd": 0.0, "Cb": None, "cbinfo": None,
                      "chave": None}

    lamx = KLx / p["rx"] if KLx > 0 else 0.0
    lamy = KLy / p["ry"] if KLy > 0 else 0.0
    lam = max(lamx, lamy)
    lim_l = 200.0 if comprimida else (300.0 if tracionada else None)

    # flecha
    if pj.get("Lflecha"):
        Lref, orig_lref = float(pj["Lflecha"]), "vão de referência informado"
    elif pj["balanco"]:
        Lref, orig_lref = 2 * L_cm / 100, "balanço: L = 2 × comprimento da barra"
    else:
        Lref, orig_lref = L_cm / 100, "comprimento da barra"
    dlim = Lref / pj["flecha"] * 1000
    fl = {"eta": 0.0, "combo": "—", "x": 0.0, "d": 0.0}
    for nome, r, _ in els:
        xs, U, Vv = core.deformada_barra(modelo, r, bid, 60)
        if pj["flecha_modo"] == "absoluta":
            k = max(range(len(xs)), key=lambda i: abs(Vv[i]))
            x, d = xs[k], Vv[k] * 1000
        elif pj["balanco"]:
            base = Vv[0] if abs(Vv[0]) <= abs(Vv[-1]) else Vv[-1]
            k = max(range(len(xs)), key=lambda i: abs(Vv[i] - base))
            x, d = xs[k], (Vv[k] - base) * 1000
        else:
            x, d = core.deslocamento_max_barra(modelo, r, bid)
        if abs(d) / dlim > fl["eta"]:
            fl = {"eta": abs(d) / dlim, "combo": nome, "x": x, "d": d}

    eM = "x" if forte else "y"
    v.itens = [
        {"nome": "Esbeltez λ", "artigo": "5.3.4 / 5.2.8.1", "eta": lam / lim_l if lim_l else 0.0,
         "combo": "—", "x": None, "sd": lam, "rd": lim_l or 0.0, "un": ""},
        {"nome": "Tração Nt", "artigo": "5.2", "un": "kN", **ext["Nt"]},
        {"nome": "Compressão Nc", "artigo": "5.3", "un": "kN", **ext["Nc"]},
        {"nome": f"Flexão M{eM}", "artigo": "5.4.2 / Anexo G", "un": "kN·m", **ext["M"]},
        {"nome": f"Cortante V{'y' if forte else 'x'}", "artigo": "5.4.3", "un": "kN", **ext["V"]},
        {"nome": f"Interação N+M{eM}", "artigo": "5.5.1.2", "un": "", **ext["NM"]},
        {"nome": "Flecha (ELS)", "artigo": "Anexo C", "eta": fl["eta"], "combo": fl["combo"], "x": fl["x"],
         "sd": abs(fl["d"]), "rd": dlim, "un": "mm"},
    ]
    esbelta = any(dd[2].get("esbelta") for dd in cache.values())
    v.eta_max = max(it["eta"] for it in v.itens)
    v.governa = max(v.itens, key=lambda it: it["eta"])["nome"]
    v.status = "OK" if v.eta_max <= 1.0 + 1e-9 and not esbelta else "FALHA"
    if esbelta:
        v.motivo = "Alma esbelta (λ > λr) — Anexo H não implementado; barra marcada como não aprovada."

    # =================================================================== MEMORIAL
    M = Memorial()
    lin = lambda it: n(it["x"], 3) if it.get("x") is not None else "—"   # noqa: E731

    M.h1("1. DADOS DA BARRA")
    M.tab(["Item", "Valor"], [
        ["Barra / nós", f"{bid}  (N{b.ni} → N{b.nj})"],
        ["Comprimento L", f"{n(L_cm / 100, 3)} m"],
        ["Perfil", f"{s.nome}  ({'laminado' if p['laminado'] else 'soldado'})"],
        ["Giro do perfil β", f"{beta:g}° → flexão no plano em torno do eixo {'X (maior inércia)' if forte else 'Y (menor inércia)'}"],
        ["Ligação nas extremidades", b.tipo],
    ])
    M.h2("Material e coeficientes de ponderação")
    M.tab(["Grandeza", "Valor", "Unidade de cálculo"], [
        ["fy — resistência ao escoamento", f"{n(mt.fy, 1)} MPa", f"{n(mat['fy'], 3)} kN/cm²"],
        ["fu — resistência à ruptura", f"{n(mt.fu, 1)} MPa", f"{n(mat['fu'], 3)} kN/cm²"],
        ["E — módulo de elasticidade", f"{n(mt.E, 0)} MPa", f"{n(mat['E'], 0)} kN/cm²"],
        ["G — módulo de elasticidade transversal", f"{n(mt.G, 0)} MPa", f"{n(mat['G'], 0)} kN/cm²"],
        ["γa1 — escoamento, flambagem e instabilidade", n(c["ga1"]), ""],
        ["γa2 — ruptura", n(c["ga2"]), ""],
    ], ["l", "r", "r"])
    M.h2("Propriedades geométricas da seção")
    f = p["fontes"]
    M.tab(["Propriedade", "Valor", "Origem"], [
        ["d / bf", f"{n(p['d'] * 10, 1)} / {n(p['bf'] * 10, 1)} mm", "catálogo / dimensões"],
        ["tw / tf", f"{n(p['tw'] * 10)} / {n(p['tf'] * 10)} mm", "catálogo / dimensões"],
        ["r — raio de concordância", f"{n(p['r'] * 10, 1)} mm", f.get("r", "—")],
        ["h — altura da alma (cálculo)", f"{n(p['h'] * 10, 1)} mm", f.get("h", "")],
        ["Ag", f"{n(p['A'])} cm²", f.get("A, Ix, Iy", "")],
        ["Ix / Iy", f"{n(p['Ix'], 1)} / {n(p['Iy'], 1)} cm⁴", f.get("A, Ix, Iy", "")],
        ["Wx / Wy", f"{n(p['Wx'])} / {n(p['Wy'])} cm³", "Ix/(d/2) ; Iy/(bf/2)"],
        ["Zx", f"{n(p['Zx'])} cm³", f.get("Zx", "")],
        ["Zy", f"{n(p['Zy'])} cm³", f.get("Zy", "")],
        ["rx / ry", f"{n(p['rx'])} / {n(p['ry'])} cm", "√(I/Ag)"],
        ["J (It)", f"{n(p['J'], 3)} cm⁴", f.get("J", "")],
        ["Cw", f"{n(p['Cw'], 0)} cm⁶", f.get("Cw", "")],
    ], ["l", "r", "l"])
    M.h2("Comprimentos de flambagem e contenção lateral")
    M.tab(["Direção", "Eixo", "Comprimento efetivo", "Origem"], [
        ["Flexão no plano", "X" if forte else "Y", f"{n(KLp / 100, 3)} m", comp["p"][2]],
        ["Flexão fora do plano", "Y" if forte else "X", f"{n(KLf / 100, 3)} m", comp["f"][2]],
        ["Torção", "Z", f"{n(KLz / 100, 3)} m", comp["z"][2]],
        ["Contenção lateral (FLT)", "—", f"Lb = {n(Lb / 100, 3)} m", comp["b"][2]],
        ["Fator Cb", "—", "automático (5.4.2.3)" if cb_fixo is None else n(cb_fixo, 3),
         "balanço → 1,0" if pj["balanco"] else ""],
    ], ["l", "l", "r", "l"])
    M.p("Unidades nas substituições numéricas: kN e cm (tensões em kN/cm²; momentos em kN·cm, convertidos para kN·m).")
    M.h2("Combinações consideradas")
    M.p(f"Estados-limites últimos: {len(elu)} combinações ({c['comb_elu']}).  "
        f"Estados-limites de serviço: {len(els)} combinações ({c['comb_els']}).")
    for a in (aviso_elu, aviso_els):
        if a:
            M.aviso(a)
    M.p("Modelo plano (2D): momento fletor e cortante fora do plano e torção são nulos. "
        "Esforços obtidos em 41 seções ao longo da barra em cada combinação; x medido a partir do nó inicial.")

    M.h1("2. RESUMO DAS VERIFICAÇÕES")
    linhas = []
    for it in v.itens:
        if it["nome"].startswith("Flexão"):
            sd, rd = f"{n(it['sd'] / 100)} kN·m", f"{n(it['rd'] / 100)} kN·m"
        elif it["nome"].startswith("Interação"):
            sd, rd = "—", "1,000"
        elif it["nome"].startswith("Esbeltez"):
            sd, rd = n(it["sd"], 1), (f"{it['rd']:g}" if it["rd"] else "n.a.")
        else:
            sd, rd = f"{n(it['sd'])} {it['un']}", f"{n(it['rd'])} {it['un']}"
        linhas.append([it["nome"], it["artigo"], sd, rd, n(it["eta"] * 100, 1), lin(it), it["combo"],
                       "✔" if it["eta"] <= 1 + 1e-9 else "✘"])
    M.tab(["Verificação", "NBR 8800", "Solicitante", "Resistente", "η (%)", "x (m)", "Combinação", ""],
          linhas, ["l", "l", "r", "r", "r", "r", "l", "l"])
    M.p(f"Estado: {'PASSA' if v.status == 'OK' else 'NÃO PASSA'} — η máx = {n(v.eta_max * 100, 1)} % ({v.governa})")
    if v.motivo:
        M.aviso(v.motivo)

    # 3. esbeltez
    M.h1("3. LIMITAÇÃO DO ÍNDICE DE ESBELTEZ", "NBR 8800 — 5.3.4 e 5.2.8.1")
    M.eq("λx", "KxLx / rx", f"{n(KLx, 1)} / {n(p['rx'])}", n(lamx, 1))
    M.eq("λy", "KyLy / ry", f"{n(KLy, 1)} / {n(p['ry'])}", n(lamy, 1))
    M.eq("λ", "máx(λx ; λy)", None, n(lam, 1))
    if lim_l:
        M.p("Limite: 200 para barras comprimidas (5.3.4); 300 recomendado para barras tracionadas (5.2.8.1).")
        M.ver("λ / λlim", f"{n(lam, 1)} / {lim_l:g}", lam / lim_l)
    else:
        M.cond("Força axial nula em todas as combinações ELU → limitação não se aplica")

    # 4. tração
    it = ext["Nt"]
    M.h1("4. RESISTÊNCIA À FORÇA AXIAL DE TRAÇÃO", "NBR 8800 — 5.2")
    M.p("Condição: Nt,Sd ≤ Nt,Rd")
    if it["sd"] > 0:
        M.p(f"Solicitante de cálculo desfavorável: Nt,Sd = {n(it['sd'])} kN, em x = {lin(it)} m, combinação {it['combo']}.")
    else:
        M.p("A barra não é tracionada em nenhuma combinação ELU (Nt,Sd = 0).")
    M.extend(m_t)
    M.ver("Nt,Sd / Nt,Rd", f"{n(it['sd'])} / {n(NtRd)}", it["eta"])

    # 5. compressão
    it = ext["Nc"]
    M.h1("5. RESISTÊNCIA À FORÇA AXIAL DE COMPRESSÃO", "NBR 8800 — 5.3, Anexos E e F")
    M.p("Condição: Nc,Sd ≤ Nc,Rd = χ·Q·Ag·fy/γa1")
    if it["sd"] > 0:
        M.p(f"Solicitante de cálculo desfavorável: Nc,Sd = {n(it['sd'])} kN, em x = {lin(it)} m, combinação {it['combo']}.")
    else:
        M.p("A barra não é comprimida em nenhuma combinação ELU (Nc,Sd = 0).")
    M.p(f"Comprimentos: KxLx = {n(KLx, 1)} cm;  KyLy = {n(KLy, 1)} cm;  KzLz = {n(KLz, 1)} cm.")
    M.extend(m_c)
    M.ver("Nc,Sd / Nc,Rd", f"{n(it['sd'])} / {n(NcRd)}", it["eta"])

    # 6. flexão
    it = ext["M"]
    M.h1(f"6. RESISTÊNCIA AO MOMENTO FLETOR — EIXO {eM.upper()}", "NBR 8800 — 5.4.2 e Anexo G")
    M.p("Condição: MSd ≤ MRd")
    M.p(f"Solicitante de cálculo desfavorável: MSd = {n(it['sd'] / 100)} kN·m, em x = {lin(it)} m, "
        f"combinação {it['combo']}.")
    if forte and Lb > 0:
        M.h2("Fator de modificação para diagrama de momento não uniforme — Cb (5.4.2.3)")
        ci = it.get("cbinfo") or {}
        if ci.get("motivo"):
            M.cond(ci["motivo"])
            M.eq("Cb", "valor adotado", None, n(it.get("Cb") or 1.0, 3))
        elif "Mmax" in ci:
            if ci["nseg"] > 1:
                M.p(f"Barra dividida em {ci['nseg']} trechos de comprimento Lb; adotado o menor Cb "
                    f"(trecho de x = {n(ci['a'], 3)} m a {n(ci['b'], 3)} m), na combinação {it['combo']}.")
            else:
                M.p(f"Trecho entre contenções = barra inteira; combinação {it['combo']}.")
            M.tab(["Momento (valor absoluto)", "kN·m"], [
                ["Mmax — máximo no trecho", n(ci["Mmax"])], ["MA — a 1/4 do trecho", n(ci["MA"])],
                ["MB — no centro do trecho", n(ci["MB"])], ["MC — a 3/4 do trecho", n(ci["MC"])]], ["l", "r"])
            M.eq("Cb", "12,5·Mmax / (2,5·Mmax + 3·MA + 4·MB + 3·MC) · Rm ≤ 3,0",
                 f"12,5 × {n(ci['Mmax'])} / (2,5 × {n(ci['Mmax'])} + 3 × {n(ci['MA'])} + 4 × {n(ci['MB'])} + "
                 f"3 × {n(ci['MC'])}) × 1,0", n(it.get("Cb") or 1.0, 3), nota="Rm = 1,0 (seção duplamente simétrica)")
    mflex = cache.get(it.get("chave")) if it.get("chave") is not None else None
    if mflex is None:
        mflex = flexao_x(p, mat, c, Lb, 1.0) if forte else flexao_y(p, mat, c)
    M.extend(mflex[1])
    M.ver("MSd / MRd", f"{n(it['sd'] / 100)} / {n(it['rd'] / 100 if it['rd'] else mflex[0] / 100)}", it["eta"])

    # 7. cortante
    it = ext["V"]
    M.h1(f"7. RESISTÊNCIA À FORÇA CORTANTE — {'ALMA' if forte else 'MESAS'}", "NBR 8800 — 5.4.3")
    M.p("Condição: VSd ≤ VRd")
    M.p(f"Solicitante de cálculo desfavorável: VSd = {n(it['sd'])} kN, em x = {lin(it)} m, combinação {it['combo']}.")
    M.extend(m_v)
    M.ver("VSd / VRd", f"{n(it['sd'])} / {n(VRd)}", it["eta"])

    # 8. interação
    it = ext["NM"]
    M.h1("8. FORÇA AXIAL E MOMENTO FLETOR COMBINADOS", "NBR 8800 — 5.5.1.2")
    if it["combo"] != "—":
        M.p(f"Seção determinante: x = {lin(it)} m, combinação {it['combo']} (esforços simultâneos).")
        M.eq("NSd", "", None, f"{n(abs(it['N']))} kN ({'tração' if it['N'] >= 0 else 'compressão'})")
        M.eq("NRd", "Nt,Rd" if it["N"] >= 0 else "Nc,Rd", None, n(it["NRd"]), "kN")
        M.eq(f"M{eM},Sd", "", None, n(it["M"] / 100), "kN·m")
        M.eq(f"M{eM},Rd", "", None, n(it["MRd"] / 100), "kN·m")
        M.eq("NSd/NRd", f"{n(abs(it['N']))} / {n(it['NRd'])}", None, n(it["rN"], 3))
        if it["rN"] >= 0.2:
            M.cond(f"NSd/NRd = {n(it['rN'], 3)} ≥ 0,2 → equação (a)")
            M.ver(f"NSd/NRd + 8/9·(M{eM},Sd/M{eM},Rd)",
                  f"{n(it['rN'], 3)} + 8/9 × ({n(it['M'] / 100)}/{n(it['MRd'] / 100)})", it["eta"])
        else:
            M.cond(f"NSd/NRd = {n(it['rN'], 3)} < 0,2 → equação (b)")
            M.ver(f"NSd/(2·NRd) + M{eM},Sd/M{eM},Rd",
                  f"{n(it['rN'], 3)}/2 + {n(it['M'] / 100)}/{n(it['MRd'] / 100)}", it["eta"])
    else:
        M.p("Sem esforços.")

    # 9. flecha
    M.h1("9. DESLOCAMENTO MÁXIMO (FLECHA)", "NBR 8800 — Anexo C, Tabela C.1")
    modo = {"absoluta": "absoluto (deslocamento total perpendicular ao eixo da barra, em relação à posição indeformada)",
            "corda": "relativo à corda da barra (descontado o deslocamento das extremidades)"}[pj["flecha_modo"]]
    if pj["balanco"] and pj["flecha_modo"] != "absoluta":
        modo = "relativo à extremidade apoiada do balanço"
    M.p(f"Deslocamento {modo}.")
    M.eq("Lref", orig_lref, None, n(Lref, 3), "m")
    M.eq("δlim", f"Lref / {pj['flecha']:g}", f"{n(Lref * 1000, 0)} / {pj['flecha']:g}", n(dlim), "mm")
    M.p(f"Deslocamento máximo: δ = {n(abs(fl['d']))} mm, em x = {n(fl['x'], 3)} m, combinação {fl['combo']}.")
    M.ver("δ / δlim", f"{n(abs(fl['d']))} / {n(dlim)}", fl["eta"])

    M.h1("10. CONCLUSÃO")
    M.p(f"{'A barra ATENDE a todas as verificações.' if v.status == 'OK' else 'A barra NÃO ATENDE às verificações.'}  "
        f"Maior aproveitamento: η = {n(v.eta_max * 100, 1)} % — {v.governa}.")
    M.p("Referência normativa: ABNT NBR 8800 (numeração de artigos da edição 2008; conferir revisão 2024).")
    v.mem = M
    cab = ["=" * 96, f"MEMORIAL DE VERIFICAÇÃO — {v.titulo}", f"ABNT NBR 8800 — ESTADO: "
           f"{'PASSA' if v.status == 'OK' else 'NÃO PASSA'} (η máx = {n(v.eta_max * 100, 1)} %)", "=" * 96]
    v.memorial = "\n".join(cab) + M.texto()
    return v


def verificar(modelo, resultados, ids=None):
    ids = list(modelo.barras) if ids is None else ids
    return {bid: verificar_barra(modelo, resultados, bid) for bid in ids if bid in modelo.barras}
