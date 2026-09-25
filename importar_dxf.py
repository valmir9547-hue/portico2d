# -*- coding: utf-8 -*-
"""
importar_dxf.py — Importação de barras a partir de DXF (sem GUI)
================================================================
Entidades lidas: LINE, LWPOLYLINE, POLYLINE (2D/3D, cada trecho vira uma barra).
Com a biblioteca `ezdxf` instalada (pip install ezdxf), também explode blocos
(INSERT) e lê DXF binário; sem ela, usa um leitor ASCII próprio (LINE/LWPOLYLINE/POLYLINE).
Coordenada Z é ignorada (desenho no plano XY).  Arcos/círculos não são importados.

Fluxo:  ler_segmentos() → [ (layer, x1, y1, x2, y2) ]  (unidades do desenho)
        importar()      → cria nós/barras no Modelo com escala, tolerância, quebra
                          nas interseções e mapeamento layer → seção/ligação.
"""
from __future__ import annotations

import math
from collections import defaultdict


# =============================================================================
# LEITURA
# =============================================================================
def ler_segmentos(caminho):
    try:
        import ezdxf  # noqa: F401
        return _ler_ezdxf(caminho)
    except ImportError:
        return _ler_ascii(caminho)


def _ler_ezdxf(caminho):
    import ezdxf
    doc = ezdxf.readfile(caminho)
    segs = []

    def trata(e, layer=None):
        lay = layer or e.dxf.layer
        t = e.dxftype()
        if t == "LINE":
            segs.append((lay, e.dxf.start.x, e.dxf.start.y, e.dxf.end.x, e.dxf.end.y))
        elif t == "LWPOLYLINE":
            pts = [(p[0], p[1]) for p in e.get_points("xy")]
            if e.closed and len(pts) > 2:
                pts.append(pts[0])
            segs.extend((lay, *a, *b) for a, b in zip(pts, pts[1:]))
        elif t == "POLYLINE":
            pts = [(v.dxf.location.x, v.dxf.location.y) for v in e.vertices]
            if e.is_closed and len(pts) > 2:
                pts.append(pts[0])
            segs.extend((lay, *a, *b) for a, b in zip(pts, pts[1:]))
        elif t == "INSERT":
            for v in e.virtual_entities():
                # entidades no layer "0" dentro do bloco herdam o layer do INSERT
                trata(v, lay if v.dxf.layer == "0" else v.dxf.layer)

    for e in doc.modelspace():
        trata(e)
    return segs


def _ler_ascii(caminho):
    """Leitor mínimo de DXF ASCII: pares (código, valor) na seção ENTITIES."""
    with open(caminho, "r", encoding="utf-8", errors="ignore") as f:
        linhas = [l.strip() for l in f]
    pares = list(zip(linhas[0::2], linhas[1::2]))
    segs = []
    em_ent = False
    i = 0
    atual = None          # (tipo, dict de listas)
    poly = None           # POLYLINE em montagem: [layer, pts, fechado]

    def fecha(ent):
        nonlocal poly
        if ent is None:
            return
        tipo, d = ent
        lay = d.get("8", ["0"])[0]
        if tipo == "LINE":
            try:
                segs.append((lay, float(d["10"][0]), float(d["20"][0]), float(d["11"][0]), float(d["21"][0])))
            except (KeyError, ValueError):
                pass
        elif tipo == "LWPOLYLINE":
            xs, ys = d.get("10", []), d.get("20", [])
            pts = [(float(x), float(y)) for x, y in zip(xs, ys)]
            if int(d.get("70", ["0"])[0]) & 1 and len(pts) > 2:
                pts.append(pts[0])
            segs.extend((lay, *a, *b) for a, b in zip(pts, pts[1:]))
        elif tipo == "POLYLINE":
            poly = [lay, [], bool(int(d.get("70", ["0"])[0]) & 1)]
        elif tipo == "VERTEX" and poly is not None:
            try:
                poly[1].append((float(d["10"][0]), float(d["20"][0])))
            except (KeyError, ValueError):
                pass
        elif tipo == "SEQEND" and poly is not None:
            lay, pts, fechado = poly
            if fechado and len(pts) > 2:
                pts.append(pts[0])
            segs.extend((lay, *a, *b) for a, b in zip(pts, pts[1:]))
            poly = None

    while i < len(pares):
        cod, val = pares[i]
        if cod == "0" and val == "SECTION" and i + 1 < len(pares):
            em_ent = pares[i + 1][1] == "ENTITIES"
        elif cod == "0" and val == "ENDSEC":
            fecha(atual)
            atual = None
            em_ent = False
        elif em_ent and cod == "0":
            fecha(atual)
            atual = (val, defaultdict(list))
        elif em_ent and atual is not None:
            atual[1][cod].append(val)
        i += 1
    fecha(atual)
    return segs


def resumo_layers(segs):
    cont = defaultdict(int)
    for s in segs:
        cont[s[0]] += 1
    return dict(sorted(cont.items()))


def extensao(segs):
    if not segs:
        return 0.0, 0.0, 0.0, 0.0
    xs = [v for s in segs for v in (s[1], s[3])]
    ys = [v for s in segs for v in (s[2], s[4])]
    return min(xs), min(ys), max(xs), max(ys)


# =============================================================================
# GEOMETRIA — quebra nas interseções
# =============================================================================
def _quebrar(segs, tol):
    """Divide cada segmento nos pontos onde outro segmento o cruza ou toca (T, X, L)."""
    n = len(segs)
    cortes = [[0.0, 1.0] for _ in range(n)]
    caixas = [(min(s[1], s[3]) - tol, min(s[2], s[4]) - tol, max(s[1], s[3]) + tol, max(s[2], s[4]) + tol)
              for s in segs]
    # índice espacial simples por células
    cel = max(tol * 50, 1e-9)
    xs = [c[2] - c[0] for c in caixas] + [c[3] - c[1] for c in caixas]
    if xs:
        cel = max(cel, sorted(xs)[len(xs) // 2])
    grade = defaultdict(list)
    for k, (x0, y0, x1, y1) in enumerate(caixas):
        for gx in range(int(math.floor(x0 / cel)), int(math.floor(x1 / cel)) + 1):
            for gy in range(int(math.floor(y0 / cel)), int(math.floor(y1 / cel)) + 1):
                grade[(gx, gy)].append(k)
    pares = set()
    for lst in grade.values():
        for a in range(len(lst)):
            for b in range(a + 1, len(lst)):
                pares.add((lst[a], lst[b]) if lst[a] < lst[b] else (lst[b], lst[a]))
    for a, b in pares:
        A, B = caixas[a], caixas[b]
        if A[2] < B[0] or B[2] < A[0] or A[3] < B[1] or B[3] < A[1]:
            continue
        _, x1, y1, x2, y2 = segs[a]
        _, x3, y3, x4, y4 = segs[b]
        r = (x2 - x1, y2 - y1)
        s = (x4 - x3, y4 - y3)
        La, Lb = math.hypot(*r), math.hypot(*s)
        den = r[0] * s[1] - r[1] * s[0]
        if abs(den) > 1e-12 * La * Lb:           # não paralelos
            qp = (x3 - x1, y3 - y1)
            t = (qp[0] * s[1] - qp[1] * s[0]) / den
            u = (qp[0] * r[1] - qp[1] * r[0]) / den
            ta, tb = tol / La, tol / Lb
            if -ta <= t <= 1 + ta and -tb <= u <= 1 + tb:
                cortes[a].append(min(max(t, 0.0), 1.0))
                cortes[b].append(min(max(u, 0.0), 1.0))
        else:                                    # colineares: extremidades de um sobre o outro
            for (px, py), k, (ox, oy), v, L in (((x3, y3), a, (x1, y1), r, La), ((x4, y4), a, (x1, y1), r, La),
                                                ((x1, y1), b, (x3, y3), s, Lb), ((x2, y2), b, (x3, y3), s, Lb)):
                t = ((px - ox) * v[0] + (py - oy) * v[1]) / (L * L)
                dist = abs((px - ox) * v[1] - (py - oy) * v[0]) / L
                if dist <= tol and 0 < t < 1:
                    cortes[k].append(t)
    saida = []
    for (lay, x1, y1, x2, y2), ts in zip(segs, cortes):
        ts = sorted(ts)
        L = math.hypot(x2 - x1, y2 - y1)
        limpos = [ts[0]]
        for t in ts[1:]:
            if (t - limpos[-1]) * L > tol:
                limpos.append(t)
        if limpos[-1] < 1.0:
            limpos[-1] = 1.0
        for t0, t1 in zip(limpos, limpos[1:]):
            saida.append((lay, x1 + (x2 - x1) * t0, y1 + (y2 - y1) * t0, x1 + (x2 - x1) * t1, y1 + (y2 - y1) * t1))
    return saida


# =============================================================================
# IMPORTAÇÃO NO MODELO
# =============================================================================
def importar(modelo, segs, escala=1.0, mapa=None, tol=1e-3, quebrar=True, origem=True, casas=6):
    """
    segs   : [(layer, x1, y1, x2, y2)] em unidades do desenho
    escala : fator desenho → m (mm: 0.001, cm: 0.01, m: 1)
    mapa   : {layer: (nome_secao | None, rot_i, rot_j, beta)}   (None = ignorar layer)
    tol    : tolerância de fusão de nós e de quebra (m)
    origem : transladar para que o menor ponto fique em (0, 0)
    Retorna relatório (dict).
    """
    mapa = mapa or {}
    usados = [s for s in segs if mapa.get(s[0], (None,))[0]]
    rel = {"lidos": len(segs), "ignorados_layer": len(segs) - len(usados), "nulos": 0,
           "duplicados": 0, "nos_novos": 0, "barras_novas": 0}
    if not usados:
        return rel
    esc = [(l, x1 * escala, y1 * escala, x2 * escala, y2 * escala) for l, x1, y1, x2, y2 in usados]
    if origem:
        x0, y0, _, _ = extensao(esc)
        esc = [(l, x1 - x0, y1 - y0, x2 - x0, y2 - y0) for l, x1, y1, x2, y2 in esc]
    esc2 = []
    for s in esc:
        if math.hypot(s[3] - s[1], s[4] - s[2]) <= tol:
            rel["nulos"] += 1
        else:
            esc2.append(s)
    if quebrar:
        esc2 = _quebrar(esc2, tol)

    # fusão de nós por hash em grade de tamanho tol (verifica células vizinhas)
    grade = defaultdict(list)
    for n in modelo.nos.values():
        grade[(round(n.x / tol), round(n.y / tol))].append(n)

    def no_para(x, y):
        gx, gy = round(x / tol), round(y / tol)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for n in grade.get((gx + dx, gy + dy), []):
                    if math.hypot(n.x - x, n.y - y) <= tol:
                        return n
        antes = len(modelo.nos)
        n = modelo.add_no(round(x, casas), round(y, casas))
        if len(modelo.nos) > antes:
            rel["nos_novos"] += 1
        grade[(gx, gy)].append(n)
        return n

    existentes = {frozenset((b.ni, b.nj)) for b in modelo.barras.values()}
    for lay, x1, y1, x2, y2 in esc2:
        sec, ri, rj, beta = mapa[lay]
        a, b = no_para(x1, y1), no_para(x2, y2)
        if a.id == b.id:
            rel["nulos"] += 1
            continue
        k = frozenset((a.id, b.id))
        if k in existentes:
            rel["duplicados"] += 1
            continue
        bb = modelo.add_barra(a.id, b.id, sec, ri, rj)
        bb.beta = beta
        existentes.add(k)
        rel["barras_novas"] += 1
    return rel
