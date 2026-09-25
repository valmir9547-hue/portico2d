# -*- coding: utf-8 -*-
"""
portico2d_core.py — Motor de cálculo de pórticos, vigas e treliças planas (2D)
==============================================================================
Arquivo SEM dependências de interface gráfica (pode ser testado isoladamente).

Unidades internas:  kN, m, kN·m  (seções digitadas em mm, materiais em MPa)
Análise:           OpenSeesPy — elasticBeamColumn (Euler-Bernoulli), linear, 1ª ordem.

Convenções de sinais (padrão Ftool / prática brasileira):
  - Eixo local x: do nó inicial (i) para o nó final (j); eixo local y: x girado 90° anti-horário.
  - N > 0  → tração;  N < 0 → compressão.
  - V > 0  → resultante das forças à esquerda da seção aponta para +y local  (dM/dx = V).
  - M > 0  → traciona a face do lado −y local (em viga horizontal i→j: traciona embaixo).
    O diagrama de M é desenhado do lado tracionado.
  - Reações nos eixos globais (X → direita, Y → para cima, Mz anti-horário +).
  - Giro do perfil β (graus, em torno do eixo x local): I no plano = Ix·cos²β + Iy·sin²β.

Rótulas nas extremidades (liberação de momento) são modeladas com um nó auxiliar
ligado ao nó principal por equalDOF (ux, uy). Barra de treliça = rótula nas duas pontas.
Apoios elásticos (molas) são modelados com zeroLength ligado a um nó terra engastado.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict

# =============================================================================
# MATERIAIS
# =============================================================================
@dataclass
class Material:
    nome: str
    E: float = 200000.0      # MPa
    G: float = 77000.0       # MPa
    fy: float = 250.0        # MPa
    fu: float = 400.0        # MPa
    gamma: float = 78.5      # kN/m³


MATERIAIS_PADRAO = {
    "ASTM A36": Material("ASTM A36", 200000.0, 77000.0, 250.0, 400.0, 78.5),
    "ASTM A572 Gr50": Material("ASTM A572 Gr50", 200000.0, 77000.0, 345.0, 450.0, 78.5),
}

# =============================================================================
# CATÁLOGOS DE PERFIS (valores de referência — CONFERIR com o catálogo do fabricante)
# =============================================================================
# Perfis W laminados (Gerdau): nome: (d, bf, tw, tf [mm], A [cm²], Ix [cm4], Iy [cm4])
CATALOGO_W = {
    "W150x13,0":  (148, 100, 4.3, 4.9, 16.6, 635, 82),
    "W150x18,0":  (153, 102, 5.8, 7.1, 23.4, 939, 126),
    "W200x15,0":  (200, 100, 4.3, 5.2, 19.4, 1305, 87),
    "W200x19,3":  (203, 102, 5.8, 6.5, 25.1, 1686, 116),
    "W200x22,5":  (206, 102, 6.2, 8.0, 29.0, 2029, 142),
    "W200x26,6":  (207, 133, 5.8, 8.4, 34.2, 2611, 330),
    "W200x31,3":  (210, 134, 6.4, 10.2, 40.3, 3168, 410),
    "W250x17,9":  (251, 101, 4.8, 5.3, 23.1, 2291, 91),
    "W250x22,3":  (254, 102, 5.8, 6.9, 28.9, 2939, 123),
    "W250x25,3":  (257, 102, 6.1, 8.4, 32.6, 3473, 149),
    "W250x28,4":  (260, 102, 6.4, 10.0, 36.6, 4046, 178),
    "W250x32,7":  (258, 146, 6.1, 9.1, 42.1, 4937, 473),
    "W310x21,0":  (303, 101, 5.1, 5.7, 27.2, 3776, 98),
    "W310x23,8":  (305, 101, 5.6, 6.7, 30.7, 4346, 116),
    "W310x28,3":  (309, 102, 6.0, 8.9, 36.5, 5500, 158),
    "W310x32,7":  (313, 102, 6.6, 10.8, 42.1, 6570, 192),
    "W360x32,9":  (349, 127, 5.8, 8.5, 42.1, 8358, 291),
    "W360x39,0":  (353, 128, 6.5, 10.7, 50.2, 10331, 375),
    "W410x38,8":  (399, 140, 6.4, 8.8, 50.3, 12777, 404),
    "W410x46,1":  (403, 140, 7.0, 11.2, 59.2, 15690, 514),
    "W460x52,0":  (450, 152, 7.6, 10.8, 66.6, 21370, 634),
    "W530x66,0":  (525, 165, 8.9, 11.4, 83.6, 34971, 857),
}

# Cantoneiras laminadas de abas iguais: nome: (b, t) [mm]
CATALOGO_L = {
    'L 1"x1/8"': (25.4, 3.18), 'L 1.1/4"x1/8"': (31.75, 3.18),
    'L 1.1/2"x1/8"': (38.1, 3.18), 'L 1.1/2"x3/16"': (38.1, 4.76),
    'L 2"x3/16"': (50.8, 4.76), 'L 2"x1/4"': (50.8, 6.35),
    'L 2.1/2"x3/16"': (63.5, 4.76), 'L 2.1/2"x1/4"': (63.5, 6.35),
    'L 3"x1/4"': (76.2, 6.35), 'L 3"x5/16"': (76.2, 7.94),
    'L 4"x1/4"': (101.6, 6.35), 'L 4"x5/16"': (101.6, 7.94), 'L 4"x3/8"': (101.6, 9.53),
}

# Formados a frio: U: (bw, bf, t) | Ue: (bw, bf, D, t) [mm]
CATALOGO_U = {
    "U 75x40x2,00": (75, 40, 2.00), "U 100x40x2,00": (100, 40, 2.00),
    "U 100x50x2,65": (100, 50, 2.65), "U 127x50x2,65": (127, 50, 2.65),
    "U 150x50x2,65": (150, 50, 2.65), "U 200x50x3,00": (200, 50, 3.00),
}
CATALOGO_UE = {
    "Ue 75x40x15x2,00": (75, 40, 15, 2.00), "Ue 100x50x17x2,00": (100, 50, 17, 2.00),
    "Ue 100x50x17x2,65": (100, 50, 17, 2.65), "Ue 127x50x17x2,65": (127, 50, 17, 2.65),
    "Ue 150x60x20x2,00": (150, 60, 20, 2.00), "Ue 150x60x20x2,65": (150, 60, 20, 2.65),
    "Ue 200x75x25x2,65": (200, 75, 25, 2.65), "Ue 200x75x25x3,00": (200, 75, 25, 3.00),
    "Ue 250x85x25x3,00": (250, 85, 25, 3.00), "Ue 300x85x25x3,00": (300, 85, 25, 3.00),
}

FAMILIAS = {
    "W":   "Laminado — perfil W (catálogo Gerdau)",
    "I":   "Laminado/soldado — I por dimensões",
    "L":   "Laminado — cantoneira abas iguais",
    "U":   "Formado a frio — U simples",
    "Ue":  "Formado a frio — U enrijecido (Ue)",
    "GEN": "Genérica (A e I digitados)",
}
FAMILIA_FORMADO_FRIO = {"U", "Ue"}


def material_padrao_da_familia(fam: str) -> str:
    """Formado a frio → ASTM A36 ; laminado → ASTM A572 Gr50 (padrão do escritório)."""
    return "ASTM A36" if fam in FAMILIA_FORMADO_FRIO else "ASTM A572 Gr50"


# ---------- propriedades geométricas (mm) --------------------------------------------
def _props_retangulos(rects):
    """rects: lista de (b, h, xc, yc) — retângulos com lados paralelos aos eixos.
    Retorna A, xg, yg, Ix (eixo horizontal pelo CG), Iy (eixo vertical pelo CG)."""
    A = sum(b * h for b, h, _, _ in rects)
    xg = sum(b * h * xc for b, h, xc, _ in rects) / A
    yg = sum(b * h * yc for b, h, _, yc in rects) / A
    Ix = sum(b * h ** 3 / 12 + b * h * (yc - yg) ** 2 for b, h, _, yc in rects)
    Iy = sum(h * b ** 3 / 12 + b * h * (xc - xg) ** 2 for b, h, xc, _ in rects)
    return A, xg, yg, Ix, Iy


def props_U(bw, bf, t, D=0.0):
    """U / Ue formado a frio pela LINHA MÉDIA com cantos vivos (sem raio de dobra)."""
    h = bw - t                # altura da alma (linha média)
    b = bf - t / 2            # largura da mesa (linha média, da alma até a borda)
    c = max(D - t / 2, 0.0)   # enrijecedor (linha média)
    rects = [(t, h, 0.0, 0.0),                   # alma
             (b, t, b / 2, h / 2), (b, t, b / 2, -h / 2)]  # mesas
    if c > 0:
        rects += [(t, c, b, h / 2 - c / 2), (t, c, b, -h / 2 + c / 2)]
    return _props_retangulos(rects)


def props_L(b, t):
    """Cantoneira de abas iguais (sem raio de concordância). Ix = Iy (eixos paralelos às abas)."""
    rects = [(b, t, b / 2, t / 2), (t, b - t, t / 2, t + (b - t) / 2)]
    return _props_retangulos(rects)


def props_I(d, bf, tw, tf):
    rects = [(bf, tf, 0, d / 2 - tf / 2), (bf, tf, 0, -d / 2 + tf / 2),
             (tw, d - 2 * tf, 0, 0)]
    return _props_retangulos(rects)


FAMILIAS_COMPOSTAS = {"L", "U", "Ue"}     # admitem perfil duplo "costas com costas"


def _pol_area_centroide(p):
    """Área e centroide de polígono simples (fórmula de Gauss)."""
    A = cx = cy = 0.0
    for (x1, y1), (x2, y2) in zip(p, p[1:] + p[:1]):
        cr = x1 * y2 - x2 * y1
        A += cr
        cx += (x1 + x2) * cr
        cy += (y1 + y2) * cr
    A /= 2
    return A, cx / (6 * A), cy / (6 * A)


@dataclass
class Secao:
    """Seção transversal. Sistema da seção SEM rotação (β = 0):
       y = eixo vertical, no plano da estrutura (alma de W/I/U na vertical);
       z = eixo horizontal, perpendicular ao plano.
       Ix = inércia em torno do eixo horizontal (flexão no plano quando β = 0);
       Iy = inércia em torno do eixo vertical (flexão no plano quando β = 90°)."""
    nome: str
    familia: str = "W"                       # W | I | L | U | Ue | GEN
    material: str = "ASTM A572 Gr50"
    dims: dict = field(default_factory=dict) # mm ; para W: {"catalogo": nome} ; "gap" p/ composto
    eixo: str = "x"                          # (obsoleto — mantido p/ abrir arquivos antigos)
    qtd: int = 1                             # 1 = simples ; 2 = composto costas com costas
    cor: str = "#1f77b4"

    def composto(self):
        return int(self.qtd) == 2 and self.familia in FAMILIAS_COMPOSTAS

    def propriedades(self):
        """A [cm²], Ix, Iy [cm⁴] da seção simples (A1, Ix1, Iy1) e da seção usada (A, Ix, Iy)."""
        f, d = self.familia, self.dims
        e_costas = 0.0     # distância do CG do perfil simples à face das costas (mm)
        if f == "W":
            cat = CATALOGO_W[d.get("catalogo", "W200x15,0")]
            A, Ix, Iy = cat[4], cat[5], cat[6]
            obs = "valores de catálogo (conferir catálogo Gerdau)"
        elif f == "I":
            A, _, _, Ix, Iy = props_I(d["d"], d["bf"], d["tw"], d["tf"])
            A, Ix, Iy = A / 100, Ix / 1e4, Iy / 1e4
            obs = "retângulos, sem raios de concordância"
        elif f == "L":
            A, xg, _, Ix, Iy = props_L(d["b"], d["t"])
            e_costas = xg
            A, Ix, Iy = A / 100, Ix / 1e4, Iy / 1e4
            obs = "eixos geométricos (paralelos às abas), sem raio de concordância"
        elif f in ("U", "Ue"):
            A, xg, _, Ix, Iy = props_U(d["bw"], d["bf"], d["t"], d.get("D", 0.0) if f == "Ue" else 0.0)
            e_costas = xg + d["t"] / 2
            A, Ix, Iy = A / 100, Ix / 1e4, Iy / 1e4
            obs = "linha média, cantos vivos (sem raio de dobra)"
        else:  # GEN
            A, Ix, Iy = d.get("A", 10.0), d.get("I", 100.0), d.get("Iy", d.get("I", 100.0))
            obs = "valores digitados"
        out = {"A1": A, "Ix1": Ix, "Iy1": Iy, "A": A, "Ix": Ix, "Iy": Iy, "qtd": 1, "obs": obs}
        if self.composto():
            g = float(d.get("gap", 0.0))
            e = (g / 2 + e_costas) / 10          # cm
            out.update(A=2 * A, Ix=2 * Ix, Iy=2 * (Iy + A * e * e), qtd=2,
                       obs=obs + f"; composto costas c/ costas, afastamento {g:g} mm")
        return out

    def I_plano(self, beta=0.0):
        """Inércia (cm⁴) para flexão no plano com o perfil girado de β graus em torno do eixo da barra."""
        p = self.propriedades()
        b = math.radians(beta)
        return p["Ix"] * math.cos(b) ** 2 + p["Iy"] * math.sin(b) ** 2

    def A_m2(self):
        return self.propriedades()["A"] * 1e-4

    def I_m4(self, beta=0.0):
        return self.I_plano(beta) * 1e-8

    def descricao(self):
        p = self.propriedades()
        pre = "2x " if p["qtd"] == 2 and not self.nome.startswith("2x") else ""
        return (f"{pre}{self.nome} | {self.material} | A={p['A']:.2f} cm² | "
                f"Ix={p['Ix']:.1f} cm⁴ | Iy={p['Iy']:.1f} cm⁴")

    # ---------------- geometria para desenho ----------------
    def poligonos(self, beta=0.0):
        """Contorno(s) da seção em mm, com o CG na origem, girados de β (graus).
        Coordenadas (z, y): y no plano da estrutura (eixo local y da barra)."""
        f, d = self.familia, self.dims
        if f in ("W", "I"):
            if f == "W":
                dd, bf, tw, tf = CATALOGO_W[d.get("catalogo", "W200x15,0")][:4]
            else:
                dd, bf, tw, tf = d["d"], d["bf"], d["tw"], d["tf"]
            H, B, w = dd / 2, bf / 2, tw / 2
            p = [(-B, -H), (B, -H), (B, -H + tf), (w, -H + tf), (w, H - tf), (B, H - tf),
                 (B, H), (-B, H), (-B, H - tf), (-w, H - tf), (-w, -H + tf), (-B, -H + tf)]
            pols = [p]
        elif f in ("U", "Ue"):
            bw, bf, t = d["bw"], d["bf"], d["t"]
            H = bw / 2
            if f == "Ue" and d.get("D", 0) > t:
                D = d["D"]
                p = [(0, -H), (bf, -H), (bf, -H + D), (bf - t, -H + D), (bf - t, -H + t), (t, -H + t),
                     (t, H - t), (bf - t, H - t), (bf - t, H - D), (bf, H - D), (bf, H), (0, H)]
            else:
                p = [(0, -H), (bf, -H), (bf, -H + t), (t, -H + t), (t, H - t), (bf, H - t), (bf, H), (0, H)]
            pols = [p]
        elif f == "L":
            b, t = d["b"], d["t"]
            pols = [[(0, 0), (b, 0), (b, t), (t, t), (t, b), (0, b)]]
        else:
            A, I = d.get("A", 10.0) * 100, d.get("I", 100.0) * 1e4
            h = math.sqrt(12 * I / A)
            bb = A / h
            pols = [[(-bb / 2, -h / 2), (bb / 2, -h / 2), (bb / 2, h / 2), (-bb / 2, h / 2)]]
        if self.composto():
            g = float(d.get("gap", 0.0))
            pols = pols + [[(-z - g, y) for z, y in reversed(pols[0])]]
        # CG do conjunto
        At = cz = cy = 0.0
        for p in pols:
            a, zc, yc = _pol_area_centroide(p)
            At += a
            cz += a * zc
            cy += a * yc
        cz, cy = cz / At, cy / At
        cb, sb = math.cos(math.radians(beta)), math.sin(math.radians(beta))
        return [[((z - cz) * cb - (y - cy) * sb, (z - cz) * sb + (y - cy) * cb) for z, y in p] for p in pols]

    def linhas_elevacao(self, beta=0.0):
        """Para a vista extrudada: (y_min, y_max, [y das arestas internas visíveis]) em mm."""
        pols = self.poligonos(beta)
        ys = [y for p in pols for _, y in p]
        ymin, ymax = min(ys), max(ys)
        ints = set()
        for p in pols:
            for (z1, y1), (z2, y2) in zip(p, p[1:] + p[:1]):
                if abs(y1 - y2) < 1e-6 and abs(z1 - z2) > 1e-6:
                    y = round(y1, 3)
                    if ymin + 1e-3 < y < ymax - 1e-3:
                        ints.add(y)
        return ymin, ymax, sorted(ints)


# =============================================================================
# ENTIDADES DO MODELO
# =============================================================================
@dataclass
class No:
    id: int
    x: float
    y: float


@dataclass
class Barra:
    id: int
    ni: int
    nj: int
    secao: str
    rot_i: bool = False   # rótula (momento liberado) no nó inicial
    rot_j: bool = False   # rótula no nó final
    beta: float = 0.0     # giro do perfil em torno do eixo da barra (graus)
    pp: bool = True       # considerar peso próprio desta barra
    proj: dict = field(default_factory=dict)   # parâmetros de dimensionamento (flambagem, Lb, Cb, flecha)

    @property
    def tipo(self):
        if self.rot_i and self.rot_j:
            return "treliça (rotulada)"
        if self.rot_i or self.rot_j:
            return "rótula em uma extremidade"
        return "pórtico (rígida)"


TIPOS_GL = ("livre", "fixo", "mola")


@dataclass
class Apoio:
    no: int
    ux: str = "fixo"
    uy: str = "fixo"
    rz: str = "livre"
    kx: float = 0.0   # kN/m
    ky: float = 0.0   # kN/m
    kr: float = 0.0   # kN·m/rad

    def descricao(self):
        s = []
        for nome, t, k, un in (("ux", self.ux, self.kx, "kN/m"), ("uy", self.uy, self.ky, "kN/m"),
                               ("rz", self.rz, self.kr, "kN·m/rad")):
            s.append(f"{nome}:{t}" + (f"({k:g} {un})" if t == "mola" else ""))
        return " ".join(s)


@dataclass
class CasoCarga:
    nome: str
    peso_proprio: bool = False
    descricao: str = ""
    # --- classificação para combinações (NBR 8800 / NBR 8681) ---
    natureza: str = ""        # "permanente" | "variavel" | "ignorar"  ("" = ainda não classificado)
    categoria: str = ""       # categoria da Tabela 1 (permanentes) ou Tabela 2 (variáveis)
    grupo: str = ""           # grupo de exclusão mútua (ex.: "Vento"): só 1 caso do grupo por combinação
    gd: float = 1.0           # γg desfavorável (permanente)
    gf: float = 1.0           # γg favorável (permanente)
    gq: float = 1.0           # γq (variável)
    psi0: float = 0.0
    psi1: float = 0.0
    psi2: float = 0.0


@dataclass
class CargaNo:
    caso: str
    no: int
    fx: float = 0.0   # kN
    fy: float = 0.0   # kN
    mz: float = 0.0   # kN·m (anti-horário +)


DIRECOES = {"GX": "Global X", "GY": "Global Y", "LX": "Local x (axial)", "LY": "Local y (transversal)"}


@dataclass
class CargaBarra:
    caso: str
    barra: int
    tipo: str = "distribuida"   # "distribuida" (kN/m) | "pontual" (kN)
    direcao: str = "GY"
    valor: float = 0.0
    pos: float = 0.0            # pontual: distância ao nó i [m]
    projetada: bool = False     # distribuída em GX/GY por comprimento projetado


@dataclass
class Combinacao:
    nome: str
    fatores: dict = field(default_factory=dict)   # {nome_caso: fator}
    tipo: str = ""            # "ELU" | "ELS-R" | "ELS-F" | "ELS-QP" | "" (manual)
    descricao: str = ""
    auto: bool = False        # gerada automaticamente (substituída ao gerar de novo)


# =============================================================================
# MODELO
# =============================================================================
class Modelo:
    TOL = 1e-6

    def __init__(self):
        self.nos: dict[int, No] = {}
        self.barras: dict[int, Barra] = {}
        self.secoes: dict[str, Secao] = {}
        self.materiais: dict[str, Material] = {k: Material(**asdict(v)) for k, v in MATERIAIS_PADRAO.items()}
        self.apoios: dict[int, Apoio] = {}
        self.casos: dict[str, CasoCarga] = {}
        self.cargas_no: list[CargaNo] = []
        self.cargas_barra: list[CargaBarra] = []
        self.combinacoes: dict[str, Combinacao] = {}
        self.grade = {"dx": 1.0, "dy": 1.0, "snap": True}
        self.titulo = "Novo modelo"
        self.projeto = {}         # configuração global do dimensionamento (γa1, γa2, opções)

    # ---------------- nós ----------------
    def no_em(self, x, y, tol=None):
        tol = self.TOL if tol is None else tol
        for n in self.nos.values():
            if abs(n.x - x) <= tol and abs(n.y - y) <= tol:
                return n
        return None

    def add_no(self, x, y):
        n = self.no_em(x, y)
        if n:
            return n
        nid = max(self.nos, default=0) + 1
        n = No(nid, round(float(x), 9), round(float(y), 9))
        self.nos[nid] = n
        return n

    def remove_no(self, nid):
        for b in [b for b in self.barras.values() if nid in (b.ni, b.nj)]:
            self.remove_barra(b.id)
        self.nos.pop(nid, None)
        self.apoios.pop(nid, None)
        self.cargas_no = [c for c in self.cargas_no if c.no != nid]

    # ---------------- barras ----------------
    def add_barra(self, ni, nj, secao, rot_i=False, rot_j=False):
        if ni == nj:
            raise ValueError("Nó inicial e final iguais.")
        for b in self.barras.values():
            if {b.ni, b.nj} == {ni, nj}:
                raise ValueError(f"Já existe a barra {b.id} entre os nós {ni} e {nj}.")
        if secao not in self.secoes:
            raise ValueError("Cadastre/seleciona uma seção antes de lançar barras.")
        bid = max(self.barras, default=0) + 1
        b = Barra(bid, ni, nj, secao, rot_i, rot_j)
        self.barras[bid] = b
        return b

    def remove_barra(self, bid):
        self.barras.pop(bid, None)
        self.cargas_barra = [c for c in self.cargas_barra if c.barra != bid]

    # ---------------- divisão de barras ----------------
    def barra_sob_no(self, nid, tol=1e-6):
        """Barra cujo trecho INTERNO passa pelo nó (nó 'solto' sobre a barra). Retorna (barra, a) ou None."""
        n = self.nos[nid]
        for b in self.barras.values():
            if nid in (b.ni, b.nj):
                continue
            p1, p2 = self.nos[b.ni], self.nos[b.nj]
            L, c, s = self.geometria(b)
            a = (n.x - p1.x) * c + (n.y - p1.y) * s
            d = abs(-(n.x - p1.x) * s + (n.y - p1.y) * c)
            if d <= tol and tol < a < L - tol:
                return b, a
        return None

    def dividir_barra(self, bid, nid):
        """Divide a barra no nó nid (que deve estar sobre ela). A barra original fica com o trecho
        i→nó e uma nova barra recebe nó→j. Rótulas, seção, β e PP são preservados; cargas
        distribuídas são copiadas e cargas pontuais redistribuídas pela posição."""
        b = self.barras[bid]
        p1 = self.nos[b.ni]
        n = self.nos[nid]
        L, c, s = self.geometria(b)
        a = (n.x - p1.x) * c + (n.y - p1.y) * s
        novo_id = max(self.barras) + 1
        nb = Barra(novo_id, nid, b.nj, b.secao, False, b.rot_j, b.beta, b.pp, dict(b.proj))
        self.barras[novo_id] = nb
        b.nj, b.rot_j = nid, False
        novas = []
        for cg in self.cargas_barra:
            if cg.barra != bid:
                continue
            if cg.tipo == "distribuida":
                novas.append(CargaBarra(cg.caso, novo_id, cg.tipo, cg.direcao, cg.valor, 0.0, cg.projetada))
            elif cg.pos > a + 1e-9:
                cg.barra, cg.pos = novo_id, cg.pos - a
        self.cargas_barra += novas
        return nb

    def inserir_nos_nas_barras(self, tol=1e-6):
        """Divide todas as barras que têm nós soltos sobre elas. Retorna nº de divisões."""
        k = 0
        mudou = True
        while mudou:
            mudou = False
            for nid in list(self.nos):
                r = self.barra_sob_no(nid, tol)
                if r:
                    self.dividir_barra(r[0].id, nid)
                    k += 1
                    mudou = True
        return k

    def dividir_em_partes(self, bid, n):
        """Divide a barra em n trechos iguais (cria nós intermediários)."""
        b = self.barras[bid]
        p1, p2 = self.nos[b.ni], self.nos[b.nj]
        pts = [(p1.x + (p2.x - p1.x) * k / n, p1.y + (p2.y - p1.y) * k / n) for k in range(1, n)]
        atual = bid
        for x, y in pts:
            no = self.add_no(x, y)
            atual = self.dividir_barra(atual, no.id).id
        return n - 1

    def geometria(self, b: Barra):
        n1, n2 = self.nos[b.ni], self.nos[b.nj]
        dx, dy = n2.x - n1.x, n2.y - n1.y
        L = math.hypot(dx, dy)
        return L, dx / L, dy / L

    # ---------------- seções / casos ----------------
    def add_secao(self, s: Secao):
        self.secoes[s.nome] = s

    def remove_secao(self, nome):
        usados = [b.id for b in self.barras.values() if b.secao == nome]
        if usados:
            raise ValueError(f"Seção em uso nas barras {usados}.")
        self.secoes.pop(nome, None)

    def add_caso(self, nome, peso_proprio=False, descricao=""):
        if not nome:
            raise ValueError("Nome do caso vazio.")
        self.casos[nome] = CasoCarga(nome, peso_proprio, descricao)

    def remove_caso(self, nome):
        self.casos.pop(nome, None)
        self.cargas_no = [c for c in self.cargas_no if c.caso != nome]
        self.cargas_barra = [c for c in self.cargas_barra if c.caso != nome]
        for cb in self.combinacoes.values():
            cb.fatores.pop(nome, None)

    def renomear_caso(self, antigo, novo):
        if novo in self.casos:
            raise ValueError("Já existe um caso com esse nome.")
        c = self.casos.pop(antigo)
        c.nome = novo
        self.casos[novo] = c
        for cg in self.cargas_no + self.cargas_barra:
            if cg.caso == antigo:
                cg.caso = novo
        for cb in self.combinacoes.values():
            if antigo in cb.fatores:
                cb.fatores[novo] = cb.fatores.pop(antigo)

    # ---------------- validação ----------------
    def validar(self):
        erros = []
        if not self.barras:
            erros.append("Não há barras no modelo.")
        if not self.apoios:
            erros.append("Não há apoios no modelo.")
        usados = {b.ni for b in self.barras.values()} | {b.nj for b in self.barras.values()}
        soltos = [n for n in self.nos if n not in usados]
        avisos = [f"Nós sem barra (ignorados na análise): {soltos}"] if soltos else []
        for cg in self.cargas_no:
            if cg.no not in usados:
                sob = self.barra_sob_no(cg.no) if cg.no in self.nos else None
                erros.append(f"Carga do caso '{cg.caso}' no nó {cg.no}, que não está ligado a nenhuma barra"
                             + (f" (está sobre a barra {sob[0].id} sem dividi-la)." if sob else "."))
        for a in self.apoios:
            if a not in usados:
                erros.append(f"Apoio no nó {a}, que não está ligado a nenhuma barra.")
        for b in self.barras.values():
            if b.secao not in self.secoes:
                erros.append(f"Barra {b.id}: seção '{b.secao}' inexistente.")
            if self.geometria(b)[0] < 1e-6:
                erros.append(f"Barra {b.id}: comprimento nulo.")
        for c in self.cargas_barra:
            if c.barra in self.barras and c.tipo == "pontual":
                L = self.geometria(self.barras[c.barra])[0]
                if not (-1e-9 <= c.pos <= L + 1e-9):
                    erros.append(f"Carga pontual na barra {c.barra}: posição {c.pos} fora de [0, {L:.3f}] m.")
        return erros, avisos

    # ---------------- persistência ----------------
    def to_dict(self):
        return {
            "versao": 1, "titulo": self.titulo, "grade": self.grade, "projeto": self.projeto,
            "materiais": [asdict(m) for m in self.materiais.values()],
            "secoes": [asdict(s) for s in self.secoes.values()],
            "nos": [asdict(n) for n in self.nos.values()],
            "barras": [asdict(b) for b in self.barras.values()],
            "apoios": [asdict(a) for a in self.apoios.values()],
            "casos": [asdict(c) for c in self.casos.values()],
            "cargas_no": [asdict(c) for c in self.cargas_no],
            "cargas_barra": [asdict(c) for c in self.cargas_barra],
            "combinacoes": [asdict(c) for c in self.combinacoes.values()],
        }

    @classmethod
    def from_dict(cls, d):
        m = cls()
        m.titulo = d.get("titulo", "Modelo")
        m.grade.update(d.get("grade", {}))
        m.projeto = dict(d.get("projeto", {}))
        for x in d.get("materiais", []):
            m.materiais[x["nome"]] = Material(**x)
        for x in d.get("secoes", []):
            m.secoes[x["nome"]] = Secao(**x)
        for x in d.get("nos", []):
            m.nos[x["id"]] = No(**x)
        for x in d.get("barras", []):
            m.barras[x["id"]] = Barra(**x)
        for x in d.get("apoios", []):
            m.apoios[x["no"]] = Apoio(**x)
        for x in d.get("casos", []):
            m.casos[x["nome"]] = CasoCarga(**x)
        m.cargas_no = [CargaNo(**x) for x in d.get("cargas_no", [])]
        m.cargas_barra = [CargaBarra(**x) for x in d.get("cargas_barra", [])]
        for x in d.get("combinacoes", []):
            m.combinacoes[x["nome"]] = Combinacao(**x)
        # compatibilidade: arquivos antigos com seção "eixo y" → barra girada 90°
        for b in m.barras.values():
            s = m.secoes.get(b.secao)
            if s and s.eixo == "y" and not b.beta:
                b.beta = 90.0
        for s in m.secoes.values():
            s.eixo = "x"
        return m

    def salvar(self, caminho):
        with open(caminho, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, ensure_ascii=False, indent=1)

    @classmethod
    def abrir(cls, caminho):
        with open(caminho, encoding="utf-8") as f:
            return cls.from_dict(json.load(f))


# =============================================================================
# RESULTADOS
# =============================================================================
class Resultado:
    """Resultado de um caso ou combinação (superposição linear)."""

    def __init__(self, nome):
        self.nome = nome
        self.desloc = {}      # no -> [ux, uy, rz] (m, rad) global
        self.ext = {}         # barra -> ([ux,uy,rz]_i, [ux,uy,rz]_j) global, nas pontas da barra
        self.flocal = {}      # barra -> [Fx_i, Fy_i, M_i, Fx_j, Fy_j, M_j] (forças dos nós SOBRE a barra, local)
        self.cargas = {}      # barra -> {"wx": , "wy": , "pont": [(a, px, py)]}  (local)
        self.reacoes = {}     # no -> [Rx, Ry, Mz]
        self.soma_cargas = [0.0, 0.0]   # ΣFx, ΣFy de todas as cargas aplicadas (conferência de equilíbrio)

    @staticmethod
    def combinar(nome, lista):
        """lista: [(Resultado, fator)]"""
        r = Resultado(nome)
        for res, f in lista:
            for k, v in res.desloc.items():
                acc = r.desloc.setdefault(k, [0.0, 0.0, 0.0])
                for i in range(3):
                    acc[i] += f * v[i]
            for k, (di, dj) in res.ext.items():
                ai, aj = r.ext.setdefault(k, ([0.0] * 3, [0.0] * 3))
                for i in range(3):
                    ai[i] += f * di[i]
                    aj[i] += f * dj[i]
            for k, v in res.flocal.items():
                acc = r.flocal.setdefault(k, [0.0] * 6)
                for i in range(6):
                    acc[i] += f * v[i]
            for k, c in res.cargas.items():
                acc = r.cargas.setdefault(k, {"wx": 0.0, "wy": 0.0, "pont": []})
                acc["wx"] += f * c["wx"]
                acc["wy"] += f * c["wy"]
                acc["pont"] += [(a, f * px, f * py) for a, px, py in c["pont"]]
            for k, v in res.reacoes.items():
                acc = r.reacoes.setdefault(k, [0.0, 0.0, 0.0])
                for i in range(3):
                    acc[i] += f * v[i]
            r.soma_cargas = [r.soma_cargas[i] + f * res.soma_cargas[i] for i in range(2)]
        return r


class ErroAnalise(Exception):
    pass


# =============================================================================
# ANÁLISE COM OPENSEES
# =============================================================================
def _cargas_locais(modelo: Modelo, caso: str):
    """Converte cargas de barra do caso (e peso próprio) para o sistema local."""
    out = {bid: {"wx": 0.0, "wy": 0.0, "pont": []} for bid in modelo.barras}
    cc = modelo.casos.get(caso)
    if cc and cc.peso_proprio:
        for b in modelo.barras.values():
            if not b.pp:
                continue
            s = modelo.secoes[b.secao]
            gam = modelo.materiais.get(s.material, MATERIAIS_PADRAO["ASTM A36"]).gamma
            L, c, sn = modelo.geometria(b)
            gy = -gam * s.A_m2()                 # kN/m ao longo da barra, global −Y
            out[b.id]["wx"] += gy * sn
            out[b.id]["wy"] += gy * c
    for cg in modelo.cargas_barra:
        if cg.caso != caso or cg.barra not in modelo.barras:
            continue
        L, c, sn = modelo.geometria(modelo.barras[cg.barra])
        v = cg.valor
        if cg.direcao == "GX":
            gx, gy = v, 0.0
            if cg.tipo == "distribuida" and cg.projetada:
                gx *= abs(sn)
        elif cg.direcao == "GY":
            gx, gy = 0.0, v
            if cg.tipo == "distribuida" and cg.projetada:
                gy *= abs(c)
        else:
            gx = gy = None
        if gx is not None:
            lx, ly = gx * c + gy * sn, -gx * sn + gy * c
        else:
            lx, ly = (v, 0.0) if cg.direcao == "LX" else (0.0, v)
        if cg.tipo == "distribuida":
            out[cg.barra]["wx"] += lx
            out[cg.barra]["wy"] += ly
        else:
            out[cg.barra]["pont"].append((min(max(cg.pos, 0.0), L), lx, ly))
    return out


def soma_cargas_aplicadas(modelo, caso, cargas):
    """ΣFx, ΣFy (globais) das cargas nodais e de barra (inclusive peso próprio) de um caso."""
    sx = sy = 0.0
    for cg in modelo.cargas_no:
        if cg.caso == caso:
            sx += cg.fx
            sy += cg.fy
    for bid, c in cargas.items():
        L, co, si = modelo.geometria(modelo.barras[bid])
        lx = c["wx"] * L + sum(px for _, px, _ in c["pont"])
        ly = c["wy"] * L + sum(py for _, _, py in c["pont"])
        sx += lx * co - ly * si
        sy += lx * si + ly * co
    return [sx, sy]


def analisar_caso(modelo: Modelo, caso: str) -> Resultado:
    import openseespy.opensees as ops

    erros, _ = modelo.validar()
    if erros:
        raise ErroAnalise("\n".join(erros))

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    usados = sorted({b.ni for b in modelo.barras.values()} | {b.nj for b in modelo.barras.values()})
    for nid in usados:
        n = modelo.nos[nid]
        ops.node(nid, n.x, n.y)

    OFF_AUX, OFF_TERRA = 1_000_000, 2_000_000
    aux = {}                          # (bid, 'i'/'j') -> tag do nó auxiliar
    rigidas = {nid: 0 for nid in usados}
    ops.geomTransf("Linear", 1)
    for b in modelo.barras.values():
        ends = []
        for lado, nid, rot in (("i", b.ni, b.rot_i), ("j", b.nj, b.rot_j)):
            if rot:
                tag = OFF_AUX + 2 * b.id + (0 if lado == "i" else 1)
                n = modelo.nos[nid]
                ops.node(tag, n.x, n.y)
                ops.equalDOF(nid, tag, 1, 2)
                aux[(b.id, lado)] = tag
                ends.append(tag)
            else:
                rigidas[nid] += 1
                ends.append(nid)
        s = modelo.secoes[b.secao]
        E = modelo.materiais.get(s.material, MATERIAIS_PADRAO["ASTM A36"]).E * 1e3   # kPa = kN/m²
        ops.element("elasticBeamColumn", b.id, ends[0], ends[1], s.A_m2(), E, s.I_m4(b.beta), 1)

    # apoios
    n_fix = {}   # no -> nº de GL fixados
    molas = {}   # no -> lista de tags de nó terra com sua direção
    mat_tag = 1
    for nid, a in modelo.apoios.items():
        if nid not in rigidas:
            continue
        fx = [1 if t == "fixo" else 0 for t in (a.ux, a.uy, a.rz)]
        if a.rz != "fixo" and a.rz != "mola" and rigidas[nid] == 0:
            fx[2] = 1   # rotação do nó sem barra rígida é indeterminada → bloqueia
        if any(fx):
            ops.fix(nid, *fx)
            n_fix[nid] = sum(fx)
        for gl, (t, k) in enumerate(((a.ux, a.kx), (a.uy, a.ky), (a.rz, a.kr)), start=1):
            if t == "mola":
                if k <= 0:
                    raise ErroAnalise(f"Apoio no nó {nid}: rigidez da mola deve ser > 0.")
                tg = OFF_TERRA + 3 * nid + gl
                n = modelo.nos[nid]
                ops.node(tg, n.x, n.y)
                ops.fix(tg, 1, 1, 1)
                ops.uniaxialMaterial("Elastic", mat_tag, k)
                ops.element("zeroLength", tg, tg, nid, "-mat", mat_tag, "-dir", gl)
                mat_tag += 1
                molas.setdefault(nid, []).append(tg)
    # nós sem nenhuma ligação rígida e sem apoio → bloqueia rotação (indeterminada)
    for nid, nr in rigidas.items():
        if nr == 0 and nid not in modelo.apoios:
            ops.fix(nid, 0, 0, 1)
            n_fix[nid] = 1

    # cargas
    cargas = _cargas_locais(modelo, caso)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for cg in modelo.cargas_no:
        if cg.caso == caso and cg.no in rigidas:
            ops.load(cg.no, cg.fx, cg.fy, cg.mz)
    for bid, c in cargas.items():
        L = modelo.geometria(modelo.barras[bid])[0]
        if abs(c["wx"]) > 0 or abs(c["wy"]) > 0:
            ops.eleLoad("-ele", bid, "-type", "-beamUniform", c["wy"], c["wx"])
        for a, px, py in c["pont"]:
            ops.eleLoad("-ele", bid, "-type", "-beamPoint", py, a / L, px)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    n_gl = 3 * len(usados) + len(aux) - sum(n_fix.values())   # equações livres
    if n_gl <= 0:
        ops.wipe()
        return _resultado_totalmente_fixo(modelo, caso, usados, _cargas_locais(modelo, caso))
    verificar = n_gl <= 2400          # checagem de singularidade pela matriz cheia
    ops.system("FullGeneral" if verificar else "UmfPack")
    ops.test("NormDispIncr", 1e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    if ok == 0 and verificar:
        _checa_hipostatica(ops, modelo, usados, aux)
    if ok != 0:
        ops.wipe()
        raise ErroAnalise("A análise falhou: estrutura hipostática ou instável (matriz singular). "
                          "Verifique apoios, rótulas e vinculações.")

    ops.reactions()
    r = Resultado(caso)
    for nid in usados:
        r.desloc[nid] = list(ops.nodeDisp(nid))
    lim = 1e3
    for v in r.desloc.values():
        if any(abs(x) > lim or math.isnan(x) for x in v):
            ops.wipe()
            raise ErroAnalise("Deslocamentos absurdos: estrutura provavelmente hipostática (mecanismo).")
    for b in modelo.barras.values():
        ti = aux.get((b.id, "i"), b.ni)
        tj = aux.get((b.id, "j"), b.nj)
        r.ext[b.id] = (list(ops.nodeDisp(ti)), list(ops.nodeDisp(tj)))
        r.flocal[b.id] = list(ops.eleResponse(b.id, "localForce"))
        r.cargas[b.id] = cargas[b.id]
    # Reação = Σ forças das barras que chegam ao nó (inclusive pelas extremidades rotuladas,
    # ligadas ao nó por equalDOF) − cargas nodais aplicadas no próprio nó.
    # (ops.nodeReaction não enxerga as forças transmitidas pelo equalDOF das rótulas.)
    forcas_no = {}
    for b in modelo.barras.values():
        gf = ops.eleResponse(b.id, "globalForce")
        for nid, f in ((b.ni, gf[0:3]), (b.nj, gf[3:6])):
            acc = forcas_no.setdefault(nid, [0.0, 0.0, 0.0])
            for k in range(3):
                acc[k] += f[k]
    for nid, a in modelo.apoios.items():
        if nid not in rigidas:
            continue
        R = list(forcas_no.get(nid, [0.0, 0.0, 0.0]))
        for cg in modelo.cargas_no:
            if cg.caso == caso and cg.no == nid:
                R[0] -= cg.fx
                R[1] -= cg.fy
                R[2] -= cg.mz
        if a.rz == "livre":
            R[2] = 0.0
        r.reacoes[nid] = [0.0 if abs(x) < 1e-9 else x for x in R]
    r.soma_cargas = soma_cargas_aplicadas(modelo, caso, cargas)
    ops.wipe()
    return r


def _resultado_totalmente_fixo(modelo, caso, usados, cargas):
    """Nenhum GL livre (ex.: viga biengastada sem nós internos): deslocamentos nulos e
    esforços de engastamento perfeito, calculados analiticamente."""
    r = Resultado(caso)
    for nid in usados:
        r.desloc[nid] = [0.0, 0.0, 0.0]
    soma = {nid: [0.0, 0.0, 0.0] for nid in usados}
    for b in modelo.barras.values():
        L, c, s = modelo.geometria(b)
        cg = cargas[b.id]
        wx, wy = cg["wx"], cg["wy"]
        f = [-wx * L / 2, -wy * L / 2, -wy * L * L / 12, -wx * L / 2, -wy * L / 2, wy * L * L / 12]
        for a, px, py in cg["pont"]:
            bb = L - a
            f[0] += -px * bb / L
            f[3] += -px * a / L
            f[1] += -py * bb * bb * (3 * a + bb) / L ** 3
            f[4] += -py * a * a * (a + 3 * bb) / L ** 3
            f[2] += -py * a * bb * bb / L ** 2
            f[5] += py * a * a * bb / L ** 2
        r.flocal[b.id] = f
        r.cargas[b.id] = cg
        r.ext[b.id] = ([0.0] * 3, [0.0] * 3)
        for nid, (fx, fy, mz) in ((b.ni, f[0:3]), (b.nj, f[3:6])):
            soma[nid][0] += fx * c - fy * s
            soma[nid][1] += fx * s + fy * c
            soma[nid][2] += mz
    for cgn in modelo.cargas_no:
        if cgn.caso == caso and cgn.no in soma:
            soma[cgn.no][0] -= cgn.fx
            soma[cgn.no][1] -= cgn.fy
            soma[cgn.no][2] -= cgn.mz
    for nid in modelo.apoios:
        if nid in soma:
            r.reacoes[nid] = soma[nid]
    return r


def _checa_hipostatica(ops, modelo, usados, aux):
    """Menor autovalor da rigidez escalada pela diagonal ≈ 0 → mecanismo.
    Indica os nós/graus de liberdade que participam do modo de mecanismo."""
    import numpy as np
    K = np.array(ops.printA("-ret"), dtype=float)
    n = int(round(math.sqrt(K.size)))
    if n == 0:
        return
    K = K.reshape(n, n)
    d = np.abs(np.diag(K))
    if np.any(d <= 0):
        lam_min, modo = 0.0, (d <= 0).astype(float)
    else:
        Dm = 1.0 / np.sqrt(d)
        Ks = (K * Dm[:, None]) * Dm[None, :]
        Ks = 0.5 * (Ks + Ks.T)
        w, v = np.linalg.eigh(Ks)
        lam_min, modo = w[0], np.abs(v[:, 0])
    if lam_min > 1e-9:
        return
    nomes_gl = ("ux", "uy", "rz")
    tags = {nid: f"nó {nid}" for nid in usados}
    for (bid, lado), tg in aux.items():
        tags[tg] = f"barra {bid}, extremidade {lado} (rótula)"
    envolvidos = []
    lim = 0.2 * modo.max()
    for tg, rot in tags.items():
        eqs = ops.nodeDOFs(tg)
        for k, eq in enumerate(eqs):
            if tg >= 1_000_000 and k < 2:
                continue            # translações do nó de rótula = as do nó principal
            if 0 <= eq < n and modo[eq] >= lim:
                envolvidos.append(f"{rot}: {nomes_gl[k]}")
    ops.wipe()
    raise ErroAnalise("Estrutura hipostática (mecanismo) — matriz de rigidez singular.\n"
                      "Graus de liberdade envolvidos no mecanismo:\n  " +
                      "\n  ".join(envolvidos[:12]) +
                      ("\n  ..." if len(envolvidos) > 12 else "") +
                      "\nVerifique apoios, rótulas (ex.: nó com todas as barras rotuladas) e vinculações.")


def analisar(modelo: Modelo):
    """Analisa todos os casos e combinações. Retorna dict nome -> Resultado."""
    if not modelo.casos:
        raise ErroAnalise("Defina ao menos um caso de carga.")
    res = {}
    for nome in modelo.casos:
        res[nome] = analisar_caso(modelo, nome)
    for nome, cb in modelo.combinacoes.items():
        lista = [(res[c], f) for c, f in cb.fatores.items() if c in res and f != 0]
        res["COMB: " + nome] = Resultado.combinar(nome, lista)
    return res


# =============================================================================
# PÓS-PROCESSAMENTO (esforços ao longo da barra e deformada exata)
# =============================================================================
def esforcos_barra(modelo: Modelo, res: Resultado, bid: int, n: int = 40, xs=None):
    """Retorna listas xs, N, V, M ao longo da barra (convenção do cabeçalho)."""
    b = modelo.barras[bid]
    L = modelo.geometria(b)[0]
    Fx_i, Fy_i, M_i = res.flocal[bid][:3]
    c = res.cargas.get(bid, {"wx": 0.0, "wy": 0.0, "pont": []})
    wx, wy, pont = c["wx"], c["wy"], c["pont"]
    if xs is not None:
        xs = sorted(min(max(float(x), 0.0), L) for x in xs)
    else:
        pts = {L * k / n for k in range(n + 1)}
        eps = 1e-7 * L
        for a, _, _ in pont:
            pts |= {max(a - eps, 0.0), min(a + eps, L), a}
        xs = sorted(pts)
    N, V, M = [], [], []
    for x in xs:
        sp_x = sum(px for a, px, _ in pont if a < x)
        sp_y = sum(py for a, _, py in pont if a < x)
        sm = sum((x - a) * py for a, _, py in pont if a < x)
        N.append(-Fx_i - wx * x - sp_x)
        V.append(Fy_i + wy * x + sp_y)
        M.append(-M_i + Fy_i * x + wy * x * x / 2 + sm)
    return xs, N, V, M


def extremos(xs, vals):
    imax = max(range(len(vals)), key=lambda k: vals[k])
    imin = min(range(len(vals)), key=lambda k: vals[k])
    return (xs[imax], vals[imax]), (xs[imin], vals[imin])


def deformada_barra(modelo: Modelo, res: Resultado, bid: int, n: int = 30, xs=None):
    """Deslocamentos locais exatos (Euler-Bernoulli): Hermite com valores nas pontas +
    solução particular biengastada das cargas no vão. Retorna xs, u(x), v(x) locais [m]."""
    b = modelo.barras[bid]
    L, c, s = modelo.geometria(b)
    sec = modelo.secoes[b.secao]
    E = modelo.materiais.get(sec.material, MATERIAIS_PADRAO["ASTM A36"]).E * 1e3
    EA, EI = E * sec.A_m2(), E * sec.I_m4(b.beta)
    di, dj = res.ext[bid]
    u1 = di[0] * c + di[1] * s
    v1 = -di[0] * s + di[1] * c
    u2 = dj[0] * c + dj[1] * s
    v2 = -dj[0] * s + dj[1] * c
    t1, t2 = di[2], dj[2]
    cg = res.cargas.get(bid, {"wx": 0.0, "wy": 0.0, "pont": []})
    xs = [L * k / n for k in range(n + 1)] if xs is None else [min(max(float(x), 0.0), L) for x in xs]
    U, Vv = [], []
    for x in xs:
        xi = x / L
        H1 = 1 - 3 * xi ** 2 + 2 * xi ** 3
        H2 = L * (xi - 2 * xi ** 2 + xi ** 3)
        H3 = 3 * xi ** 2 - 2 * xi ** 3
        H4 = L * (-xi ** 2 + xi ** 3)
        v = H1 * v1 + H2 * t1 + H3 * v2 + H4 * t2
        u = u1 + (u2 - u1) * xi
        # particulares (biengastada / bi-fixa axial)
        v += cg["wy"] * x ** 2 * (L - x) ** 2 / (24 * EI)
        u += cg["wx"] * x * (L - x) / (2 * EA)
        for a, px, py in cg["pont"]:
            bb = L - a
            if x <= a:
                v += py * bb ** 2 * x ** 2 * (3 * a * L - (3 * a + bb) * x) / (6 * EI * L ** 3)
                u += px * x * bb / (L * EA)
            else:
                xr = L - x
                v += py * a ** 2 * xr ** 2 * (3 * bb * L - (3 * bb + a) * xr) / (6 * EI * L ** 3)
                u += px * a * (L - x) / (L * EA)
        U.append(u)
        Vv.append(v)
    return xs, U, Vv


def desloc_global_barra(modelo, res, bid, xs):
    """Deslocamentos GLOBAIS (dx, dy) em m nos pontos xs (m, a partir do nó i)."""
    b = modelo.barras[bid]
    _, c, s = modelo.geometria(b)
    xs, U, V = deformada_barra(modelo, res, bid, xs=xs)
    return xs, [c * u - s * v for u, v in zip(U, V)], [s * u + c * v for u, v in zip(U, V)]


def desloc_max_absoluto(modelo, res, bid, n=80):
    """Maior |δ| (vetor global) ao longo da barra: (x, dx, dy) em m."""
    L = modelo.geometria(modelo.barras[bid])[0]
    xs, DX, DY = desloc_global_barra(modelo, res, bid, [L * k / n for k in range(n + 1)])
    k = max(range(len(xs)), key=lambda i: math.hypot(DX[i], DY[i]))
    # refina em torno do máximo
    a, b_ = xs[max(k - 1, 0)], xs[min(k + 1, n)]
    xs2, DX2, DY2 = desloc_global_barra(modelo, res, bid, [a + (b_ - a) * j / 40 for j in range(41)])
    j = max(range(len(xs2)), key=lambda i: math.hypot(DX2[i], DY2[i]))
    return xs2[j], DX2[j], DY2[j]


def deslocamento_max_barra(modelo, res, bid):
    """Maior deslocamento transversal RELATIVO à corda (flecha) e absoluto em mm."""
    xs, U, V = deformada_barra(modelo, res, bid, 60)
    L = xs[-1]
    rel = [V[k] - (V[0] + (V[-1] - V[0]) * xs[k] / L) for k in range(len(xs))]
    k = max(range(len(rel)), key=lambda i: abs(rel[i]))
    return xs[k], rel[k] * 1e3


# =============================================================================
# MEMORIAL (.txt)
# =============================================================================
def memorial_txt(modelo: Modelo, resultados: dict | None) -> str:
    L = []
    w = L.append
    w("=" * 78)
    w(f"MEMORIAL — ANÁLISE ESTRUTURAL 2D (OpenSees)  |  {modelo.titulo}")
    w("=" * 78)
    w("Unidades: kN, m, kN·m. Análise elástica linear de 1ª ordem (elasticBeamColumn).")
    w("Convenções: N>0 tração | M>0 traciona o lado −y local (i→j) | Reações nos eixos globais.")
    w("")
    w("MATERIAIS")
    for m in modelo.materiais.values():
        w(f"  {m.nome:<16} E={m.E:.0f} MPa  G={m.G:.0f} MPa  fy={m.fy:.0f} MPa  fu={m.fu:.0f} MPa  γ={m.gamma} kN/m³")
    w("")
    w("SEÇÕES")
    for s in modelo.secoes.values():
        p = s.propriedades()
        w(f"  {s.nome:<22} [{FAMILIAS.get(s.familia, s.familia)}]")
        w(f"      {s.descricao()}")
        w(f"      Perfil simples: A={p['A1']:.2f} cm²  Ix={p['Ix1']:.1f} cm⁴  Iy={p['Iy1']:.1f} cm⁴  — {p['obs']}")
    w("")
    w("NÓS")
    w(f"  {'Nó':>4} {'X (m)':>10} {'Y (m)':>10}")
    for n in sorted(modelo.nos.values(), key=lambda n: n.id):
        w(f"  {n.id:>4} {n.x:>10.3f} {n.y:>10.3f}")
    w("")
    w("BARRAS")
    w(f"  {'Barra':>5} {'Ni':>4} {'Nj':>4} {'L (m)':>8} {'β(°)':>6} {'PP':>4} {'I plano cm⁴':>12}  {'Tipo':<26} Seção")
    for b in sorted(modelo.barras.values(), key=lambda b: b.id):
        s = modelo.secoes[b.secao]
        w(f"  {b.id:>5} {b.ni:>4} {b.nj:>4} {modelo.geometria(b)[0]:>8.3f} {b.beta:>6.1f} {'sim' if b.pp else 'não':>4} "
          f"{s.I_plano(b.beta):>12.1f}  {b.tipo:<26} {b.secao}")
    w("")
    w("APOIOS")
    for a in modelo.apoios.values():
        w(f"  Nó {a.no}: {a.descricao()}")
    w("")
    w("CASOS DE CARGA")
    for c in modelo.casos.values():
        w(f"  [{c.nome}]" + ("  + peso próprio automático" if c.peso_proprio else "") +
          (f"  — {c.descricao}" if c.descricao else ""))
        if c.natureza == "permanente":
            w(f"      Permanente: {c.categoria}  (γg = {c.gd:g} desf. / {c.gf:g} fav.)")
        elif c.natureza == "variavel":
            w(f"      Variável: {c.categoria}  (γq = {c.gq:g}; ψ0 = {c.psi0:g}, ψ1 = {c.psi1:g}, ψ2 = {c.psi2:g})"
              + (f"  grupo exclusivo: {c.grupo}" if c.grupo else ""))
        for cg in modelo.cargas_no:
            if cg.caso == c.nome:
                w(f"      Nó {cg.no}: Fx={cg.fx:g} kN  Fy={cg.fy:g} kN  Mz={cg.mz:g} kN·m")
        for cg in modelo.cargas_barra:
            if cg.caso == c.nome:
                if cg.tipo == "distribuida":
                    w(f"      Barra {cg.barra}: distribuída {cg.valor:g} kN/m  dir. {DIRECOES[cg.direcao]}"
                      + (" (proj.)" if cg.projetada else ""))
                else:
                    w(f"      Barra {cg.barra}: pontual {cg.valor:g} kN  dir. {DIRECOES[cg.direcao]}  a={cg.pos:g} m do nó i")
    w("")
    if modelo.combinacoes:
        w("COMBINAÇÕES")
        for cb in modelo.combinacoes.values():
            w(f"  {cb.nome:<8} " + (f"[{cb.tipo}] " if cb.tipo else "[manual] ") +
              " + ".join(f"{f:g}·[{c}]" for c, f in cb.fatores.items() if f))
        w("")
    if not resultados:
        return "\n".join(L)
    for nome, r in resultados.items():
        w("-" * 78)
        w(f"RESULTADOS — {nome}")
        w("-" * 78)
        w("  Reações de apoio")
        w(f"  {'Nó':>4} {'Rx (kN)':>11} {'Ry (kN)':>11} {'Mz (kN·m)':>11}")
        sx = sy = 0.0
        for nid, R in sorted(r.reacoes.items()):
            w(f"  {nid:>4} {R[0]:>11.3f} {R[1]:>11.3f} {R[2]:>11.3f}")
            sx += R[0]
            sy += R[1]
        w(f"  {'Σ':>4} {sx:>11.3f} {sy:>11.3f}")
        w("  Deslocamentos nodais")
        w(f"  {'Nó':>4} {'ux (mm)':>10} {'uy (mm)':>10} {'rz (rad)':>12}")
        for nid, d in sorted(r.desloc.items()):
            w(f"  {nid:>4} {d[0]*1e3:>10.3f} {d[1]*1e3:>10.3f} {d[2]:>12.3e}")
        w("  Esforços nas barras (extremidades e máximos)")
        w(f"  {'Barra':>5} {'N_i':>9} {'V_i':>9} {'M_i':>9} {'N_j':>9} {'V_j':>9} {'M_j':>9} {'M_max':>9} {'@x':>6} {'M_min':>9} {'@x':>6} {'flecha mm':>9}")
        for bid in sorted(modelo.barras):
            xs, N, V, M = esforcos_barra(modelo, r, bid)
            (xM, Mmax), (xm, Mmin) = extremos(xs, M)
            _, fl = deslocamento_max_barra(modelo, r, bid)
            w(f"  {bid:>5} {N[0]:>9.2f} {V[0]:>9.2f} {M[0]:>9.2f} {N[-1]:>9.2f} {V[-1]:>9.2f} {M[-1]:>9.2f} "
              f"{Mmax:>9.2f} {xM:>6.2f} {Mmin:>9.2f} {xm:>6.2f} {fl:>9.2f}")
        w("")
    return "\n".join(L)
