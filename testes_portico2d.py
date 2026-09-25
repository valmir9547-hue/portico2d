# -*- coding: utf-8 -*-
"""Testes de regressão do motor portico2d_core contra soluções analíticas clássicas."""
import math
from portico2d_core import (Modelo, Secao, Apoio, CargaNo, CargaBarra, Combinacao,
                            analisar, esforcos_barra, extremos, deformada_barra,
                            props_U, props_L, memorial_txt)

OK = []


def check(nome, calc, ref, tol=1e-3):
    err = abs(calc - ref) / max(abs(ref), 1e-9)
    st = "OK " if err <= tol else "FALHOU"
    OK.append(err <= tol)
    print(f"  [{st}] {nome:<48} calc={calc:>12.5f}  ref={ref:>12.5f}")


def secao_gen(m, A=50.0, I=5000.0):
    m.add_secao(Secao("GEN", "GEN", "ASTM A36", {"A": A, "I": I}))


def viga_biapoiada():
    print("1) Viga biapoiada L=6 m, q=10 kN/m + P=20 kN no meio")
    m = Modelo(); secao_gen(m)
    n1, n2 = m.add_no(0, 0), m.add_no(6, 0)
    b = m.add_barra(n1.id, n2.id, "GEN")
    m.apoios[1] = Apoio(1, "fixo", "fixo", "livre")
    m.apoios[2] = Apoio(2, "livre", "fixo", "livre")
    m.add_caso("G")
    m.cargas_barra += [CargaBarra("G", b.id, "distribuida", "GY", -10.0),
                       CargaBarra("G", b.id, "pontual", "GY", -20.0, 3.0)]
    r = analisar(m)["G"]
    xs, N, V, M = esforcos_barra(m, r, b.id)
    (xM, Mmax), _ = extremos(xs, M)
    check("Ry1 = qL/2 + P/2", r.reacoes[1][1], 40.0)
    check("Mmax = qL²/8 + PL/4", Mmax, 45.0 + 30.0)
    check("V(0) = 40", V[0], 40.0)
    E, I = 200e6, 5000e-8
    xs2, U, Vd = deformada_barra(m, r, b.id, 60)
    fmid = Vd[30]
    ref = -(5 * 10 * 6 ** 4 / (384 * E * I) + 20 * 6 ** 3 / (48 * E * I))
    check("flecha no meio (mm)", fmid * 1e3, ref * 1e3)


def balanco():
    print("2) Balanço L=3 m, P=10 kN na ponta, barra inclinada 30°")
    m = Modelo(); secao_gen(m)
    a = math.radians(30)
    n1, n2 = m.add_no(0, 0), m.add_no(3 * math.cos(a), 3 * math.sin(a))
    b = m.add_barra(n1.id, n2.id, "GEN")
    m.apoios[1] = Apoio(1, "fixo", "fixo", "fixo")
    m.add_caso("Q")
    m.cargas_no.append(CargaNo("Q", 2, 0.0, -10.0, 0.0))
    r = analisar(m)["Q"]
    xs, N, V, M = esforcos_barra(m, r, b.id)
    Lh = 3 * math.cos(a)
    check("Reação Mz engaste = P·Lh", r.reacoes[1][2], 10 * Lh)
    check("M(0) = −P·Lh (traciona em cima)", M[0], -10 * Lh)
    check("N = −P·sen30 (compressão)", N[0], -10 * math.sin(a))


def trelica_simples():
    print("3) Treliça triangular (vão 4 m, altura 2 m), P=10 kN no topo")
    m = Modelo(); secao_gen(m, 10.0, 100.0)
    n1, n2, n3 = m.add_no(0, 0), m.add_no(4, 0), m.add_no(2, 2)
    b1 = m.add_barra(1, 3, "GEN", True, True)
    b2 = m.add_barra(3, 2, "GEN", True, True)
    b3 = m.add_barra(1, 2, "GEN", True, True)
    m.apoios[1] = Apoio(1, "fixo", "fixo", "livre")
    m.apoios[2] = Apoio(2, "livre", "fixo", "livre")
    m.add_caso("Q")
    m.cargas_no.append(CargaNo("Q", 3, 0.0, -10.0, 0.0))
    r = analisar(m)["Q"]
    Nd = esforcos_barra(m, r, b1.id)[1][0]
    Nb = esforcos_barra(m, r, b3.id)[1][0]
    Mb = max(abs(x) for x in esforcos_barra(m, r, b1.id)[3])
    check("N diagonal = −P/(2·sen45)", Nd, -10 / (2 * math.sin(math.radians(45))))
    check("N banzo inferior = +P/2·cot45", Nb, 5.0)
    check("M na barra de treliça ≈ 0", Mb + 1.0, 1.0)


def portico_mola_comb():
    print("4) Viga biengastada / apoio em mola / combinação")
    m = Modelo(); secao_gen(m)
    m.add_no(0, 0); m.add_no(5, 0); m.add_no(10, 0)
    b = m.add_barra(1, 2, "GEN"); m.add_barra(2, 3, "GEN")
    m.apoios[1] = Apoio(1, "fixo", "fixo", "fixo")
    m.apoios[3] = Apoio(3, "fixo", "fixo", "fixo")
    m.add_caso("G"); m.add_caso("Q")
    for i in (1, 2):
        m.cargas_barra.append(CargaBarra("G", i, "distribuida", "GY", -12.0))
        m.cargas_barra.append(CargaBarra("Q", i, "distribuida", "GY", -4.0))
    m.combinacoes["ELU"] = Combinacao("ELU", {"G": 1.4, "Q": 1.5})
    res = analisar(m)
    xs, N, V, M = esforcos_barra(m, res["G"], 1)
    check("Biengastada L=10 (nó no meio): M engaste = −qL²/12", M[0], -12 * 100 / 12)
    check("M no meio = +qL²/24", M[-1], 12 * 100 / 24)
    q = 1.4 * 12 + 1.5 * 4
    xs, N, V, M = esforcos_barra(m, res["COMB: ELU"], 1)
    check("Combinação ELU: M engaste", M[0], -q * 100 / 12)
    # biengastada SEM nó interno (nenhum GL livre → solução analítica)
    m3 = Modelo(); secao_gen(m3); m3.add_no(0, 0); m3.add_no(6, 0); m3.add_barra(1, 2, "GEN")
    m3.apoios[1] = Apoio(1, "fixo", "fixo", "fixo"); m3.apoios[2] = Apoio(2, "fixo", "fixo", "fixo")
    m3.add_caso("P"); m3.cargas_barra.append(CargaBarra("P", 1, "pontual", "GY", -30.0, 2.0))
    r3 = analisar(m3)["P"]
    xs, N, V, M = esforcos_barra(m3, r3, 1)
    check("Biengastada s/ nós livres, P em a=2: M_A = −Pab²/L²", M[0], -30 * 2 * 16 / 36)
    check("  reação Ry_A = Pb²(3a+b)/L³", r3.reacoes[1][1], 30 * 16 * (6 + 4) / 216)
    # mola: balanço apoiado em mola na ponta
    m2 = Modelo(); secao_gen(m2)
    m2.add_no(0, 0); m2.add_no(4, 0); m2.add_barra(1, 2, "GEN")
    m2.apoios[1] = Apoio(1, "fixo", "fixo", "fixo")
    k = 1000.0
    m2.apoios[2] = Apoio(2, "livre", "mola", "livre", ky=k)
    m2.add_caso("P"); m2.cargas_no.append(CargaNo("P", 2, 0, -10, 0))
    r = analisar(m2)["P"]
    EI = 200e6 * 5000e-8
    kb = 3 * EI / 4 ** 3
    check("Reação na mola = P·k/(k+3EI/L³)", r.reacoes[2][1], 10 * k / (k + kb))


def peso_proprio_secoes():
    print("5) Peso próprio e propriedades de seções")
    m = Modelo()
    m.add_secao(Secao("W", "W", "ASTM A572 Gr50", {"catalogo": "W200x26,6"}))
    m.add_no(0, 0); m.add_no(6, 0); m.add_barra(1, 2, "W")
    m.apoios[1] = Apoio(1); m.apoios[2] = Apoio(2, "livre", "fixo")
    m.add_caso("PP", peso_proprio=True)
    r = analisar(m)["PP"]
    check("Ry total = γ·A·L", r.reacoes[1][1] + r.reacoes[2][1], 78.5 * 34.2e-4 * 6)
    A, _, _, Ix, _ = props_U(150, 60, 2.65, 20)
    t = 2.65
    check("Ue150x60x20x2,65: A (linha média) cm²", A / 100,
          t * ((150 - t) + 2 * (60 - t / 2) + 2 * (20 - t / 2)) / 100)
    A, _, yg, Ix, _ = props_L(50.8, 6.35)
    check("L 2\"x1/4\": A cm²", A / 100, 6.35 * (2 * 50.8 - 6.35) / 100)
    print(memorial_txt(m, analisar(m))[:300] + " ...")


def giro_e_pp():
    print("6) Giro do perfil (β), peso próprio desligado, composto e geometria")
    m = Modelo()
    m.add_secao(Secao("W", "W", "ASTM A572 Gr50", {"catalogo": "W200x26,6"}))
    m.add_no(0, 0); m.add_no(3, 0); m.add_no(0, 2); m.add_no(3, 2)
    m.add_barra(1, 2, "W")                       # β = 0
    b2 = m.add_barra(3, 4, "W"); b2.beta = 90.0  # β = 90 → Iy
    m.apoios[1] = Apoio(1, "fixo", "fixo", "fixo"); m.apoios[3] = Apoio(3, "fixo", "fixo", "fixo")
    m.add_caso("P"); m.cargas_no += [CargaNo("P", 2, 0, -1, 0), CargaNo("P", 4, 0, -1, 0)]
    r = analisar(m)["P"]
    E = 200e6
    check("β=0: flecha = PL³/(3E·Ix)", r.desloc[2][1], -27 / (3 * E * 2611e-8))
    check("β=90: flecha = PL³/(3E·Iy)", r.desloc[4][1], -27 / (3 * E * 330e-8))
    # peso próprio desligado em uma barra
    m.add_caso("PP", peso_proprio=True); b2.pp = False
    r = analisar(m)["PP"]
    check("PP só na barra com pp=True", r.reacoes[1][1] + r.reacoes[3][1], 78.5 * 34.2e-4 * 3)
    # composto 2L com afastamento: confere com integração do polígono
    s = Secao("2L", "L", "ASTM A36", {"b": 50.8, "t": 6.35, "gap": 8.0}, qtd=2)
    p = s.propriedades()
    pols = s.poligonos(0)
    Iy_pol = 0.0
    for q in pols:          # Iy = ∫z² dA  (Green)
        for (z1, y1), (z2, y2) in zip(q, q[1:] + q[:1]):
            cr = z1 * y2 - z2 * y1
            Iy_pol += (z1 * z1 + z1 * z2 + z2 * z2) * cr / 12
    check("2L 2\"x1/4\" gap 8: Iy composto (cm⁴)", p["Iy"], abs(Iy_pol) / 1e4)
    check("2L: Ix composto = 2·Ix1", p["Ix"], 2 * p["Ix1"])
    ymin, ymax, ints = Secao("W", "W", dims={"catalogo": "W200x26,6"}).linhas_elevacao(0)
    check("Elevação W200 β=0: altura = d", ymax - ymin, 207.0)
    check("Elevação W200 β=0: 2 arestas internas (mesas)", len(ints), 2)
    ymin, ymax, ints = Secao("W", "W", dims={"catalogo": "W200x26,6"}).linhas_elevacao(90)
    check("Elevação W200 β=90: altura = bf", ymax - ymin, 133.0)


def combinacoes_e_dxf():
    print("7) Combinações NBR 8800 e importação DXF")
    import combinacoes_nbr8800 as nbr, importar_dxf as dxf, os, tempfile
    m = Modelo()
    for n in ["G — permanente", "Q — sobrecarga", "V1 — vento 0°", "V2 — vento 90°"]:
        m.add_caso(n)
    nbr.garantir_classificacao(m)
    L = nbr.gerar(m, elu=True, rara=True, frequente=True, quase_perm=True)
    elu = {c.descricao: c for c in L if c.tipo == "ELU"}
    check("Nº de ELU (G;Q;V1|V2 exclusivos; c/ perm. favoráveis)", len(elu), 15)
    c = [c for c in L if c.tipo == "ELU" and set(c.fatores) == {"G — permanente", "Q — sobrecarga", "V1 — vento 0°"}
         and c.fatores["G — permanente"] == 1.25][0]
    check("1,25G + 1,5Q + 1,4·0,6·V1: fator do vento", c.fatores["V1 — vento 0°"], 0.84)
    check("Nenhuma ELU com V1 e V2 juntos", sum(1 for c in L if "V1 — vento 0°" in c.fatores
                                                  and "V2 — vento 90°" in c.fatores), 0)
    fr = [c for c in L if c.tipo == "ELS-F" and "V1 — vento 0°" in c.fatores and "Q — sobrecarga" in c.fatores][0]
    check("Frequente: G + ψ1·V1 + ψ2·Q → ψ2 da cobertura", fr.fatores["Q — sobrecarga"], 0.6)
    # DXF ASCII mínimo (LINE contínua + LWPOLYLINE), em mm
    txt = ["0", "SECTION", "2", "ENTITIES",
           "0", "LINE", "8", "BANZO", "10", "0", "20", "0", "11", "6000", "21", "0",
           "0", "LWPOLYLINE", "8", "DIAG", "90", "3", "70", "0", "10", "0", "20", "0",
           "10", "3000", "20", "1500", "10", "6000", "20", "0",
           "0", "LINE", "8", "MONT", "10", "3000", "20", "0", "11", "3000", "21", "1500",
           "0", "ENDSEC", "0", "EOF"]
    p = os.path.join(tempfile.gettempdir(), "t.dxf")
    open(p, "w").write("\n".join(txt))
    segs = dxf._ler_ascii(p)
    m2 = Modelo(); secao_gen(m2)
    rel = dxf.importar(m2, segs, 0.001, {k: ("GEN", True, True, 0) for k in ("BANZO", "DIAG", "MONT")}, tol=1e-3)
    check("DXF: banzo contínuo quebrado no montante → 5 barras", rel["barras_novas"], 5)
    check("DXF: 4 nós", len(m2.nos), 4)


def nos_sobre_barra():
    print("8) Nós lançados sobre barra existente (cargas nodais)")
    m = Modelo(); secao_gen(m)
    m.add_no(0, 0); m.add_no(10, 0); m.add_barra(1, 2, "GEN")
    m.apoios[1] = Apoio(1); m.apoios[2] = Apoio(2, "livre", "fixo")
    m.add_caso("G")
    m.cargas_barra += [CargaBarra("G", 1, "distribuida", "GY", -1.0), CargaBarra("G", 1, "pontual", "GY", -3.0, 8.0)]
    m.add_no(6.4, 0); m.add_no(8.8, 0)
    m.cargas_no += [CargaNo("G", 3, 0, -2.38, 0), CargaNo("G", 4, 0, -2.38, 0)]
    check("Carga em nó solto é acusada (não ignorada)", len(m.validar()[0]), 2)
    check("Divisão automática: 2 divisões", m.inserir_nos_nas_barras(), 2)
    r = analisar(m)["G"]
    check("Ry1 = q·L/2 + ΣP·b/L (cargas nodais + pontual realocada)", r.reacoes[1][1],
          5 + 2.38 * 3.6 / 10 + 2.38 * 1.2 / 10 + 3 * 2 / 10)
    m2 = Modelo(); secao_gen(m2); m2.add_no(0, 0); m2.add_no(6, 0); m2.add_barra(1, 2, "GEN", True, False)
    check("Dividir em 3 partes → 3 barras", m2.dividir_em_partes(1, 3) + 1, len(m2.barras))
    check("Rótula inicial preservada só no 1º trecho", sum(b.rot_i for b in m2.barras.values()), 1)


def nbr8800_baseline_cype():
    print("9) Dimensionamento NBR 8800 — baseline memorial de referência W310x21,0 (L = 1,50 m)")
    import verificacao_nbr8800 as vf
    from portico2d_core import Material
    t = 9.80665                      # tf → kN
    m = Modelo()
    m.materiais["A36 ref"] = Material("A36 ref", 2038736 * 0.0980665, 784913 * 0.0980665, 2548.42 * 0.0980665, 400)
    m.add_secao(Secao("W310x21,0", "W", "A36 ref",
                      {"catalogo": "W310x21,0", "Zx": 292.0, "Zy": 31.0, "J": 3.27, "Cw": 21628.0}))
    m.projeto = {"h_alma": "simplificado", "sigma_qa": "fy"}     # mesmas hipóteses do memorial de referência
    c = vf.cfg(m)
    p = vf.props_I(m.secoes["W310x21,0"], c)
    mt = m.materiais["A36 ref"]
    mat = {"E": mt.E / 10, "G": mt.G / 10, "fy": mt.fy / 10, "fu": mt.fu / 10}
    tol = 3e-3                       # catálogo A = 27,20 cm² × referência 27,21 cm²
    check("Wx (cm³)", p["Wx"], 249.24, tol)
    check("Nt,Rd (tf)", vf.tracao(p, mat, c, vf.PROJ_PADRAO)[0] / t, 63.039, tol)
    Nc, _, d = vf.compressao(p, mat, c, 150, 150, 0)
    check("Ne = Ney (tf)", d["Ne"] / t, 87.640, tol)
    check("Q = Qs·Qa", d["Q"], 0.885, tol)
    check("χ", d["chi"], 0.746, tol)
    check("Nc,Rd (tf)", Nc / t, 41.625, tol)
    check("Mx,Rd (tf·m) — Lb = 0", vf.flexao_x(p, mat, c, 0, 1.0)[0] / 100 / t, 6.765, tol)
    check("My,Rd (tf·m)", vf.flexao_y(p, mat, c)[0] / 100 / t, 0.450, tol)
    check("Vy,Rd alma (tf)", vf.cortante(p, mat, c, True)[0] / t, 21.480, tol)
    check("Vx,Rd mesas (tf)", vf.cortante(p, mat, c, False)[0] / t, 16.005, tol)
    # FLT: W200x15, L = 6 m, carga uniforme → Cb = 1,136
    xs = [6 * k / 40 for k in range(41)]
    M = [x * (6 - x) for x in xs]
    check("Cb carga uniforme (biapoiada)", vf.calcular_cb(xs, M, 6, 6, False)[0], 12.5 / 11.0, 1e-3)


def kl_manual_e_flecha_absoluta():
    print("10) KL manual em barra dividida por nós de carga e flecha absoluta")
    import verificacao_nbr8800 as vf
    m = Modelo()
    m.add_secao(Secao("W", "W", "ASTM A572 Gr50", {"catalogo": "W310x32,7"}))
    for x in (0, 2, 4, 6):
        m.add_no(x, 0)
    for i in (1, 2, 3):
        m.add_barra(i, i + 1, "W")
    m.apoios[1] = Apoio(1); m.apoios[4] = Apoio(4, "livre", "fixo")
    m.add_caso("Q")
    m.cargas_no += [CargaNo("Q", 2, 0, -10, 0), CargaNo("Q", 3, 0, -10, 0)]
    res = analisar(m)
    b = m.barras[2]
    b.proj = {"Lf": 6.0, "Lb": 6.0, "Lz": 6.0, "flecha_modo": "absoluta", "Lflecha": 6.0}
    comp, L = vf.comprimentos(m, b)
    check("KL fora do plano manual = 6 m (barra de 2 m)", comp["f"][0] / 100, 6.0)
    check("Lb manual = 6 m", comp["b"][0] / 100, 6.0)
    v = vf.verificar_barra(m, res, 2)
    fl = [it for it in v.itens if it["nome"].startswith("Flecha")][0]
    E, I, P, a, Lt = 200e6, 6570e-8, 10.0, 2.0, 6.0
    ref = P * a * (3 * Lt ** 2 - 4 * a ** 2) / (24 * E * I) * 1000   # 2 cargas simétricas, flecha no meio
    check("Flecha absoluta no meio do vão (mm)", fl["sd"], ref, 2e-3)
    check("Limite com Lref = 6 m: 6000/250 (mm)", fl["rd"], 24.0)
    comp2, _ = vf.comprimentos(m, m.barras[1])
    check("Barra sem KL manual continua K×L = 2 m", comp2["f"][0] / 100, 2.0)


def reacoes_apoios_rotulados():
    print("11) Reações em apoios que só recebem barras rotuladas (equilíbrio global)")
    m = Modelo(); secao_gen(m, 10.0, 100.0)
    m.add_no(0, 0); m.add_no(4, 0); m.add_no(2, 2)
    for a, b in ((1, 3), (3, 2), (1, 2)):
        m.add_barra(a, b, "GEN", True, True)
    m.apoios[1] = Apoio(1); m.apoios[2] = Apoio(2, "livre", "fixo"); m.add_caso("Q")
    m.cargas_no.append(CargaNo("Q", 3, 3, -10, 0))
    r = analisar(m)["Q"]
    check("Treliça: Ry2 = (10·2 + 3·2)/4", r.reacoes[2][1], 6.5)
    check("Treliça: Rx1 = −3", r.reacoes[1][0], -3.0)
    # caso do usuário: viga em balanço atirantada por diagonal rotulada até apoio
    m = Modelo(); secao_gen(m)
    for p in [(0, 0), (1.6, 0), (3.1, 0), (4.4, 0), (4.4, 2.3)]:
        m.add_no(*p)
    for a, b in ((1, 2), (2, 3), (3, 4)):
        m.add_barra(a, b, "GEN")
    m.add_barra(2, 5, "GEN", True, True)
    m.apoios[4] = Apoio(4, "fixo", "fixo", "livre"); m.apoios[5] = Apoio(5, "fixo", "fixo", "livre")
    m.add_caso("G", peso_proprio=True)
    m.cargas_no += [CargaNo("G", 1, 0, -7.57, 0), CargaNo("G", 3, 0, -7.57, 0)]
    r = analisar(m)["G"]
    peso = sum(78.5 * m.secoes["GEN"].A_m2() * m.geometria(b)[0] for b in m.barras.values())
    check("ΣRy = cargas + peso próprio", sum(R[1] for R in r.reacoes.values()), 2 * 7.57 + peso)
    check("ΣRx = 0", sum(R[0] for R in r.reacoes.values()) + 1.0, 1.0)
    check("Apoio do tirante recebe reação (≠ 0)", abs(r.reacoes[5][1]) > 1.0, True)


if __name__ == "__main__":
    viga_biapoiada(); balanco(); trelica_simples(); portico_mola_comb(); peso_proprio_secoes(); giro_e_pp(); combinacoes_e_dxf(); nos_sobre_barra(); nbr8800_baseline_cype(); kl_manual_e_flecha_absoluta(); reacoes_apoios_rotulados()
    print(f"\n{sum(OK)}/{len(OK)} verificações OK")
