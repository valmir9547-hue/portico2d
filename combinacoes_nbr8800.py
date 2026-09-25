# -*- coding: utf-8 -*-
"""
combinacoes_nbr8800.py — Geração automática de combinações de ações (sem GUI)
==============================================================================
Base normativa: ABNT NBR 8800 (seção 4.7 — combinações de ações; Tabela 1 — γf;
Tabela 2 — fatores de combinação ψ0 e de redução ψ1, ψ2), em conformidade com a
NBR 8681.  Os valores abaixo reproduzem as Tabelas 1 e 2 da NBR 8800; CONFIRA-OS
com o texto da edição 2024 em uso no escritório — todos são editáveis por caso.

ELU — combinações últimas normais (e especiais/de construção, com os γ próprios):
    Fd = Σ γg,i·FGi,k  +  γq1·FQ1,k  +  Σ(j≥2) γqj·ψ0j·FQj,k
ELS — combinações de serviço:
    quase permanente:  Σ FGi,k + Σ ψ2j·FQj,k
    frequente:         Σ FGi,k + ψ1·FQ1,k + Σ(j≥2) ψ2j·FQj,k
    rara:              Σ FGi,k + FQ1,k + Σ(j≥2) ψ1j·FQj,k

Regras de montagem adotadas:
  • Cada ação variável é, por vez, a principal (FQ1).
  • Casos com o mesmo "grupo" são mutuamente exclusivos (ex.: Vento 0°, Vento 90°,
    Vento sucção): no máximo um caso do grupo entra em cada combinação.
  • Ações variáveis favoráveis não são consideradas: as secundárias entram ou não
    (todas as possibilidades), se a opção "secundárias opcionais" estiver ligada.
  • Permanentes: todos com γ desfavorável; opcionalmente, uma segunda família com
    todos os permanentes com γ favorável (situação de sucção de vento, tombamento).
"""
from __future__ import annotations

import itertools

from portico2d_core import Combinacao, CasoCarga

# ---------------------------------------------------------------------------
# Tabela 1 — Valores dos coeficientes de ponderação das ações (γf = γf1·γf3)
# permanentes: categoria -> {tipo_comb: (γ desfavorável, γ favorável)}
# ---------------------------------------------------------------------------
PERMANENTES = {
    "Peso próprio de estruturas metálicas":
        {"normal": (1.25, 1.00), "especial": (1.15, 1.00)},
    "Peso próprio de estruturas pré-moldadas":
        {"normal": (1.30, 1.00), "especial": (1.20, 1.00)},
    "Peso próprio de estruturas moldadas no local e de elementos construtivos industrializados; empuxos permanentes":
        {"normal": (1.35, 1.00), "especial": (1.25, 1.00)},
    "Peso próprio de elementos construtivos industrializados com adições in loco":
        {"normal": (1.40, 1.00), "especial": (1.30, 1.00)},
    "Peso próprio de elementos construtivos em geral e equipamentos":
        {"normal": (1.50, 1.00), "especial": (1.40, 1.00)},
    "Ações permanentes indiretas (recalques, retração)":
        {"normal": (1.20, 0.00), "especial": (1.20, 0.00)},
}

# variáveis: categoria -> (chave de γq, ψ0, ψ1, ψ2)
GAMA_Q = {  # Tabela 1 — ações variáveis
    "temperatura": {"normal": 1.20, "especial": 1.00},
    "vento":       {"normal": 1.40, "especial": 1.20},
    "truncadas":   {"normal": 1.20, "especial": 1.10},
    "demais":      {"normal": 1.50, "especial": 1.30},
}
VARIAVEIS = {  # Tabela 2 — ψ0, ψ1, ψ2
    "Sobrecarga — locais sem predominância de pesos/equipamentos fixos nem elevadas concentrações de pessoas":
        ("demais", 0.5, 0.4, 0.3),
    "Sobrecarga — locais com predominância de pesos/equipamentos fixos ou elevadas concentrações de pessoas":
        ("demais", 0.7, 0.6, 0.4),
    "Sobrecarga — bibliotecas, arquivos, depósitos, oficinas, garagens e sobrecargas em coberturas":
        ("demais", 0.8, 0.7, 0.6),
    "Vento — pressão dinâmica nas estruturas em geral":
        ("vento", 0.6, 0.3, 0.0),
    "Temperatura — variações uniformes":
        ("temperatura", 0.6, 0.5, 0.3),
    "Cargas móveis — passarelas de pedestres":
        ("demais", 0.6, 0.4, 0.3),
    "Cargas móveis — vigas de rolamento de pontes rolantes":
        ("demais", 1.0, 0.8, 0.5),
    "Cargas móveis — pilares e elementos que suportam vigas de rolamento":
        ("demais", 0.7, 0.6, 0.4),
    "Ações truncadas":
        ("truncadas", 0.8, 0.7, 0.6),
}

NATUREZAS = ("permanente", "variavel", "ignorar")


def classificar_por_nome(caso: CasoCarga):
    """Sugere natureza/categoria/grupo a partir do nome (G, PP, Q, SC, V, T...)."""
    n = caso.nome.strip().lower()
    if n.startswith(("g", "pp", "perm", "cp")) or "permanente" in n or "peso" in n:
        return "permanente", "Peso próprio de estruturas metálicas", ""
    if n.startswith("v") or "vento" in n:
        return "variavel", "Vento — pressão dinâmica nas estruturas em geral", "Vento"
    if n.startswith("t") or "temperat" in n:
        return "variavel", "Temperatura — variações uniformes", "Temperatura"
    return "variavel", "Sobrecarga — bibliotecas, arquivos, depósitos, oficinas, garagens e sobrecargas em coberturas", ""


def aplicar_categoria(caso: CasoCarga, natureza, categoria, tipo_comb="normal"):
    """Preenche γ e ψ do caso a partir das tabelas normativas."""
    caso.natureza, caso.categoria = natureza, categoria
    if natureza == "permanente" and categoria in PERMANENTES:
        caso.gd, caso.gf = PERMANENTES[categoria][tipo_comb]
        caso.gq = caso.psi0 = caso.psi1 = caso.psi2 = 0.0
    elif natureza == "variavel" and categoria in VARIAVEIS:
        chave, p0, p1, p2 = VARIAVEIS[categoria]
        caso.gq = GAMA_Q[chave][tipo_comb]
        caso.psi0, caso.psi1, caso.psi2 = p0, p1, p2
        caso.gd = caso.gf = 0.0


def garantir_classificacao(modelo, tipo_comb="normal"):
    for c in modelo.casos.values():
        if not c.natureza:
            nat, cat, grp = classificar_por_nome(c)
            aplicar_categoria(c, nat, cat, tipo_comb)
            c.grupo = grp


def _f(v):
    return f"{v:.3f}".rstrip("0").rstrip(".").replace(".", ",")


def _descr(fatores):
    return " + ".join(f"{_f(f)}·{c}" for c, f in fatores.items())


def gerar(modelo, tipo_comb="normal", elu=True, rara=False, frequente=False, quase_perm=False,
          perm_favoravel=True, secundarias_opcionais=True, limite=2000):
    """Retorna lista de Combinacao (auto=True). Não altera o modelo."""
    perm = [c for c in modelo.casos.values() if c.natureza == "permanente"]
    var = [c for c in modelo.casos.values() if c.natureza == "variavel"]
    grupos = {}
    for c in var:
        grupos.setdefault(c.grupo.strip() or f"__{c.nome}", []).append(c)

    def opcoes_secundarias(excluir_grupo):
        listas = []
        for g, cs in grupos.items():
            if g == excluir_grupo:
                continue
            listas.append(([None] if secundarias_opcionais else []) + cs)
        return itertools.product(*listas) if listas else [()]

    def grupo_de(c):
        return c.grupo.strip() or f"__{c.nome}"

    saida, vistos = [], set()

    def add(tipo, fat, prefixo):
        fat = {k: round(v, 6) for k, v in fat.items() if abs(v) > 1e-12}
        if not fat:
            return
        chave = (tipo, tuple(sorted(fat.items())))
        if chave in vistos:
            return
        vistos.add(chave)
        n = sum(1 for c in saida if c.tipo == tipo) + 1
        saida.append(Combinacao(f"{prefixo}{n:02d}", fat, tipo, _descr(fat), True))
        if len(saida) > limite:
            raise ValueError(f"Mais de {limite} combinações — use grupos de exclusão ou desligue "
                             "'secundárias opcionais'.")

    if elu:
        familias = [("desf", {c.nome: c.gd for c in perm})]
        if perm_favoravel and any(abs(c.gf - c.gd) > 1e-9 for c in perm):
            familias.append(("fav", {c.nome: c.gf for c in perm}))
        for fam, base in familias:
            if fam == "desf":
                add("ELU", dict(base), "ELU")                   # só permanentes
            for q1 in var:
                for sec in opcoes_secundarias(grupo_de(q1)):
                    fat = dict(base)
                    fat[q1.nome] = q1.gq
                    for qj in sec:
                        if qj is not None:
                            fat[qj.nome] = qj.gq * qj.psi0
                    add("ELU", fat, "ELU")
    base_s = {c.nome: 1.0 for c in perm}
    if quase_perm:
        for sec in opcoes_secundarias(None):
            fat = dict(base_s)
            for qj in sec:
                if qj is not None:
                    fat[qj.nome] = qj.psi2
            add("ELS-QP", fat, "QP")
    for ativo, tipo, pref, f1, fj in ((frequente, "ELS-F", "FREQ", "psi1", "psi2"),
                                      (rara, "ELS-R", "RARA", None, "psi1")):
        if not ativo:
            continue
        for q1 in var:
            for sec in opcoes_secundarias(grupo_de(q1)):
                fat = dict(base_s)
                fat[q1.nome] = getattr(q1, f1) if f1 else 1.0
                for qj in sec:
                    if qj is not None:
                        fat[qj.nome] = getattr(qj, fj)
                add(tipo, fat, pref)
        if not var:
            add(tipo, dict(base_s), pref)
    return saida


def aplicar(modelo, combinacoes, substituir_auto=True):
    """Insere as combinações no modelo (remove as geradas anteriormente, se pedido)."""
    if substituir_auto:
        for k in [k for k, c in modelo.combinacoes.items() if c.auto]:
            modelo.combinacoes.pop(k)
    for c in combinacoes:
        nome, i = c.nome, 1
        while nome in modelo.combinacoes:
            i += 1
            nome = f"{c.nome}_{i}"
        c.nome = nome
        modelo.combinacoes[nome] = c
