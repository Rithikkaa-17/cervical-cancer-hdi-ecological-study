"""Fig. 1: analytical pipeline. The spatial dependence analysis is the primary
analysis; the temporal, HIV/data-quality and mediation analyses are supporting checks.
Run from the scripts/figures folder (or anywhere): writes fig01_analytical_pipeline.png
into the project's figures/ folder."""
import os
import graphviz

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                   "figures", "fig01_analytical_pipeline")

dot = graphviz.Digraph(format="png")
dot.attr(rankdir="TB", bgcolor="white", nodesep="0.35", ranksep="0.4", dpi="220")
dot.attr("node", shape="box", style="rounded,filled", fontname="Helvetica", fontsize="12",
         fillcolor="#eaf0fb", color="#2E5090", penwidth="1.5", margin="0.18,0.12")
dot.attr("edge", color="#555555", penwidth="1.3", arrowsize="0.8")

dot.node("A", "Data Collection\n(GLOBOCAN 2024/2022, UNDP HDI,\nWHO/UNICEF, UNAIDS, OWID)")
dot.node("B", "Merge on ISO-3 Code\n(186 -> 176 countries)")
dot.node("C", "Bivariate Analysis\n(Spearman correlation, HDI tiers)")
dot.node("D", "Multivariable OLS\n(HDI + HPV vax + Screening +\nSmoking [+ HIV], VIF screening)")

# Primary analysis: highlighted
dot.node("S", "PRIMARY ANALYSIS: Spatial Dependence (Section V-A)\n"
              "Moran's I on OLS residuals -> robust LM tests -> ML spatial lag model\n"
              "Direct / spillover / total HDI impacts\n"
              "9 specifications compared (Table V): k-NN k=3,5,8,10;\n"
              "inverse distance; distance band; Queen (+/- islands); spatial error",
         fillcolor="#fde9d9", color="#C0504D", penwidth="2.4", fontsize="12.5",
         fontname="Helvetica-Bold", margin="0.25,0.15")

dot.node("E", "Supporting Robustness Checks", shape="plaintext", fontsize="13", fontname="Helvetica-Bold")
dot.attr("node", fillcolor="#dff0e5", color="#4a9c6d", fontsize="10.5")
dot.node("E1", "Temporal Panel (V-B)\n(2022 vs 2024,\npaired n=156 / 175)")
dot.node("E2", "HIV + Data-Quality\nStratification (V-C)\n(n=131, n=142)")
dot.node("E3", "Mediation Analysis (V-D)\n(HDI -> HIV -> ASMR)")

dot.attr("node", fillcolor="#eaf0fb", color="#2E5090", fontsize="12")
dot.node("F", "Synthesis & Reporting\n(Sections IV-VI)")

dot.edge("A", "B")
dot.edge("B", "C")
dot.edge("C", "D")
dot.edge("D", "S", color="#C0504D", penwidth="2.0")
dot.edge("S", "E")
for n in ["E1", "E2", "E3"]:
    dot.edge("E", n)
    dot.edge(n, "F")

dot.render(OUT, cleanup=True)
print("saved", OUT + ".png")
