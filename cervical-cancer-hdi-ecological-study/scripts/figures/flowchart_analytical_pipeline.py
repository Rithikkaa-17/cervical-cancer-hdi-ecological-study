import graphviz

dot = graphviz.Digraph(format="png")
dot.attr(rankdir="TB", bgcolor="white", nodesep="0.35", ranksep="0.4", dpi="220")
dot.attr("node", shape="box", style="rounded,filled", fontname="Helvetica", fontsize="12",
         fillcolor="#eaf0fb", color="#2E5090", penwidth="1.5", margin="0.18,0.12")
dot.attr("edge", color="#555555", penwidth="1.3", arrowsize="0.8")

dot.node("A", "Data Collection\n(GLOBOCAN 2024/2022, UNDP HDI,\nWHO/UNICEF, UNAIDS, OWID)")
dot.node("B", "Merge on ISO-3 Code\n(186 -> 176 countries)")
dot.node("C", "Bivariate Analysis\n(Spearman correlation, HDI tiers)")
dot.node("D", "Multivariable OLS\n(HDI + HPV vax + Screening +\nSmoking [+ HIV], VIF screening)")

dot.node("E", "Robustness Checks", shape="plaintext", fontsize="13", fontname="Helvetica-Bold")

dot.attr("node", fillcolor="#dff0e5", color="#4a9c6d", fontsize="10.5")
dot.node("E1", "Temporal Panel\n(2022 vs 2024,\npaired n=156)")
dot.node("E2", "HIV + Data-Quality\nStratification\n(n=131, n=142)")
dot.node("E3", "Mediation Analysis\n(HDI -> HIV -> ASMR)")
dot.node("E4", "Spatial-Lag\nCorrection\n(k-NN, k=3,5,8,10)")

dot.attr("node", fillcolor="#eaf0fb", color="#2E5090", fontsize="12")
dot.node("F", "Synthesis & Reporting\n(Sections IV-VI)")

dot.edge("A", "B")
dot.edge("B", "C")
dot.edge("C", "D")
dot.edge("D", "E")
for n in ["E1","E2","E3","E4"]:
    dot.edge("E", n)
    dot.edge(n, "F")

dot.render("fig_pipeline_gv", cleanup=True)
print("saved fig_pipeline_gv.png")
