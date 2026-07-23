import graphviz

dot = graphviz.Digraph(format="png")
dot.attr(rankdir="TB", bgcolor="white", nodesep="0.5", ranksep="0.45", dpi="240", splines="ortho")
dot.attr("node", shape="box", style="filled,rounded", fontname="Helvetica", fontsize="12",
         fillcolor="#eaf0fb", color="#2E5090", penwidth="1.6", margin="0.22,0.16")
dot.attr("edge", color="#444444", penwidth="1.4", arrowsize="0.85")

# Main vertical spine
dot.node("N0", "GLOBOCAN 2024 countries/territories\nwith cervical cancer estimates\n(n = 186)", width="3.6")
dot.node("N1", "Core analytic sample\n(HDI available)\n(n = 176)", width="3.6")
dot.node("N2", "Multivariable base-model sample\n(HDI + HPV vaccination + screening)\n(n = 175)", width="3.6")
dot.node("N3", "Full-covariate model sample\n(+ female smoking prevalence)\n(n = 156)", width="3.6")
dot.node("N4", "HIV / mediation / spatial\nanalytic sample\n(+ HIV prevalence)\n(n = 131)", width="3.6")

# Exclusion boxes (side branches)
dot.attr("node", shape="box", style="filled,rounded", fillcolor="#fbeaea", color="#c0392b", fontsize="10.5", margin="0.18,0.12")
dot.node("X0", "Excluded (n = 10):\nsmall territories without UNDP HDI\n(e.g. Puerto Rico, Guadeloupe, Guam),\nNorth Korea, GLOBOCAN world-\naggregate row", width="3.0")
dot.node("X1", "Excluded (n = 1):\nmissing HPV vaccination coverage\nor screening-program status", width="3.0")
dot.node("X2", "Excluded (n = 19):\nmissing female smoking prevalence\n(WHO GHO)", width="3.0")
dot.node("X3", "Excluded (n = 25):\nmissing HIV prevalence\n(UNAIDS)", width="3.0")

# Registry-quality side sample (separate branch off N1)
dot.attr("node", shape="box", style="filled,rounded", fillcolor="#dff0e5", color="#4a9c6d", fontsize="10.5", margin="0.18,0.12")
dot.node("N5", "Registry-quality stratification sample\n(death-registration completeness\navailable, Karlinsky 2024)\n(n = 142)", width="3.2")

# Rank pairs to align exclusion boxes beside their corresponding spine node
with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("N0"); s.node("X0")
with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("N1"); s.node("X1"); s.node("N5")
with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("N2"); s.node("X2")
with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("N3"); s.node("X3")

dot.edge("N0", "N1")
dot.edge("N1", "N2")
dot.edge("N2", "N3")
dot.edge("N3", "N4")
dot.edge("N0", "X0")
dot.edge("N1", "X1")
dot.edge("N2", "X2")
dot.edge("N3", "X3")
dot.edge("N1", "N5")

dot.render("fig_strobe_v2", cleanup=True)
print("saved")
