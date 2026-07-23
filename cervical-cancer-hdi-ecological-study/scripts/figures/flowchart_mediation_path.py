import graphviz

dot = graphviz.Digraph(format="png")
dot.attr(rankdir="LR", bgcolor="white", nodesep="0.6", ranksep="0.9", dpi="220")
dot.attr("node", shape="ellipse", style="filled", fontname="Helvetica", fontsize="13",
         fillcolor="#eaf0fb", color="#2E5090", penwidth="1.6", width="1.6", height="0.9")
dot.attr("edge", color="#333333", penwidth="1.5", fontname="Helvetica", fontsize="11.5")

dot.node("HDI", "HDI\n(Development)")
dot.node("HIV", "HIV Prevalence\n(Mediator)", fillcolor="#dff0e5", color="#4a9c6d")
dot.node("ASMR", "log(ASMR)\n(Mortality)")

dot.edge("HDI", "HIV", label="  a = -5.55\n  (p = 0.005)")
dot.edge("HIV", "ASMR", label="  b = +0.081\n  (p < 0.001)")
dot.edge("HDI", "ASMR", label="  c' = -5.66 (direct, p < 0.001)\n  c = -6.34 (total, p < 0.001)", constraint="false")

dot.render("fig_mediation_path", cleanup=True)
print("saved")
