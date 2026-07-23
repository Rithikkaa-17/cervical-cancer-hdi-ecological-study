import graphviz

dot = graphviz.Digraph(format="png")
dot.attr(rankdir="LR", bgcolor="white", nodesep="0.45", ranksep="0.75", dpi="230")
dot.attr("node", shape="ellipse", style="filled", fontname="Helvetica", fontsize="11.5",
         fillcolor="#eaf0fb", color="#2E5090", penwidth="1.5")
dot.attr("edge", color="#333333", penwidth="1.3", fontname="Helvetica", fontsize="9.5")

dot.node("HDI", "HDI\n(Development)", fillcolor="#dbe6fa", width="1.4", height="0.9")
dot.node("HIV", "HIV\nPrevalence", fillcolor="#dff0e5", color="#4a9c6d")
dot.node("SMK", "Smoking\nPrevalence", fillcolor="#fff3d6", color="#c99a1e")
dot.node("VAX", "HPV Vax.\nCoverage", fillcolor="#fbeaea", color="#c0392b")
dot.node("SCR", "Screening\nProgram", fillcolor="#fbeaea", color="#c0392b")
dot.node("ASIR", "log(ASIR)\nIncidence", fillcolor="#eaf0fb", width="1.3", height="0.9")
dot.node("ASMR", "log(ASMR)\nMortality", fillcolor="#eaf0fb", width="1.3", height="0.9")

# Direct hypothesized causal paths (solid)
dot.edge("HDI", "ASIR", label="direct −", color="#2E5090", penwidth="1.8")
dot.edge("HDI", "ASMR", label="direct −\n(treatment access)", color="#2E5090", penwidth="1.8")
dot.edge("HDI", "HIV", label="−", color="#4a9c6d")
dot.edge("HIV", "ASMR", label="+ (biological\ncofactor)", color="#4a9c6d")
dot.edge("HDI", "SMK", label="+ (confound)", color="#c99a1e", style="dashed")
dot.edge("SMK", "ASIR", label="ns after\nadjustment", color="#c99a1e", style="dashed")

# Reverse-causation / policy-targeting paths (dashed, red) - the counter-intuitive finding
dot.edge("ASIR", "VAX", label="reverse causation /\npolicy targeting", color="#c0392b", style="dashed")
dot.edge("ASIR", "SCR", label="reverse causation /\npolicy targeting", color="#c0392b", style="dashed")
dot.edge("VAX", "ASIR", label="+ (spurious,\nSection VI)", color="#c0392b", style="dotted")
dot.edge("SCR", "ASIR", label="+ (spurious,\nSection VI)", color="#c0392b", style="dotted")

dot.render("fig_dag", cleanup=True)
print("saved")
