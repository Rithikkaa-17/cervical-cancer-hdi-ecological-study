import graphviz

dot = graphviz.Digraph(format="png")
dot.attr(rankdir="TB", bgcolor="white", nodesep="0.4", ranksep="0.4", dpi="230")
dot.attr("node", shape="box", style="filled,rounded", fontname="Helvetica", fontsize="11.5",
         fillcolor="#eaf0fb", color="#2E5090", penwidth="1.5", margin="0.2,0.14")
dot.attr("edge", color="#444444", penwidth="1.3", fontname="Helvetica", fontsize="10")

dot.node("A", "Fit OLS: log(ASIR/ASMR) ~\nHDI + HPV vax + Screening +\nSmoking [+ HIV]")
dot.node("B", "Build k-NN spatial weights\nmatrix (k = 3, 5, 8, 10)\nfrom country centroids")
dot.node("C", "Compute global Moran's I\non OLS residuals\n(999 permutations)")

dot.node("D", "Moran's I significant?\n(p < 0.05)", shape="diamond", style="filled",
         fillcolor="#fff3d6", color="#c99a1e", fontsize="11", width="2.2", height="1.3")

dot.node("E1", "Spatial independence\nnot rejected —\nOLS retained as final model",
         fillcolor="#dff0e5", color="#4a9c6d")

dot.node("F", "Run Lagrange Multiplier\ndiagnostics\n(LM-lag vs. LM-error,\nrobust versions)")

dot.node("G", "Which robust LM\nis significant?", shape="diamond", style="filled",
          fillcolor="#fff3d6", color="#c99a1e", fontsize="11", width="2.2", height="1.3")

dot.node("H1", "Robust LM-lag\nsignificant\n(this study's result)", fillcolor="#eaf0fb")
dot.node("H2", "Robust LM-error\nsignificant")
dot.node("H3", "Both / neither\nclearly dominant")

dot.node("I1", "Fit maximum-likelihood\nSpatial Lag model (ML_Lag)\nSection V-D result", fillcolor="#dbe6fa", color="#2E5090", penwidth="2")
dot.node("I2", "Fit maximum-likelihood\nSpatial Error model (ML_Error)")
dot.node("I3", "Report both specifications;\nfavor theoretically-motivated\nmodel")

dot.edge("A", "B")
dot.edge("B", "C")
dot.edge("C", "D")
dot.edge("D", "E1", label="  No")
dot.edge("D", "F", label="  Yes")
dot.edge("F", "G")
dot.edge("G", "H1", label="  Lag")
dot.edge("G", "H2", label="  Error")
dot.edge("G", "H3", label="  Ambiguous")
dot.edge("H1", "I1")
dot.edge("H2", "I2")
dot.edge("H3", "I3")

dot.render("fig_spatial_decision", cleanup=True)
print("saved")
