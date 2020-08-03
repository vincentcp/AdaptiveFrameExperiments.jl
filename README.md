# AdaptiveFrameExperiments.jl
A julia package containing the code to create the figures from <a href="https://arxiv.org/abs/2004.11317"><i>On the adaptive spectral approximation of functions using redundant sets and frames</i></a>
by Vincent Coppé and Daan Huybrechs.

Run the notebooks to create the figures yourself, its prerequisit is a Julia 1.3 kernel.
Two notebooks are presented that create similar figures. The notebook <strong>adaptivity figures, no AZ</strong> uses truncated singular value decompositions to solve all systems. This is computational heavy, but the most general approach. The notebook  <strong>adaptivity figures</strong> uses the AZ algorithm to speed up the experiments. More detail on the AZ algorithm is provided in <a href="https://arxiv.org/abs/1912.03648 "><i>The AZ algorithm for least squares systems with a known incomplete generalized inverse</i></a>. To appear in SIAM Journal on Matrix Analysis and Applications.
