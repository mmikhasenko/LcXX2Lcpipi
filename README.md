# Spin-Averaged Distribution for $\Lambda_c^{**+} \to \Lambda_c^+ \pi^+ \pi^-$ Decays

This repository contains a Julia notebook for computing spin-averaged angular functions for three-body decays of excited $\Lambda_c^{**+}$ states, exploring various quantum number hypotheses.

## Key Features:
- **Spin Factors Evaluation:** Based on conventions from the Dalitz plot decomposition for three-body decays.
- **Symbolic Implementation:** Utilizes [ThreeBodyDecay.jl](https://github.com/mmikhasenko/ThreeBodyDecay.jl) and [SymbolicThreeBodyDecay.jl](https://github.com/mmikhasenko/SymbolicThreeBodyDecays.jl).
- **Intensity Distribution:** Decay intensity as a bilinear form on isobars lineshapes.
- **Isobar Resonances:** Six subchannel resonances included.
- **LS Couplings and Helicity Treatment:** Detailed mathematical descriptions and symbolic calculations.

## Installation
1. **Install Julia:** Download and install Julia from the [official website](https://julialang.org/downloads/).
2. **Install Pluto:** Open Julia and install Pluto by entering the following in the Julia REPL:
   ```julia
   ] add Pluto
   ```

## Rerunning the Notebook

> [!NOTE]
> The computation takes about 30'

- **Run with Julia:**
  ```bash
  julia notebooks/symbolic_expressions.jl
  ```
- **Run with Pluto:**
  1. Start Pluto in Julia:
     ```julia
     using Pluto
     Pluto.run()
     ```
  2. Open the notebook from the Pluto interface.


## Output File Structure

The output of the notebook is a structured JSON file containing an array of symbolic expressions. These expressions are exported in multiple formats to facilitate various uses. The structure of the JSON file is as follows:

1. **["lineshapes"]:** A collection of expressions representing the lineshapes of intermediate resonances, provided in SymPy format.

2. **["intensity_latex"]:** The decay intensity distribution is presented in LaTeX format for easy interpretation and rendering in documents or presentations.

3. **["spin_averaged_crossings"]:**
    - **["latex"]:** The spin-averaged crossings are formatted in LaTeX, offering a mathematical representation that is clear and ready for publication.
    - **["sympy"]:** Each element of the spin-averaged matrix is represented in SymPy format, allowing for easy integration into symbolic computation environments.
    - **["ccode"]:** The elements are also available in C code format, facilitating their use in computational simulations or numerical analysis.

The json file is stored in the `data/` folder.
The latest computation is labeled with the highest version number, `data/LcXX2Lcpipi.v[?].json`.