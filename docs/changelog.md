# 2.0.0

- Breaks backwards compatibility, but hopefully a useful starting point for development of systems with cross sections that vary along their length.
- Changes the mechanical state vector from `[Xp; dXdt; M; δstart; δend; Lfilm_start; Lfilm_end]` to `[ζ_L; ζ_R; Q_B; M; Af_L; Af_R; Lfilm_start; Lfilm_end]`.
    - `Xp` and `dXdt` previously stored pairs of interface positions and velocities. `ζ_L` and `ζ_R` now store volume coordinates, `Q_B` stores bulk volume flow, and `Af_L` and `Af_R` store film areas.
        - The relationship between these two state vectors for a uniform tube is given quite explicitly in the `reparse_u_v1_as_u_v2()` function.
    - Temperature profiles and their separators are unchanged.
- Adds the beginning of a `VariableTube` and `VariablePhpSystem`, although propagation hasn't been implemented yet.
- Adds mappings between arclength and volume coordinates for variable cross sections.

- I've done my best to keep functionality the same (even in cases where there might have been bugs).
    - Unit tests have been tweaked to use the new state vector while testing the same behavior and physics.
    - Despite that, the exact state trajectories from before/after this change converge only in the limit of small timesteps. This is because of a nonlinear change of variables from $\delta_f$ to $A_f$ (for a uniform tube, all the other variable changes are linear). RK4 produces a slightly different trajectory under this transformation, which can change when callbacks fire.

- Known issues to be addressed in future versions:
    - Initialization calculates charge ratio differently than literature (and the charge ratio function does).
    - Potential indexing issue in `openloopdynamics`.
    - Readability and optimization of `Thermomodel.jl` and the callbacks. Current integration makes minimal changes instead of clean code.
    - Remeshing should keep a minimum of two points per slug
    - I want to revisit boiling model