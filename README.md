# llopt-hawt
Design problem for Horizontal Axis Wind Turbine blades, using Lifting Line theory and Induction Factors

-----------------------

30th anniversary edition 
* ported to GNU Fortran (from HP-UX Fortran and Fortran PowerStation 4)
* resolved the dependency on IMSL (a commercial library)
* started a post-processing/visualization script
* adding some of the original references 


# References
1. Mishkevich V. G., A New Approach to Lifting Line Theory: Hub and Duct Effects, Journal of Ship Research, Vol. 50, No. 2, June 2006, pp. 138-146

2. Rand O., Rosen A. Efficient Method for Calculating the Axial Velocities Induced Along Rotating Blades by Trailing Helical Vortices. Journal of Aircraft, Vol. 21, June 1984

3. Chiu Y. D., Peters D. A. Numerical Solutions of Induced Velocities by Semi-Infinite Tip Vortex Lines, Journal of Aircraft, Vol. 25, Aug. 1988

4. Dantchev O. D., Design of HAWT Blades using Lifting Line Theory, M.Sc. Thesis, TU - Sofia, 1996

------------------------

| Symbol | Description |
|---|---|
| $D$ | impeller diameter |
| $R$ | impeller radius |
| $Z$ | number of impeller blades |
| $n$ | rotational frequency, s⁻¹ |
| $\omega$ | angular velocity of rotation |
| $V$ | undisturbed flow velocity |
| $\bar{r}_H$ | dimensionless hub radius, $\bar{r}_H = r_H / R$ |
| $\Gamma$ | circulation of the velocity around the profile of the cylindrical section |
| $G$ | dimensionless circulation $G = \Gamma / (DV)$ |
| $u_a, u_t, u_r$ | induced velocities |
| $i_a, i_t$ | induction factors |
| $J$ | advance ratio, $J = V / nD$ |
| $H$ | helix pitch, $H = 2\pi V / \omega$ |
| $\lambda_I$ | inductive advance coefficient, $\lambda_I = \dfrac{V + u_a}{r\omega - u_t}$ (upper sign for propeller) |
| $\lambda_t$ | advance coefficient $\lambda_t = J / \pi$ |
| $T$ | thrust |
| $Q$ | torque |
| $K_T, K_Q, K_P$ | coefficients of: thrust, torque and power |
| $C_T, C_P, C_Q$ | load coefficients for: thrust, power, torque |
| $F_P$ or $A_0$ | turbine (propeller) disk area |
| $B$ | hydrodynamic characteristic of the section, $B = cC_L / 2$ |
| $c$ | chord length of the airfoil section |
| $\beta$ | advance angle, $\beta = \arctan(\lambda_t / \bar{r})$ |
| $\beta_I$ | inductive advance angle, $\beta_I = \arctan(\lambda_I / \bar{r})$ |
| $\varphi$ | setting angle |
| $\alpha$ | geometric angle of attack |
| $\alpha_0$ | zero-lift angle of attack |
| $\alpha_a$ | aerodynamic angle of attack |
| $C_L$ | lift coefficient of the airfoil section |
| $C_D$ | drag coefficient of the airfoil section |
| $\bar{C}$ | relative camber of the airfoil section |
| $\bar{\delta}$ | relative thickness of the airfoil section |
| $\varepsilon$ | inverse quality of the profile, $\varepsilon = C_D / C_L$ |
| $Rn_S$ | Reynolds number, $V_R c / \nu$ |
| $x$ | $x = \dfrac{\bar{r} - \bar{r}_H}{1 - \bar{r}_H} \cdot 2 - 1$ |
| $\chi_j$ | angular displacement of the $j^{th}$ blade |
| $\chi_j$ | $\chi_j = \dfrac{2\pi j}{Z}, \quad j = 0, 1, \ldots, Z-1$ |
| $\rho$ | fluid density |
| $\nu$ | kinematic viscosity of the fluid |
| ЛРК | impeller blade |
| ВТХО or HAWT | horizontal axis wind turbine |

------------------------

## Output Files

### `planform.csv`
Spanwise planform distribution table for aerodynamic solvers. Each row is a cylindrical section along the dimensionless radius $r/R$.

| Column | Symbol | Description |
|--------|--------|-------------|
| `r/R` | $r/R$ | Dimensionless radial station |
| `r [m]` | $r$ | Radial distance from hub center (meters) |
| `chord [m]` | $c$ | Section chord length, leading edge to trailing edge (meters) |
| `twist [deg]` | $\beta$ | Geometric twist angle of the section relative to the plane of rotation (degrees) |
| `thickness` | $\delta$ | Maximum thickness-to-chord ratio ($t/c$) |
| `camber` | $\delta_c$ | Maximum camber-to-chord ratio |
| `Cl` | $C_l$ | Design section lift coefficient |
| `Cd` | $C_d$ | Design section drag coefficient |

Profiles are custom blends of NACA 66 thickness distribution and NACA $a=0.8$ mean line, scaled per station by the thickness and camber values.

### `sections.dat`
Gnuplot-compatible file with all section profiles in normalized coordinates ($x/c$, $y/c$). Each spanwise station is a separate data block delimited by blank lines, with a comment line showing $r/R$. Plot all sections with:
```
plot 'sections.dat' with lines
```

### `section_XX.dat`
One file per spanwise station in Selig airfoil format, compatible with QBlade, XFOIL, XFLR5, and OpenVSP. The first line is a header with $r/R$, $t/c$, and camber; remaining lines are $x/c$, $y/c$ coordinates running from upper trailing edge around the leading edge to lower trailing edge. Reference the station files from `planform.csv` using the `Airfoil` column convention of the target solver.

### `blade.plt`
Twisted, physical-scale section profiles for visualizing the actual blade shape. Each station's coordinates are rotated by the local twist angle $\beta$ and scaled to the physical chord length. Gnuplot-compatible with blank-line-separated blocks. Plot with:
```
set size ratio -1; plot 'blade.plt' with lines notitle
```

### `delta`
Per-station aerodynamic data: $r/R$, $C_l$, $C_d$, $r$, $c/R$, camber, thickness, twist (deg).

### `prof`
Raw airfoil coordinates at each station: chordwise position and upper/lower surface ordinates scaled to physical dimensions.

### `g`
Dimensionless bound circulation $G(x)$ and its derivative across the span.

