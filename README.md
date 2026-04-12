# llopt-hawt
Design problem for Horizontal Axis Wind Turbine blades, using Lifting Line theory and Induction Factors

-----------------------

# References
1. Mishkevich V. G., A New Approach to Lifting Line Theory: Hub and Duct Effects, Journal of Ship Research, Vol. 50, No. 2, June 2006, pp. 138-146

2. Rand O., Rosen A. Efficient Method for Calculating the Axial Velocities Induced Along Rotating Blades by Trailing Helical Vortices. Journal of Aircraft, Vol. 21, June 1984

3. Chiu Y. D., Peters D. A. Numerical Solutions of Induced Velocities by Semi-Infinite Tip Vortex Lines, Journal of Aircraft, Vol. 25, Aug. 1988

4. Dantchev O. D., Design of HAWT Blades using Lifting Line Theory, M.Sc. Thesis, TU - Sofia, 1996

------------------------

Z = number of blades

$ \Gamma $ = bound circulation

G = dimensionless bound circulation

$ x, r, \theta $ = cylindrical coordinates 

V = wind speed

$ u_{a}, u_{t} $ = induced axial and tangental velocity, respectively 

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

### `delta`
Per-station aerodynamic data: $r/R$, $C_l$, $C_d$, $r$, $c/R$, camber, thickness, twist (deg).

### `prof`
Raw airfoil coordinates at each station: chordwise position and upper/lower surface ordinates scaled to physical dimensions.

### `g`
Dimensionless bound circulation $G(x)$ and its derivative across the span.

