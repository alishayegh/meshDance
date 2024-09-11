# meshDance

## Description
A code to read a scalar field on a particular mesh, move the mesh nodes while preserving topology, and calculate the field on the new mesh.

## Theory
Writing the integral transport equation for a variable $T$ leads to

$$
\frac{d}{dt}\int_{AR}\rho T \mathrm{d}V = \int_{AR} \rho \frac{\mathrm{d}T}{\mathrm{d}t}\mathrm{d}V + \int_{AR}\rho T n_i ({u_b}_i - {u}_i)\mathrm{d}S
$$

where $AR$ is an arbitrary region, $u_i$ is the material velocity, and ${u_{b_i}}$ is the boundary velocity for the arbitrary region. 
<!--:#and $0$ on the RHS stands for
$\int_{AR} \rho \frac{\mathrm{d}T}{\mathrm{d}t}\mathrm{d}V$, which is $0$, because we are solving for a fixed time, i.e., $t$ does not exist in our problem. -->
The discretised form reads  

$$
\frac{{(TV)}^N-{(TV)}^o}{\Delta t} = (\frac{\mathrm{d}T}{\mathrm{d}t})_PV_P + \sum_f T S_i ({u_b}_i - {u}_i)
$$

where superscripts $\mathrm{N}$ and $\mathrm{O}$ stand for new and old. Since we are looking at a snapshot (fixed time) while we move the mesh, we have $(\frac{\mathrm{d}T}{\mathrm{d}t})_P = 0$ (no field generation or destruction) and ${u}_i = \mathbf{0}$
(again, since it is a fixed point in time). In other words, `time' variable in the mesh motion is not the physical time.
Therefore, the implicit form for T reads  

$$
\frac{{(TV)}^N-{(TV)}^o}{\Delta t} = \sum_f T S_i{u_b}_i
$$

and in OpenFOAM notation,

```c++
fvm::ddt(T) == fvm::div(mesh.phi(), T) 
```
Also, solving this explicitly for $T^N$ reads  

$$
{{T}^N} = \frac{1}{V^N}[{(TV)}^o + [\sum_f T S_i {u_b}_i]{\Delta t}]
$$

## Example
![](https://github.com/alishayegh/meshDance/blob/master/figures/ALE_cavity_meshOnly.gif)  

![](https://github.com/alishayegh/meshDance/blob/master/figures/ALE_cavity_fieldOnly.gif)

## Dependencies
foam-extend-4.0

## Signature
Maalik, Maxwell Corner, Cavendish Lab, Cambridge, 040924  
ali@tensorfields.com
