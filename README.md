# meshDance

## Description
A code to read a scalar field on a particular mesh, move the mesh nodes while preserving topology, and calculate the field on the new mesh.

## Theory
Writing the integral transport equation for a variable $T$ leads to

$$
\frac{d}{dt}\int_{AR}\rho T \mathrm{d}V = \int_{AR} \rho \frac{\mathrm{d}T}{\mathrm{d}t}\mathrm{d}V + \int_{AR}\rho T n_i ({u_b}_i - {u_b})\mathrm{d}S
$$

where $AR$ is an arbitrary region, $u_i$ is the material velocity, and ${u_{b_i}}$ is the boundary velocity for the arbitrary region. 
<!--:#and $0$ on the RHS stands for
$\int_{AR} \rho \frac{\mathrm{d}T}{\mathrm{d}t}\mathrm{d}V$, which is $0$, because we are solving for a fixed time, i.e., $t$ does not exist in our problem. -->
The discretised form reads  

$$
\frac{{(TV)}^N-{(TV)}^o}{\Delta t} = (\frac{\mathrm{d}T}{\mathrm{d}t})_PV_P + \sum_f T S_i ({u_b}_i - {u}_i)
$$

where superscripts $\mathrm{N}$ and $\mathrm{O}$ stand for new and old. Solving this explicitly for $T^N$ reads  

$$
{{T}^N} = \frac{1}{V^N}[{(TV)}^o + [(\frac{\mathrm{d}T}{\mathrm{d}t})_PV_P + \sum_f T S_i ({u_b}_i - {u}_i)]{\Delta t}]
$$

The implicit form reads  

```c++
fvm::ddt(T) == fvc::ddt(T) + fvm::div(meshPhi, T) - fvm::div(phi, T)
```

## Example
![](https://github.com/alishayegh/meshDance/blob/master/figures/ALE_cavity_fieldOnly.gif)  

![](https://github.com/alishayegh/meshDance/blob/master/figures/ALE_cavity_fieldOnly.gif)

## Dependencies
foam-extend-4.0

## Signature
Maalik, Maxwell Corner, Cavendish Lab, Cambridge, 040924  
ali@tensorfields.com
