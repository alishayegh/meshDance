# meshDance

## Description
A code to read a scalar field on a particular mesh, move the mesh nodes while preserving topology, and calculate the field on the new mesh.

## Theory
Writing the integral transport equation for a variable $T$ leads to

$$
\frac{d}{dt}\int_{AR}\rho T \mathrm{d}V = 0 + \int_{AR}\rho T n_i (u_i - {u_b}_i)\mathrm{d}S
$$

where $AR$ is an arbitrary region, $u_i$ is the material velocity, ${u_b}_i$ is the boundary velocity for the arbitrary region, and $0$ on the RHS stands for
$\int_{AR} \rho \frac{\mathrm{d}T}{\mathrm{d}t}\mathrm{d}V$, which is $0$, because we are solving for a fixed time, i.e., $t$ does not exist in our problem.

## Dependencies
foam-extend-4.0

## Signature
Maalik, Maxwell corner, Cavendish Lab, Cambridge, 040924  
ali@tensorfields.com
