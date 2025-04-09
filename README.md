# [Diffuse Like an Egyptian](https://www.youtube.com/watch?v=Cv6tuzHUuuk)

This code computes the MSD of a LAMMPS trajectory. It also provided a few C libraries for 


Mean Squared Displacement is:

$$MSD(\Delta t) = \frac{1}{N} \sum_{i=1}^{N} \left[ \vec{x}(t+\Delta t) - \vec{x}(t)\right]$$

To improve our esteem, we can average over time:

$$\langle MSD(\Delta t) \rangle_{t} = \frac{1}{N_T} \sum_{t=t_0}^{t_{N_T}}\left( \frac{1}{N} \sum_{i=1}^{N} \left[ \vec{x}(t+\Delta t) - \vec{x}(t)\right] \right)$$


