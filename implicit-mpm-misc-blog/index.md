# Implicit MPM Misc


# Some Other Things for Implicit MPM

In the last blog, we talked about some basic things for implicit MPM. But there are also some details we should mention in this continued blog.

## Hessian Matrix

In the last blog, we defined $E(\mathbf{u}) = \frac{1}{2}(\mathbf{u}-\mathbf{u^\*})\mathbf{M}(\mathbf{u}-\mathbf{u^\*}) + \Phi(\mathbf{x}^n+\Delta t \mathbf{u})$ as the energy for a system, our goal is to minimize it. We know the gradient of it is :

$$
\nabla_\mathbf{u}E = \mathbf{M}(\mathbf{u}-\mathbf{u}^\*) + \Delta t \frac{\partial \Phi}{\partial\mathbf{x}}|_{(\mathbf{x}^n+\Delta t \mathbf{u})} = \mathbf{M}(\mathbf{u}-\mathbf{u}^\*) - \Delta t \mathbf{f}(\mathbf{x}^n+\Delta t \mathbf{u})
$$

However, in some optimization methods, we have to calculate the Hessian matrix $\mathbf{H}$ for $E$. $\mathbf{H}$ is a really big but sparse matrix in this problem, so sometimes we can save it in the memory using designed data structures. We wouldn't talk about the details of storing a sparse matrix, but it's necessary to show how to calculate it. Again, compute the derivative $\nabla^2 E$ : 

$$
\nabla^2_\mathbf{u}E = \mathbf{M} + \Delta t^2 \frac{\partial^2\Phi}{(\partial\mathbf{x})^2}(\mathbf{x}^n+\Delta t \mathbf{u})
$$

How to compute the $\frac{\partial^2\Phi}{(\partial\mathbf{x})^2}$ ? You can find a lot of articles proving how to compute $\frac{\partial\Phi}{\partial\mathbf{x}} = -\mathbf{f}$ in MPM :
$$
\frac{\partial\Phi}{\partial\mathbf{x}\_i} = \sum_p V_p^0 \mathbf{P}(\mathbf{F}\')\mathbf{F}^T \nabla\omega_{ip}
$$

The $\mathbf{x}_i$ here means we are computing the entries $\\{i\times3, i\times3+1, i\times3+2\\}$ of $\frac{\partial\Phi}{\partial\mathbf{x}}$ . $\mathbf{F}\'$ here means "new deformation gradient", which is different from $\mathbf{F}$. Use the similar techniques to compute the derivative of this equation, we have:

$$
\frac{\partial^2\Phi}{\partial\mathbf{x}\_i\partial\mathbf{x}\_j} = \sum_p V_p^0 \frac{\partial\mathbf{P}}{\partial\mathbf{F}}(\mathbf{F}\')(\mathbf{F}^T\nabla\omega_{jp})(\mathbf{F}^T\nabla\omega_{ip})^T
$$

For each particle, it will generate a $81\times81$ hessian matrix, which is a submatrix of the whole $\nabla^2_\mathbf{u}E$, because it is connected to $27$ grids in 3D simulation. From the equation above, we know how to compute a $3\times3$ submatrix of that $81\times 81$ matrix because only 2 grids are evolved. Finally, the code for computing hessian matrix should be like:

``` python
for p in particles:
    for g1 in 27 grids:
        for g2 in 27 grids:
            H = compute_small_hessian(p, g1, g2) # 3 by 3 matrix
            total_H(g1*3 + i, g2*3 + j) += H(i, j) # add it to the whole hessian matrix
```

## Line Search & Projected Newton 
In Newton Method, we solve $\mathbf{H}\Delta\mathbf{x} = -\nabla E$ in each step and use $\Delta\mathbf{x}$ to update. But there are two problems: 
* This method is not guaranteed to converge;
* Solved $\Delta\mathbf{x}$ is not guaranteed to reduce the energy $E$ in each step.

To eliminate these two problems, we should use *Line Search* method to make it converge and *Projected Newton* to make sure energy is getting lower in each step. For details of these two methods, you can read Chapter 3 of this online book [Physics-Based Simulation](https://phys-sim-book.github.io/).
