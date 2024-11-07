# Physics 427 Homework 9

__Due 11:59pm Wednesday 11/13/2024__

In addition to this homework, remember Project 2 is due in 2 weeks. If you haven't started working on it, please do! 

## 1. Time-dependent Schrödinger equation in 1D (25 points)

In this problem, we will try to solve the time-dependent Schrödinger equation in 1D in an infinite potential well. The 1D Schrödinger equation is:

$$
i\hbar \frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \psi}{\partial x^2} + V(x)\psi
$$

where $\psi(x, t)$ is the wavefunction, $V(x)$ is the potential, and $m$ is the mass of the particle. Using a procedure similar to HW6, we can make this equation dimensionless:

$$
i\frac{\partial \psi}{\partial t} = -\frac{1}{2}\frac{\partial^2 \psi}{\partial x^2} + V(x)\psi
$$

For the purpose of this homework, we will simply take $V = 0$ in our finite domain. Consider a domain of $x\in [-10, 10]$. There is an infinite potential barrier at the boundaries of the domain, therefore the boundary conditions are $\psi(-10, t) = \psi(10, t) = 0$. This boundary condition is automatically built-in in the numerical method, so you don't need to do anything special to implement it. Use a grid resolution of $N_x = 1000$.

We start with a Gaussian wave packet that is centered at $x = 0$ and traveling to the right at speed $v$:

$$
\psi(x, 0) = \exp\left(-\frac{x^2}{2\sigma^2} + ivx\right)
$$

where $\sigma$ is the width of the wave packet. Set $\sigma = 0.5$ and $v = 10.0$.

We will use the Crank-Nicolson method to solve this equation. The Crank-Nicolson method utilizes the following discretization of the Schrödinger equation:

$$
\left(1 + \frac{1}{2}i\hat{H}\Delta t\right)\psi^{n+1}_j = \left(1 - \frac{1}{2}i\hat{H}\Delta t\right)\psi^n_j,
$$

where $\hat{H}$ is the Hamiltonian operator, $\hat{H} = -\partial_x^2/2 + V$, and $\psi^n_j$ is the value of $\psi$ at the $j$-th grid point at time step $n$. I'm again forced to use $j$ for index since $i$ is reserved for the imaginary unit. The Hamiltonian operator in our case is simply a 2nd order spatial derivative, and it's approximated by a central difference:

$$
\hat{H}\psi^n_j = -\frac{1}{2}\frac{\psi^n_{j+1} - 2\psi^n_j + \psi^n_{j-1}}{\Delta x^2}.
$$

The Crank-Nicolsom method is an implicit method, and it requires solving a linear system at each time step. Fortunately, the linear system is tri-diagonal, and we can use the LU solver you wrote in Homework 7 to solve it.

Write a C++ file `problem1.cpp` to solve the time-dependent Schrödinger equation using the method above. Use $\Delta t = 10^{-3}$ and simulate the wave function for a total of $t = 1.0$. Create an output every time $t$ elapses by $0.01$, and write the values of $\psi$ in a csv file. Since the wave function is complex, you will need to include the header file `<complex>` and use `std::vector<std::complex<double>>` instead of `std::vector<double>` as the type of your arrays. You also need to use the `tri_diagonal<std::complex<double>>` version of the class you implemented before. To help facilitate implementing the initial condition, you can include a line before the `main` function:
```cpp
using namespace std::complex_literals;
```
to use the `i` literal for complex numbers. When assigning the initial condition, you can then write something like the following:
```cpp
psi[j] = std::exp(-0.5 * (x * x) / (sigma * sigma) + 1.0i * v * x);
```
Here `1.0i` means the imaginary unit $i$. This helps a lot with the readability of the code.

When writing output, remember to write both the real and imaginary parts of $\psi$ to the csv file. You can do it with something like the following:
```cpp
for (int i = 0; i < Nx; i++) {
  double x = x0 + i * dx;
  output_file << x << "," << std::real(psi[i]) << "," << std::imag(psi[i]) << std::endl;
}
```
Plot the outputs and create a movie of the evolution of $\psi$ in time. Remember to plot $|\psi|^2 = (\mathrm{Re}\ \psi)^2 + (\mathrm{Im}\ \psi)^2$ instead of $\psi$, since $|\psi|^2$ is the probability of finding the particle at position $x$, while $\psi$ itself doesn't have a clear physical meaning. Commit the resulting `problem1.mp4` to the repository. DO NOT COMMIT CSV FILES, FIGURES, OR BINARY FILES TO THE REPOSITORY!

## 2. Solving the Poisson Equation in 2D (25 points)

In this problem, we will solve the Poisson equation in 2D:

$$
\nabla^2 \phi = \frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = -\rho(x, y)
$$

Your goal is to implement the Successive Overrelaxation (SOR) method to solve this differential equation with a given boundary condition and source term $\rho$. SOR is an iterative method to solve linear systems of equations that result from the discretization of a PDE. For the Poisson equation, the SOR scheme is given by:

$$
\phi_{i,j}^{n+1} = \phi_{i,j}^n + \frac{\omega}{4} \left(\phi_{i+1,j}^n + \phi_{i-1,j}^{n} + \phi_{i,j+1}^n + \phi_{i,j-1}^{n} + \rho_{i,j}\Delta^2 - 4 \phi_{i,j}^n\right)
$$

where we have assumed that $\Delta x = \Delta y = \Delta$. $\omega$ is the overrelaxation parameter, and the scheme is stable only when $0 < \omega < 2$. When $\omega = 1$, this scheme is called Jacobi's method. You can also use the Gauss-Seidel construction where you alternate between $\phi^n$ and $\phi^{n+1}$ in the bracket.

In this homework, we will adopt a box of size $L_x \times L_y = 1.0 \times 1.0$ with homogeneous Neumann boundary conditions on all sides:

$$
\begin{align}
\frac{\partial \phi}{\partial x} &= 0\quad \mathrm{when}\quad x = 0\mathrm{\ or\ }1 \\
\frac{\partial \phi}{\partial y} &= 0\quad \mathrm{when}\quad y = 0\mathrm{\ or\ }1
\end{align}
$$

Use ghost cells to enforce the boundary conditions (see Lecture 19.5). I suggest a resolution of $1024\times 1024$. The source term $\rho$ is given by:

$$
\rho(x, y) = \begin{cases}
1 & \mathrm{if}\quad 0.3 < x < 0.7 \; \mathrm{\ and\ }\; 0.65 < y < 0.7 \\
-1 & \mathrm{if}\quad 0.3 < x < 0.7 \; \mathrm{\ and\ }\; 0.3 < y < 0.35 \\
0 & \mathrm{otherwise}
\end{cases}
$$

This charge distribution describes two uniformly charged plates, one positive and one negative, separated by some distance. The solution should be similar to the electric potential surrounding a capacitor.

In `problem2.cpp`, implement the SOR algorithm and solve the Poisson equation with the given boundary conditions and source term. Experiment with your overrelaxation parameter $\omega$ so that it converges as quickly as possible.

To measure convergence, define the residual as:

$$
R = \sqrt{\sum_{i,j}\left(\phi_{i+1,j}^n + \phi_{i-1,j}^{n+1} + \phi_{i,j+1}^n + \phi_{i,j-1}^{n+1} - \rho_{i,j}\Delta^2 - 4 \phi_{i,j}^n\right)^2}
$$

Iterate until $R < 10^{-8}$. If your implementation does not converge in 10,000 iterations, you might want to adjust your overrelaxation parameter $\omega$. Your code should print out the step number and the residual at each step, in the following format:
    
    0 [R at step 0]
    1 [R at step 1]
    2 [R at step 2]
    ...

Save your program output to `problem2.txt` using the shell redirection:

    ./a.out > problem2.txt

At the end of the calculation, write the values of the solution $\phi$ to a csv file. Plot the value of $\phi$ in python using `matplotlib`. In the same plot, visualize the electric field using [`streamplot`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html). The result should look like the plot shown in Lecture 21. If you see that the electric field is not becoming parallel to the boundaries, but pointing into them, then it means you have implemented the wrong boundary condition. You can compute $E_x$ and $E_y$ in python using the following:

```python
Ex = -(np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)) / (2.0 * dx)
Ey = -(np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)) / (2.0 * dy)
```

Remember to label your $x$ and $y$ axes correctly. Save the plot as `problem2.png` and commit it to the homework repository. DO NOT COMMIT THE CSV OUTPUT.


