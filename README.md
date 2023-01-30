# Computational Physics

Code to accompany lectures.
See the `Demos` directory for code from in-class demos.

## Week 1: Projectile motion using the Euler method

- `proj_euler.m`.
  - Solves the projectile motion problem using the Euler method.

## Week 2: Projectile motion using the Euler method

## Week 3: A General Form of ODEs and the Simple Pendulum

The exponential problem, dx/dt = x:
- `exp_euler.m`
  - Solve the exponential problem using Euler.
- `exp_rk4.m`
  - Solve the exponential problem using RK4.
- `rhs_exp.m`
  - "My first right-hand side" function :baby:
- `exp_battle`
  - Animate a step-by-step comparison between Euler and RK4.

The simple pendulum:
- `pend_rk4.m`
  - Animates the RK4 solution to the simple pendulum.
- `rhs_pend.m`
  - The RHS function for the simple pendulum.

## Week 4: 1D Diffusion

- `diffusion_ftcs.m`
  - Solves the 1D diffusion equation for an initial spike profile with Dirichlet conditions using FTCS, in a matrix formulation.

## Week 5: 1D Advection

- `advection.m`
  - Solves the advection problem for a Gaussian pulse intial condition using the FTCS and Lax methods.
