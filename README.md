# Numerical Analysis of Differential Equations using Matlab
Group 25

**Final Paper:** 

## Abstract
In the first part, different Runge-Kutta methods are used to model an RC circuit. In the second part, said methods are applied to RLC circuits. The third and last part implements the Relaxation method in order to model solutions to the Laplace equation.

During part one, we observed a low-pass filter usually used in NAB filters. The circuit cuts off any signals above 1.592kHz. The first Exercise was to use three methods; Heun’s, midpoint and custom method. Once we had the numerically approximated solutions to the ODE, we compared our results to the exact solutions for error analysis, expressing the error between the two as a function of the chosen step-size and thus, we showed that the order of the global truncation error of the numerical solution is O(h2). We found that the filter behaved as expected for input signals of low frequency but for higher frequency inputs it acted as an integrator. For each input waveform, a transient and steady state output could be observed.

In part two we came across an RLC circuit acting as an band-pass filter. The purpose of this part was to approximate the output signal using classic 4th order Runge-Kutta method, for a variety of input signals.

In the third part, we applied the Relaxation method to numerically solve Laplace’s equation on the Cartesian plane and model it. We analysed limitations of the relaxation method, and further explored how different resolutions, maximum errors and boundary conditions impact the performance of our implementation, both in terms of its accuracy and complexity. Furthermore, we analyzed how Successive Over-relaxation can improve the performance of our algorithm and how the optimal relaxation factor depends on the boundary conditions of the system. We found that a relaxation factor of 1.94 works consistently well for all systems, and that the optimal relaxation factor is lower for more complex boundary conditions.
