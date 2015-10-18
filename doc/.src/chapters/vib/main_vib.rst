.. Automatically generated Sphinx-extended reStructuredText file from DocOnce source
   (https://github.com/hplgit/doconce/)

.. Document title:

Finite difference methods for vibration problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:Authors: Hans Petter Langtangen
:Date: Oct 17, 2015

.. Externaldocuments: /home/hpl/vc/decay-book/doc/src/book/book

.. Note: **PRELIMINARY VERSION** (expect typos)

.. !split

.. 2DO:

.. _undamped -> _simple everywhere

.. Long time integration by adaptive RK: will that improve the

.. phase error? Do experiments where we measure the wavelength

.. and plot it as function of time. Can we vectorize the

.. max/min pt computation?

Vibration problems lead to differential equations with solutions that
oscillate in time, typically in a damped or undamped sinusoidal
fashion.  Such solutions put certain demands on the numerical methods
compared to other phenomena whose solutions are monotone or very smooth.
Both the frequency and amplitude of the oscillations need to be
accurately handled by the numerical schemes. Most of the reasoning and
specific building blocks introduced in the forthcoming text can be
reused to construct sound methods for partial differential equations
of wave nature in multiple spatial dimensions.

[**hpl 1**: Need to discuss errors also for the damped and nonlinear models. At least the frequency errors must be illustrated here as well and investigated numerically, either in text or exercises.]

.. _vib:model1:

Finite difference discretization
================================

Many of the numerical challenges faced when computing oscillatory
solutions to ODEs and PDEs can be captured by the very simple ODE
:math:`u^{\prime\prime} + u =0`. This ODE is thus chosen as our starting
point for method development, implementation, and analysis.

A basic model for vibrations
----------------------------

.. index:: vibration ODE

.. index:: oscillations

.. index:: mechanical vibrations

A system that vibrates without damping and external forcing
can be described by the ODE problem

.. math::
   :label: vib:ode1
        
        u^{\prime\prime} + \omega^2u = 0,\quad u(0)=I,\ u^{\prime}(0)=0,\ t\in (0,T]
        {\thinspace .}
        
        

Here, :math:`\omega` and :math:`I` are given constants.
The exact solution of :eq:`vib:ode1` is

.. index:: period (of oscillations)

.. index:: frequency (of oscillations)

.. index:: Hz (unit)

.. math::
   :label: vib:ode1:uex
        
        u(t) = I\cos (\omega t)
        {\thinspace .}
        
        

That is, :math:`u` oscillates with constant amplitude :math:`I` and
angular frequency :math:`\omega`.
The corresponding period of oscillations (i.e., the time between two
neighboring peaks in the cosine function) is :math:`P=2\pi/\omega`.
The number of periods per second
is :math:`f=\omega/(2\pi)` and measured in the unit Hz.
Both :math:`f` and :math:`\omega` are referred to as frequency, but :math:`\omega`
is more precisely named *angular frequency*, measured in rad/s.

In vibrating mechanical systems modeled by :eq:`vib:ode1`, :math:`u(t)`
very often represents a position or a displacement of a particular
point in the system. The derivative :math:`u^{\prime}(t)` then has the
interpretation of velocity, and :math:`u^{\prime\prime}(t)` is the associated
acceleration.  The model :eq:`vib:ode1` is not only
applicable to vibrating mechanical systems, but also to oscillations
in electrical circuits.

.. _vib:ode1:fdm:

A centered finite difference scheme
-----------------------------------

To formulate a finite difference method for the model
problem  :eq:`vib:ode1` we follow the `four steps <http://tinyurl.com/opdfafk/pub/sphinx-decay/main_decay.html#the-forward-euler-scheme>`__ explained in [Ref1]_.

.. index::
   single: mesh; finite differences

.. index:: mesh function

Step 1: Discretizing the domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The domain is discretized by
introducing a uniformly partitioned time mesh.
The points in the mesh are :math:`t_n=n\Delta t`, :math:`n=0,1,\ldots,N_t`,
where :math:`\Delta t = T/N_t` is the constant length of the time steps.
We introduce a mesh function :math:`u^n` for :math:`n=0,1,\ldots,N_t`, which
approximates the exact solution at the mesh points. The mesh
function will be computed from algebraic equations derived from
the differential equation problem.

Step 2: Fulfilling the equation at discrete time points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ODE is to be satisfied at each mesh point:

.. math::
   :label: vib:ode1:step2
        
        u^{\prime\prime}(t_n) + \omega^2u(t_n) = 0,\quad n=1,\ldots,N_t
        {\thinspace .}
        
        

.. index:: centered difference

.. index::
   single: finite differences; centered

Step 3: Replacing derivatives by finite differences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The derivative :math:`u^{\prime\prime}(t_n)` is to be replaced by a finite
difference approximation. A common second-order accurate approximation
to the second-order derivative is

.. math::
   :label: vib:ode1:step3
        
        u^{\prime\prime}(t_n) \approx \frac{u^{n+1}-2u^n + u^{n-1}}{\Delta t^2}
        {\thinspace .}
        
        

Inserting :eq:`vib:ode1:step3` in :eq:`vib:ode1:step2`
yields

.. math::
   :label: vib:ode1:step3b
        
        \frac{u^{n+1}-2u^n + u^{n-1}}{\Delta t^2} = -\omega^2 u^n
        {\thinspace .}
        
        

We also need to replace the derivative in the initial condition by
a finite difference. Here we choose a centered difference, whose
accuracy is similar to the centered difference we used for :math:`u^{\prime\prime}`:

.. math::
   :label: vib:ode1:step3c
        
        \frac{u^1-u^{-1}}{2\Delta t} = 0
        
        {\thinspace .}
        

Step 4: Formulating a recursive algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To formulate the computational algorithm, we assume that we
have already computed :math:`u^{n-1}` and :math:`u^n` such that :math:`u^{n+1}` is the
unknown value, which we can readily solve for:

.. math::
   :label: vib:ode1:step4
        
        u^{n+1} = 2u^n - u^{n-1} - \Delta t^2\omega^2 u^n
        {\thinspace .}
        
        

The computational algorithm is simply to apply :eq:`vib:ode1:step4`
successively for :math:`n=1,2,\ldots,N_t-1`. This numerical scheme sometimes
goes under the name
Stormer's
method or `Verlet integration <http://en.wikipedia.org/wiki/Verlet_integration>`__.

Computing the first step
~~~~~~~~~~~~~~~~~~~~~~~~

We observe that :eq:`vib:ode1:step4` cannot be used for :math:`n=0` since
the computation of :math:`u^1` then involves the undefined value :math:`u^{-1}`
at :math:`t=-\Delta t`. The discretization of the initial condition
then comes to our rescue: :eq:`vib:ode1:step3c` implies :math:`u^{-1} = u^1`
and this relation can be combined with :eq:`vib:ode1:step4`
for :math:`n=1` to yield a value for :math:`u^1`:

.. math::
         u^1 = 2u^0 - u^{1} - \Delta t^2 \omega^2 u^0,

which reduces to

.. math::
   :label: vib:ode1:step4b
        
        u^1 = u^0 - \frac{1}{2} \Delta t^2 \omega^2 u^0
        {\thinspace .}
        
        

:ref:`vib:exer:step4b:alt` asks you to perform an alternative derivation
and also to generalize the initial condition to :math:`u^{\prime}(0)=V\neq 0`.

The computational algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The steps for solving :eq:`vib:ode1` becomes

 1. :math:`u^0=I`

 2. compute :math:`u^1` from :eq:`vib:ode1:step4b`

 3. for :math:`n=1,2,\ldots,N_t-1`:

   1. compute :math:`u^{n+1}` from :eq:`vib:ode1:step4`

The algorithm is more precisely expressed directly in Python:

.. code-block:: python

        t = linspace(0, T, Nt+1)  # mesh points in time
        dt = t[1] - t[0]          # constant time step
        u = zeros(Nt+1)           # solution
        
        u[0] = I
        u[1] = u[0] - 0.5*dt**2*w**2*u[0]
        for n in range(1, Nt):
            u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n]


.. admonition:: Remark on using ``w`` for :math:`\omega`

   In the code, we use ``w`` as the symbol for :math:`\omega`.
   The reason is that this author prefers ``w`` for readability
   and comparison with the mathematical :math:`\omega` instead of
   the full word ``omega`` as variable name.




Operator notation
~~~~~~~~~~~~~~~~~

We may write the scheme using a compact difference notation
(see also 
`examples <http://tinyurl.com/opdfafk/pub/sphinx-decay/main_decay.html#compact-operator-notation-for-finite-differences>`__ in [Ref1]_).
The difference :eq:`vib:ode1:step3` has the operator
notation :math:`[D_tD_t u]^n` such that we can write:

.. math::
   :label: vib:ode1:step4:op
        
        [D_tD_t u  + \omega^2 u = 0]^n
        {\thinspace .}
        
        

Note that :math:`[D_tD_t u]^n` means applying a central difference with step :math:`\Delta t/2` twice:

.. math::
         [D_t(D_t u)]^n = \frac{[D_t u]^{n+\frac{1}{2}} - [D_t u]^{n-\frac{1}{2}}}{\Delta t}

which is written out as

.. math::
        
        \frac{1}{\Delta t}\left(\frac{u^{n+1}-u^n}{\Delta t} - \frac{u^{n}-u^{n-1}}{\Delta t}\right) = \frac{u^{n+1}-2u^n + u^{n-1}}{\Delta t^2}
        {\thinspace .}
        

The discretization of initial conditions can in the operator notation
be expressed as

.. math::
   :label: _auto1
        
        [u = I]^0,\quad [D_{2t} u = 0]^0,
        
        

where the operator :math:`[D_{2t} u]^n` is defined as

.. math::
   :label: _auto2
        
        [D_{2t} u]^n = \frac{u^{n+1} - u^{n-1}}{2\Delta t}
        {\thinspace .}
        
        

.. _vib:impl1:

Implementation          (1)
===========================

.. _vib:impl1:solver:

Making a solver function
------------------------

The algorithm from the previous section is readily translated to
a complete Python function for computing and returning
:math:`u^0,u^1,\ldots,u^{N_t}` and :math:`t_0,t_1,\ldots,t_{N_t}`, given the
input :math:`I`, :math:`\omega`, :math:`\Delta t`, and :math:`T`:

.. code-block:: python

        import numpy as np
        import matplotlib.pyplot as plt
        
        def solver(I, w, dt, T):
            """
            Solve u'' + w**2*u = 0 for t in (0,T], u(0)=I and u'(0)=0,
            by a central finite difference method with time step dt.
            """
            dt = float(dt)
            Nt = int(round(T/dt))
            u = np.zeros(Nt+1)
            t = np.linspace(0, Nt*dt, Nt+1)
        
            u[0] = I
            u[1] = u[0] - 0.5*dt**2*w**2*u[0]
            for n in range(1, Nt):
                u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n]
            return u, t

We do a simple ``from module import *`` to make the code as close as
possible to MATLAB, although good programming habits would prefix
the ``numpy`` and ``matplotlib`` calls by (abbreviations of) the module
name.

[**hpl 2**: Refer to right section in decay book for prefix discussion.]

A function for plotting the numerical and the exact solution is also
convenient to have:

.. code-block:: python

        def u_exact(t, I, w):
            return I*np.cos(w*t)
        
        def visualize(u, t, I, w):
            plt.plot(t, u, 'r--o')
            t_fine = np.linspace(0, t[-1], 1001)  # very fine mesh for u_e
            u_e = u_exact(t_fine, I, w)
            plt.hold('on')
            plt.plot(t_fine, u_e, 'b-')
            plt.legend(['numerical', 'exact'], loc='upper left')
            plt.xlabel('t')
            plt.ylabel('u')
            dt = t[1] - t[0]
            plt.title('dt=%g' % dt)
            umin = 1.2*u.min();  umax = -umin
            plt.axis([t[0], t[-1], umin, umax])
            plt.savefig('tmp1.png');  plt.savefig('tmp1.pdf')

A corresponding main program calling these functions for a simulation
of a given number of periods (``num_periods``) may take the form

.. code-block:: python

        I = 1
        w = 2*pi
        dt = 0.05
        num_periods = 5
        P = 2*pi/w    #  one period
        T = P*num_periods
        u, t = solver(I, w, dt, T)
        visualize(u, t, I, w, dt)

Adjusting some of the input parameters via the command line can be
handy. Here is a code segment using the ``ArgumentParser`` tool in
the ``argparse`` module to define option value (``--option value``)
pairs on the command line:

.. code-block:: python

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('--I', type=float, default=1.0)
        parser.add_argument('--w', type=float, default=2*pi)
        parser.add_argument('--dt', type=float, default=0.05)
        parser.add_argument('--num_periods', type=int, default=5)
        a = parser.parse_args()
        I, w, dt, num_periods = a.I, a.w, a.dt, a.num_periods

Such parsing of the command line is explained in more detailed in
 the
"section on user interfaces": "..." in [Ref1]_.

[**hpl 3**: Fix reference to web document.]

A typical execution goes like

.. code-block:: text

        Terminal> python vib_undamped.py --num_periods 20 --dt 0.1

Computing :math:`u^{\prime}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In mechanical vibration applications one is often interested in
computing the velocity :math:`v(t)=u^{\prime}(t)` after :math:`u(t)` has been computed.
This can be done by a central difference,

.. math::
   :label: _auto3
        
        v(t_n)=u^{\prime}(t_n) \approx v^n = \frac{u^{n+1}-u^{n-1}}{2\Delta t} = [D_{2t}u]^n
        {\thinspace .}
        
        

This formula applies for all inner mesh points, :math:`n=1,\ldots,N_t-1`.
For :math:`n=0`, :math:`v(0)` is given by the initial condition on :math:`u^{\prime}(0)`,
and for :math:`n=N_t` we can use a one-sided, backward difference:

.. math::
         v^n=[D_t^-u]^n = \frac{u^{n} - u^{n-1}}{\Delta t}{\thinspace .}

Typical (scalar) code is

.. code-block:: python

        v = np.zeros_like(u)  # or v = np.zeros(len(u))
        # Use central difference for internal points
        for i in range(1, len(u)-1):
            v[i] = (u[i+1] - u[i-1])/(2*dt)
        # Use initial condition for u'(0) when i=0
        v[0] = 0
        # Use backward difference at the final mesh point
        v[-1] = (u[-1] - u[-2])/dt

We can get rid of the loop, which is slow for large :math:`N_t`, by
vectorizing the central difference. The above code segment
goes as follows in its vectorized version:

.. code-block:: python

        v = np.zeros_like(u)
        v[1:-1] = (u[2:] - u[:-2])/(2*dt)  # central difference
        v[0] = 0                           # boundary condition u'(0)
        v[-1] = (u[-1] - u[-2])/dt         # backward difference

.. _vib:ode1:verify:

Verification          (1)
-------------------------

Manual calculation
~~~~~~~~~~~~~~~~~~

The simplest type of verification, which is also instructive for understanding
the algorithm, is to compute :math:`u^1`, :math:`u^2`, and :math:`u^3`
with the aid of a calculator
and make a function for comparing these results with those from the ``solver``
function. The ``test_three_steps`` function in
the file `vib_undamped.py <http://tinyurl.com/nm5587k/vib/vib_undamped.py>`__
shows the details how we use the hand calculations to test the code:

.. code-block:: python

        def test_three_steps():
            from math import pi
            I = 1;  w = 2*pi;  dt = 0.1;  T = 1
            u_by_hand = np.array([1.000000000000000,
                                  0.802607911978213,
                                  0.288358920740053])
            u, t = solver(I, w, dt, T)
            diff = np.abs(u_by_hand - u[:3]).max()
            tol = 1E-14
            assert diff < tol

Testing very simple solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Constructing test problems where the exact solution is constant or
linear helps initial debugging and verification as one expects any
reasonable numerical method to reproduce such solutions to machine
precision.  Second-order accurate methods will often also reproduce a
quadratic solution.  Here :math:`[D_tD_tt^2]^n=2`, which is the exact
result. A solution :math:`u=t^2` leads to :math:`u^{\prime\prime}+\omega^2 u=2 + (\omega
t)^2\neq 0`.  We must therefore add a source in the equation: :math:`u^{\prime\prime} +
\omega^2 u = f` to allow a solution :math:`u=t^2` for :math:`f=(\omega t)^2`.  By
simple insertion we can show that the mesh function :math:`u^n = t_n^2` is
also a solution of the discrete equations.  :ref:`vib:exer:undamped:verify:linquad` asks you to carry out all
details to show that linear and quadratic solutions are solutions
of the discrete equations. Such results are very useful for debugging
and verification. You are strongly encouraged to do this problem now!

Checking convergence rates
~~~~~~~~~~~~~~~~~~~~~~~~~~

Empirical computation of convergence rates
yields a good method for verification. The method and its computational
are explained in detail for a simple ODE model in the section on `computing convergence rates <http://hplgit.github.io/decay-book/doc/pub/book/sphinx/._book006.html#computing-convergence-rates>`__
in [Ref1]_. Readers not familiar with the concept should
look up this reference before proceeding.

In the present problem, computing convergence rates means that we must

 * perform :math:`m` simulations with halved time steps: :math:`\Delta t_i=2^{-i}\Delta t_0`, :math:`i=0,\ldots,m-1`,

 * compute the :math:`L^2` norm of the error,
   :math:`E_i=\sqrt{\Delta t_i\sum_{n=0}^{N_t-1}(u^n-{u_{\small\mbox{e}}}(t_n))^2}` in each case,

 * estimate the convergence rates :math:`r_i` based on two consecutive
   experiments :math:`(\Delta t_{i-1}, E_{i-1})` and :math:`(\Delta t_{i}, E_{i})`,
   assuming :math:`E_i=C(\Delta t_i)^{r}` and :math:`E_{i-1}=C(\Delta t_{i-1})^{r}`.
   From these equations it follows that
   :math:`r = \ln (E_{i-1}/E_i)/\ln (\Delta t_{i-1}/\Delta t_i)`. Since this :math:`r`
   will vary with :math:`i`, we equip it with an index and call it :math:`r_{i-1}`,
   where :math:`i` runs from :math:`1` to :math:`m-1`.

The computed rates :math:`r_0,r_1,\ldots,r_{m-2}` hopefully converges to a
number, which hopefully is 2, the right one, in the present
problem. The convergence of the rates demands that the time steps
:math:`\Delta t_i` are sufficiently small for the error model :math:`E_i=(\Delta t_i)^r`
to be valid.

All the implementational details of computing the sequence
:math:`r_0,r_1,\ldots,r_{m-2}` appear below.

.. code-block:: python

        def convergence_rates(m, solver_function, num_periods=8):
            """
            Return m-1 empirical estimates of the convergence rate
            based on m simulations, where the time step is halved
            for each simulation.
            solver_function(I, w, dt, T) solves each problem, where T
            is based on simulation for num_periods periods.
            """
            from math import pi
            w = 0.35; I = 0.3       # just chosen values
            P = 2*pi/w              # period
            dt = P/30               # 30 time step per period 2*pi/w
            T = P*num_periods
        
            dt_values = []
            E_values = []
            for i in range(m):
                u, t = solver_function(I, w, dt, T)
                u_e = u_exact(t, I, w)
                E = np.sqrt(dt*np.sum((u_e-u)**2))
                dt_values.append(dt)
                E_values.append(E)
                dt = dt/2
        
            r = [np.log(E_values[i-1]/E_values[i])/
                 np.log(dt_values[i-1]/dt_values[i])
                 for i in range(1, m, 1)]
            return r

The expected convergence rate is 2, because we have used
a second-order finite
difference approximations :math:`[D_tD_tu]^n` to the ODE and a
second-order finite difference formula for the initial condition for
:math:`u^{\prime}`. Other theoretical error measures also points to
:math:`r=2`.

In the present problem, when :math:`\Delta t_0` corresponds to 30 time steps
per period, the returned ``r`` list has all its values equal to 2.00
(if rounded to two decimals). This amazing result means that all
:math:`\Delta t_i` values are well into the asymptotic regime where the
error model :math:`E_i = C(\Delta t_i)^r` is valid.

We can now construct a test function that computes convergence rates
and checks that the final (and usually the best) estimate is sufficiently
close to 2. Here, a rough tolerance of 0.1 is enough. This unit test
goes like

.. code-block:: python

        def test_convergence_rates():
            r = convergence_rates(m=5, solver_function=solver, num_periods=8)
            # Accept rate to 1 decimal place
            tol = 0.1
            assert abs(r[-1] - 2.0) < tol

The complete code appears in the file ``vib_undamped.py``.

Scaled model
------------

[**hpl 4**: Need reference to scaling book and maybe also decay book.]

It is advantageous to use dimensionless variables in simulations, because
fewer parameters need to be set. The present problem is made dimensionless
by introducing dimensionless variables :math:`\bar t = t/t_c` and :math:`\bar u = u/u_c`,
where :math:`t_c` and :math:`u_c` are characteristic scales for :math:`t` and :math:`u`,
respectively. The scaled ODE problem reads

.. math::
         \frac{u_c}{t_c^2}\frac{d^2\bar u}{d\bar t^2} + u_c\bar u = 0,\quad
        u_c\bar u(0) = I,\ \frac{u_c}{t_c}\frac{d\bar u}{d\bar t}(0)=0{\thinspace .}

A common choice is to take :math:`t_c` as one period of
the oscillations, :math:`t_c = 2\pi/w`, and :math:`u_c=I`.
This gives the dimensionless model

.. math::
   :label: vib:ode1:model:scaled
        
        \frac{d^2\bar u}{\bar t^2} + 4\pi^2 \bar u = 0,\quad \bar u(0)=1,\ 
        \bar u^{\prime}(0)=0{\thinspace .}
        
        

Observe that there are no physical parameters in :eq:`vib:ode1:model:scaled`!
We can therefore perform
a single numerical simulation :math:`\bar u(\bar t)` and afterwards
recover any :math:`u(t; \omega, I)` by

.. math::
         u(t;\omega, I) = u_c\bar u(t/t_c) = I\bar u(omega t/(2\pi)){\thinspace .}

We can easily check this assertion: the solution of the scaled problem
is :math:`\bar u(\bar t) = \cos(2\pi\bar t)`. The formula for :math:`u` in terms
of :math:`\bar u` gives :math:`u = I\cos(\omega t)`, which is nothing but the solution
of the original problem with dimensions.

The scaled model can by run by calling ``solver(I=1, w=2*pi, dt, T)``.
Each period is now 1 and ``T`` simply counts the number of periods.
Choosing ``dt`` as ``1./M`` gives ``M`` time steps per period.

.. _vib:ode1:longseries:

Long time simulations
=====================

Figure :ref:`vib:ode1:2dt` shows a comparison of the exact and numerical
solution for the scaled model :eq:`vib:ode1:model:scaled` with
:math:`\Delta t=0.1, 0.05`.
From the plot we make the following observations:

 * The numerical solution seems to have correct amplitude.

 * There is a angular frequency error which is reduced by reducing the time step.

 * The total angular frequency error grows with time.

By angular frequency error we mean that the numerical angular frequency differs
from the exact :math:`\omega`. This is evident by looking
at the peaks of the numerical solution: these have incorrect
positions compared with the peaks of the exact cosine solution. The
effect can be mathematical expressed by writing the numerical solution
as :math:`I\cos\tilde\omega t`, where :math:`\tilde\omega` is not exactly
equal to :math:`\omega`. Later, we shall mathematically
quantify this numerical angular frequency :math:`\tilde\omega`.

.. _vib:ode1:2dt:

.. figure:: fig-vib/vib_freq_err1.png
   :width: 800

   *Effect of halving the time step*

Using a moving plot window
--------------------------

In vibration problems it is often of interest to investigate the system's
behavior over long time intervals. Errors in the angular frequency accumulate
and become more visible as time grows. We can investigate long
time series by introducing a moving plot window that can move along with
the :math:`p` most recently computed periods of the solution. The
`SciTools <https://github.com/hplgit/scitools>`__ package contains
a convenient tool for this: ``MovingPlotWindow``. Typing
``pydoc scitools.MovingPlotWindow`` shows a demo and a description of its use.
The function below utilizes the moving plot window and is in fact
called by the ``main`` function the ``vib_undamped`` module
if the number of periods in the simulation exceeds 10.

.. code-block:: python

        def visualize_front(u, t, I, w, savefig=False, skip_frames=1):
            """
            Visualize u and the exact solution vs t, using a
            moving plot window and continuous drawing of the
            curves as they evolve in time.
            Makes it easy to plot very long time series.
            Plots are saved to files if savefig is True.
            Only each skip_frames-th plot is saved (e.g., if
            skip_frame=10, only each 10th plot is saved to file;
            this is convenient if plot files corresponding to
            different time steps are to be compared).
            """
            import scitools.std as st
            from scitools.MovingPlotWindow import MovingPlotWindow
            from math import pi
        
            # Remove all old plot files tmp_*.png
            import glob, os
            for filename in glob.glob('tmp_*.png'):
                os.remove(filename)
        
            P = 2*pi/w  # one period
            umin = 1.2*u.min();  umax = -umin
            dt = t[1] - t[0]
            plot_manager = MovingPlotWindow(
                window_width=8*P,
                dt=dt,
                yaxis=[umin, umax],
                mode='continuous drawing')
            frame_counter = 0
            for n in range(1,len(u)):
                if plot_manager.plot(n):
                    s = plot_manager.first_index_in_plot
                    st.plot(t[s:n+1], u[s:n+1], 'r-1',
                            t[s:n+1], I*cos(w*t)[s:n+1], 'b-1',
                            title='t=%6.3f' % t[n],
                            axis=plot_manager.axis(),
                            show=not savefig) # drop window if savefig
                    if savefig and n % skip_frames == 0:
                        filename = 'tmp_%04d.png' % frame_counter
                        st.savefig(filename)
                        print 'making plot file', filename, 'at t=%g' % t[n]
                        frame_counter += 1
                plot_manager.update(n)

We run the scaled problem (the default values for the command-line arguments
``--I`` and ``--w`` correspond to the scaled problem) for 40 periods with 20
time steps per period:

.. code-block:: text

        Terminal> python vib_undamped.py --dt 0.05 --num_periods 40

The moving plot window is invoked, and we can follow the numerical and exact
solutions as time progresses. From this demo we see that
the angular frequency error is small in the beginning, but it becomes more
prominent with time. A new run with :math:`\Delta t=0.1` (i.e., only 10 time steps per period)
clearly shows that the phase errors become significant even earlier
in the time series, deteriorating the solution further.

.. _vib:ode1:anim:

Making animations
-----------------

.. index:: making movies

.. index:: animation

.. index:: WebM (video format)

.. index:: Ogg (video format)

.. index:: MP4 (video format)

.. index:: Flash (video format)

.. index:: video formats

Producing standard video formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``visualize_front`` function stores all the plots in
files whose names are numbered:
``tmp_0000.png``, ``tmp_0001.png``, ``tmp_0002.png``,
and so on. From these files we may make a movie. The Flash
format is popular,

.. code-block:: text

        Terminal> ffmpeg -r 12 -i tmp_%04d.png -c:v flv movie.flv

The ``ffmpeg`` program can be replaced by the ``avconv`` program in
the above command if desired (but at the time of this writing it seems
to be more momentum in the ``ffmpeg`` project).
The ``-r`` option should come first and
describes the number of frames per second in the movie. The
``-i`` option describes the name of the plot files.
Other formats can be generated by changing the video codec
and equipping the video file with the right extension:

======  ============================  
Format       Codec and filename       
======  ============================  
Flash   ``-c:v flv movie.flv``        
MP4     ``-c:v libx264 movie.mp4``    
WebM    ``-c:v libvpx movie.webm``    
Ogg     ``-c:v libtheora movie.ogg``  
======  ============================  

.. index:: HTML5 video tag

The video file can be played by some video player like ``vlc``, ``mplayer``,
``gxine``, or ``totem``, e.g.,

.. code-block:: text

        Terminal> vlc movie.webm

A web page can also be used to play the movie. Today's standard is
to use the HTML5 ``video`` tag:

.. code-block:: html

        <video autoplay loop controls
               width='640' height='365' preload='none'>
        <source src='movie.webm'  type='video/webm; codecs="vp8, vorbis"'>
        </video>

Modern browsers do not support all of the video formats.
MP4 is needed to successfully play the videos on Apple devices
that use the Safari browser.
WebM is the preferred format for Chrome, Opera, Firefox, and Internet
Explorer v9+. Flash was a popular format, but older browsers that
required Flash can play MP4. All browsers that work with Ogg can also
work with WebM. This means that to have a video work in all browsers,
the video should be available in the MP4 and WebM formats.
The proper HTML code reads

.. code-block:: html

        <video autoplay loop controls
               width='640' height='365' preload='none'>
        <source src='movie.mp4'   type='video/mp4;
         codecs="avc1.42E01E, mp4a.40.2"'>
        <source src='movie.webm'  type='video/webm;
         codecs="vp8, vorbis"'>
        </video>

The MP4 format should appear first to ensure that Apple devices will
load the video correctly.


.. admonition:: Caution: number the plot files correctly

   To ensure that the individual plot frames are shown in correct order,
   it is important to number the files with zero-padded numbers
   (0000, 0001, 0002, etc.). The printf format ``%04d`` specifies an
   integer in a field of width 4, padded with zeros from the left.
   A simple Unix wildcard file specification like ``tmp_*.png``
   will then list the frames in the right order. If the numbers in the
   filenames were not zero-padded, the frame ``tmp_11.png`` would appear
   before ``tmp_2.png`` in the movie.




Paying PNG files in a web browser
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: scitools movie command

The ``scitools movie`` command can create a movie player for a set
of PNG files such that a web browser can be used to watch the movie.
This interface has the advantage that the speed of the movie can
easily be controlled, a feature that scientists often appreciate.
The command for creating an HTML with a player for a set of
PNG files ``tmp_*.png`` goes like

.. code-block:: text

        Terminal> scitools movie output_file=vib.html fps=4 tmp_*.png

The ``fps`` argument controls the speed of the movie ("frames per second").

To watch the movie, load the video file ``vib.html`` into some browser, e.g.,

.. code-block:: text

        Terminal> google-chrome vib.html  # invoke web page

Clicking on ``Start movie`` to see the result. Moving this movie to
some other place requires moving ``vib.html`` *and all the PNG files*
``tmp_*.png``:

.. code-block:: text

        Terminal> mkdir vib_dt0.1
        Terminal> mv tmp_*.png vib_dt0.1
        Terminal> mv vib.html vib_dt0.1/index.html

Making animated GIF files
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``convert`` program from the ImageMagick software suite can be
used to produce animated GIF files from a set of PNG files:

.. code-block:: text

        Terminal> convert -delay 25 tmp_vib*.png tmp_vib.gif

The ``-delay`` option needs an argument of the delay between each frame,
measured in 1/100 s, so 4 frames/s here gives 25/100 s delay.
Note, however, that in this particular example
with :math:`\Delta t=0.05` and 40 periods,
making an animated GIF file out of
the large number of PNG files is a very heavy process and not
considered feasible. Animated GIFs are best suited for animations with
not so many frames and where you want to see each frame and play them
slowly.

[**hpl 5**: Combine two simulations side by side!]

Using Bokeh to compare graphs
-----------------------------

Instead of a moving plot frame, one can use tools that allows panning
by the mouse. For example, we can show four periods of a signal in
a plot and then scroll with the mouse through the rest of the
simulation. The `Bokeh <http://bokeh.pydata.org/en/latest/docs/quickstart.html>`__ plotting library offers such tools, but the plot must be displayed in
a web browser. The documentation of Bokeh is excellent, so here we just
show how the library can be used to compare a set of :math:`u` curves corresponding
to long time simulations.

Imagine we have performed experiments for a set of :math:`\Delta t` values.
We want each curve, together with the exact solution, to appear in
a plot, and then arrange all plots in a grid-like fashion:

.. figure:: fig-vib/bokeh_gridplot1.png
   :width: 800

Furthermore, we want the axis to couple such that if we move into
the future in one plot, all the other plots follows (note the
displaced :math:`t` axes!):

.. figure:: fig-vib/bokeh_gridplot2.png
   :width: 800

A function for creating a Bokeh plot, given a list of ``u`` arrays
and corresponding ``t`` arrays, from different simulations, described
compactly in a list of strings ``legends``, takes the following form:

.. code-block:: python

        def bokeh_plot(u, t, legends, I, w, t_range, filename):
            """
            Make plots for u vs t using the Bokeh library.
            u and t are lists (several experiments can be compared).
            legens contain legend strings for the various u,t pairs.
            """
            if not isinstance(u, (list,tuple)):
                u = [u]  # wrap in list
            if not isinstance(t, (list,tuple)):
                t = [t]  # wrap in list
            if not isinstance(legends, (list,tuple)):
                legends = [legends]  # wrap in list
        
            import bokeh.plotting as plt
            plt.output_file(filename, mode='cdn', title='Comparison')
            # Assume that all t arrays have the same range
            t_fine = np.linspace(0, t[0][-1], 1001)  # fine mesh for u_e
            tools = 'pan,wheel_zoom,box_zoom,reset,'\ 
                    'save,box_select,lasso_select'
            u_range = [-1.2*I, 1.2*I]
            font_size = '8pt'
            p = []  # list of plot objects
            # Make the first figure
            p_ = plt.figure(
                width=300, plot_height=250, title=legends[0],
                x_axis_label='t', y_axis_label='u',
                x_range=t_range, y_range=u_range, tools=tools,
                title_text_font_size=font_size)
            p_.xaxis.axis_label_text_font_size=font_size
            p_.yaxis.axis_label_text_font_size=font_size
            p_.line(t[0], u[0], line_color='blue')
            # Add exact solution
            u_e = u_exact(t_fine, I, w)
            p_.line(t_fine, u_e, line_color='red', line_dash='4 4')
            p.append(p_)
            # Make the rest of the figures and attach their axes to
            # the first figure's axes
            for i in range(1, len(t)):
                p_ = plt.figure(
                    width=300, plot_height=250, title=legends[i],
                    x_axis_label='t', y_axis_label='u',
                    x_range=p[0].x_range, y_range=p[0].y_range, tools=tools,
                    title_text_font_size=font_size)
                p_.xaxis.axis_label_text_font_size = font_size
                p_.yaxis.axis_label_text_font_size = font_size
                p_.line(t[i], u[i], line_color='blue')
                p_.line(t_fine, u_e, line_color='red', line_dash='4 4')
                p.append(p_)
        
            # Arrange all plots in a grid with 3 plots per row
            grid = [[]]
            for i, p_ in enumerate(p):
                grid[-1].append(p_)
                if (i+1) % 3 == 0:
                    # New row
                    grid.append([])
            plot = plt.gridplot(grid, toolbar_location='left')
            plt.save(plot)
            plt.show(plot)

A particular example using the ``bokeh_plot`` function appears below.

.. code-block:: python

        def demo_bokeh():
            """Solve a scaled ODE u'' + u = 0."""
            from math import pi
            w = 1.0        # Scaled problem (frequency)
            P = 2*np.pi/w  # Period
            num_steps_per_period = [5, 10, 20, 40, 80]
            T = 40*P       # Simulation time: 40 periods
            u = []         # List of numerical solutions
            t = []         # List of corresponding meshes
            legends = []
            for n in num_steps_per_period:
                dt = P/n
                u_, t_ = solver(I=1, w=w, dt=dt, T=T)
                u.append(u_)
                t.append(t_)
                legends.append('# time steps per period: %d' % n)
            bokeh_plot(u, t, legends, I=1, w=w, t_range=[0, 4*P],
                       filename='tmp.html')

Using a line-by-line ascii plotter
----------------------------------

Plotting functions vertically, line by line, in the terminal window
using ascii characters only is a simple, fast, and convenient
visualization technique for long time series. Note that the time
axis then is positive downwards on the screen.
The tool
``scitools.avplotter.Plotter`` makes it easy to create such plots:

.. code-block:: python

        def visualize_front_ascii(u, t, I, w, fps=10):
            """
            Plot u and the exact solution vs t line by line in a
            terminal window (only using ascii characters).
            Makes it easy to plot very long time series.
            """
            from scitools.avplotter import Plotter
            import time
            from math import pi
            P = 2*pi/w
            umin = 1.2*u.min();  umax = -umin
        
            p = Plotter(ymin=umin, ymax=umax, width=60, symbols='+o')
            for n in range(len(u)):
                print p.plot(t[n], u[n], I*cos(w*t[n])), \ 
                      '%.1f' % (t[n]/P)
                time.sleep(1/float(fps))

The call ``p.plot`` returns a line of text, with the :math:`t` axis marked and
a symbol ``+`` for the first function (``u``) and ``o`` for the second
function (the exact solution). Here we append to this text
a time counter reflecting how many periods the current time point
corresponds to. A typical output (:math:`\omega =2\pi`, :math:`\Delta t=0.05`)
looks like this:

.. code-block:: text

                                      |                       o+      14.0
                                      |                      + o      14.0
                                      |                  +    o       14.1
                                      |             +     o           14.1
                                      |     +        o                14.2
                                     +|       o                       14.2
                             +        |                               14.2
                      +       o       |                               14.3
                 +     o              |                               14.4
              +   o                   |                               14.4
             +o                       |                               14.5
             o +                      |                               14.5
              o    +                  |                               14.6
                  o      +            |                               14.6
                       o        +     |                               14.7
                              o       | +                             14.7
                                      |        +                      14.8
                                      |       o       +               14.8
                                      |              o     +          14.9
                                      |                   o   +       14.9
                                      |                       o+      15.0

.. _vib:ode1:empirical:

Empirical analysis of the solution
----------------------------------

For oscillating functions like those in Figure :ref:`vib:ode1:2dt` we may
compute the amplitude and frequency (or period) empirically.
That is, we run through the discrete solution points :math:`(t_n, u_n)` and
find all maxima and minima points. The distance between two consecutive
maxima (or minima) points can be used as estimate of the local period,
while half the difference between the :math:`u` value at a maximum and a nearby
minimum gives an estimate of the local amplitude.

The local maxima are the points where

.. math::
   :label: _auto4
        
        u^{n-1} < u^n > u^{n+1},\quad n=1,\ldots,N_t-1,
        
        

and the local minima are recognized by

.. math::
   :label: _auto5
        
        u^{n-1} > u^n < u^{n+1},\quad n=1,\ldots,N_t-1
        {\thinspace .}
        
        

In computer code this becomes

.. code-block:: python

        def minmax(t, u):
            minima = []; maxima = []
            for n in range(1, len(u)-1, 1):
                if u[n-1] > u[n] < u[n+1]:
                    minima.append((t[n], u[n]))
                if u[n-1] < u[n] > u[n+1]:
                    maxima.append((t[n], u[n]))
            return minima, maxima

Note that the two returned objects are lists of tuples.

Let :math:`(t_i, e_i)`, :math:`i=0,\ldots,M-1`, be the sequence of all
the :math:`M` maxima points, where :math:`t_i`
is the time value and :math:`e_i` the corresponding :math:`u` value.
The local period can be defined as :math:`p_i=t_{i+1}-t_i`.
With Python syntax this reads

.. code-block:: python

        def periods(maxima):
            p = [extrema[n][0] - maxima[n-1][0]
                 for n in range(1, len(maxima))]
            return np.array(p)

The list ``p`` created by a list comprehension is converted to an array
since we probably want to compute with it, e.g., find the corresponding
frequencies ``2*pi/p``.

Having the minima and the maxima, the local amplitude can be
calculated as the difference between two neighboring minimum and
maximum points:

.. code-block:: python

        def amplitudes(minima, maxima):
            a = [(abs(maxima[n][1] - minima[n][1]))/2.0
                 for n in range(min(len(minima),len(maxima)))]
            return np.array(a)

The code segments are found in the file `vib_empirical_analysis.py <http://tinyurl.com/nm5587k/vib/vib_empirical_analysis.py>`__.

Since ``a[i]`` and ``p[i]`` correspond to
the :math:`i`-th amplitude estimate and the :math:`i`-th period estimate, respectively,
it is most convenient to visualize the ``a`` and ``p`` values with the
index ``i`` on the horizontal axis.
(There is no unique time point associated with either of these estimate
since values at two different time points were used in the
computations.)

In the analysis of very long time series, it is advantageous to
compute and plot ``p`` and ``a`` instead of :math:`u` to get an impression of
the development of the oscillations. Let us do this for the scaled
problem and :math:`\Delta t=0.1, 0.05, 0.01`.
A ready-made function

.. code-block:: python

        plot_empirical_freq_and_amplitude(u, t, I, w)

computes the empirical amplitudes and periods, and creates a plot
where the amplitudes and angular frequencies
are visualized together with the exact amplitude ``I``
and the exact angular frequency ``w``. We can make a little program
for creating the plot:

.. code-block:: python

        from vib_undamped import solver, plot_empirical_freq_and_amplitude
        from math import pi
        dt_values = [0.1, 0.05, 0.01]
        u_cases = []
        t_cases = []
        for dt in dt_values:
            # Simulate scaled problem for 40 periods
            u, t = solver(I=1, w=2*pi, dt=dt, T=40)
            u_cases.append(u)
            t_cases.append(t)
        plot_empirical_freq_and_amplitude(u_cases, t_cases, I=1, w=2*pi)

Figure :ref:`vib:ode1:fig:freq_ampl` shows the result: we clearly see that
lowering :math:`\Delta t` improves the angular frequency significantly, while the
amplitude seems to be more accurate.
The lines with
:math:`\Delta t=0.01`, corresponding to 100 steps per period, can hardly be
distinguished from the exact values. The next section shows how we
can get mathematical insight into why amplitudes are good and frequencies
are more inaccurate.

.. _vib:ode1:fig:freq_ampl:

.. figure:: fig-vib/empirical_ampl_freq.png
   :width: 800

   *Empirical amplitude and angular frequency for three cases of time steps*

.. Use it for very long time integration of CN! And of RK4!

.. _vib:ode1:analysis:

Analysis of the numerical scheme
================================

Deriving a solution of the numerical scheme
-------------------------------------------

After having seen the phase error grow with time in the previous
section, we shall now quantify this error through mathematical
analysis.  The key tool in the analysis will be to establish an exact
solution of the discrete equations.  The difference equation
:eq:`vib:ode1:step4` has constant coefficients and is
homogeneous. Such equations are known to have solutions on the form
:math:`u^n=CA^n`, where :math:`A` is some number
to be determined from the difference equation and :math:`C` is found as the
initial condition (:math:`C=I`).  Recall that :math:`n` in :math:`u^n` is a
superscript labeling the time level, while :math:`n` in :math:`A^n` is an
exponent.

With oscillating functions as solutions, the algebra will
be considerably simplified if we seek an :math:`A` on the form

.. math::
         A=e^{i\tilde\omega \Delta t},

and solve for the numerical frequency :math:`\tilde\omega` rather than
:math:`A`. Note that :math:`i=\sqrt{-1}` is the imaginary unit. (Using a
complex exponential function gives simpler arithmetics than working
with a sine or cosine function.)
We have

.. math::
        
        A^n = e^{i\tilde\omega \Delta t\, n}=e^{i\tilde\omega t} =
        \cos (\tilde\omega t) + i\sin(\tilde \omega t)
        {\thinspace .}
        

The physically relevant numerical solution can
be taken as the real part of this complex expression.

The calculations go as

.. math::
        
        [D_tD_t u]^n &= \frac{u^{n+1} - 2u^n + u^{n-1}}{\Delta t^2}\\ 
        &= I\frac{A^{n+1} - 2A^n + A^{n-1}}{\Delta t^2}\\ 
        &= \frac{I}{\Delta t^{2}}(e^{i\tilde\omega(t+\Delta t)} - 2e^{i\tilde\omega t} + e^{i\tilde\omega(t-\Delta t)})\\ 
        &= Ie^{i\tilde\omega t}\frac{1}{\Delta t^2}\left(e^{i\tilde\omega\Delta t} + e^{i\tilde\omega(-\Delta t)} - 2\right)\\ 
        &= Ie^{i\tilde\omega t}\frac{2}{\Delta t^2}\left(\cosh(i\tilde\omega\Delta t) -1 \right)\\ 
        &= Ie^{i\tilde\omega t}\frac{2}{\Delta t^2}\left(\cos(\tilde\omega\Delta t) -1 \right)\\ 
        &= -Ie^{i\tilde\omega t}\frac{4}{\Delta t^2}\sin^2(\frac{\tilde\omega\Delta t}{2})
        

The last line follows from the relation
:math:`\cos x - 1 = -2\sin^2(x/2)` (try ``cos(x)-1`` in
`wolframalpha.com <http://www.wolframalpha.com>`__ to see the formula).

The scheme :eq:`vib:ode1:step4`
with :math:`u^n=Ie^{i\omega\tilde\Delta t\, n}` inserted now gives

.. math::
   :label: _auto6
        
        -Ie^{i\tilde\omega t}\frac{4}{\Delta t^2}\sin^2(\frac{\tilde\omega\Delta t}{2})
        + \omega^2 Ie^{i\tilde\omega t} = 0,
        
        

which after dividing by :math:`Ie^{i\tilde\omega t}` results in

.. math::
   :label: _auto7
        
        \frac{4}{\Delta t^2}\sin^2(\frac{\tilde\omega\Delta t}{2}) = \omega^2
        {\thinspace .}
        
        

The first step in solving for the unknown :math:`\tilde\omega` is

.. math::
         \sin^2(\frac{\tilde\omega\Delta t}{2}) = \left(\frac{\omega\Delta t}{2}\right)^2
        {\thinspace .}
        

Then, taking the square root, applying the inverse sine function, and
multiplying by :math:`2/\Delta t`, results in

.. math::
   :label: vib:ode1:tildeomega
        
        \tilde\omega = \pm \frac{2}{\Delta t}\sin^{-1}\left(\frac{\omega\Delta t}{2}\right)
        {\thinspace .}
        
        

The first observation of :eq:`vib:ode1:tildeomega` tells that
there is a phase error since the numerical frequency :math:`\tilde\omega`
never equals the exact frequency :math:`\omega`. But how good is
the approximation :eq:`vib:ode1:tildeomega`? That is, what
is the error :math:`\omega - \tilde\omega` or :math:`\tilde\omega/\omega`?
Taylor series expansion
for small :math:`\Delta t` may give an expression that is easier to understand
than the complicated function in :eq:`vib:ode1:tildeomega`:

.. code-block:: ipy

        >>> from sympy import *
        >>> dt, w = symbols('dt w')
        >>> w_tilde_e = 2/dt*asin(w*dt/2)
        >>> w_tilde_series = w_tilde_e.series(dt, 0, 4)
        >>> print w_tilde_series
        w + dt**2*w**3/24 + O(dt**4)

This means that

.. See vib_symbolic.py for computations with sympy

.. math::
   :label: vib:ode1:tildeomega:series
        
        \tilde\omega = \omega\left( 1 + \frac{1}{24}\omega^2\Delta t^2\right)
        + {\mathcal{O}(\Delta t^4)}
        {\thinspace .}
        
        

The error in the numerical frequency is of second-order in
:math:`\Delta t`, and the error vanishes as :math:`\Delta t\rightarrow 0`.
We see that :math:`\tilde\omega > \omega` since the term :math:`\omega^3\Delta t^2/24 >0`
and this is by far the biggest term in the series expansion for small
:math:`\omega\Delta t`. A numerical frequency that is too large gives an oscillating
curve that oscillates too fast and therefore "lags behind" the exact
oscillations, a feature that can be seen in the left plot in Figure
:ref:`vib:ode1:2dt`.

Figure :ref:`vib:ode1:tildeomega:plot` plots the discrete frequency
:eq:`vib:ode1:tildeomega`
and its approximation :eq:`vib:ode1:tildeomega:series` for :math:`\omega =1` (based
on the program `vib_plot_freq.py <http://tinyurl.com/nm5587k/vib/vib_plot_freq.py>`__).
Although :math:`\tilde\omega` is a function of :math:`\Delta t` in
:eq:`vib:ode1:tildeomega:series`,
it is misleading to think of :math:`\Delta t` as the important
discretization parameter. It is the product :math:`\omega\Delta t` that is
the key discretization parameter. This quantity reflects the
*number of time steps per period* of the oscillations.
To see this, we set :math:`P=N_P\Delta t`, where :math:`P` is the length of
a period, and :math:`N_P` is the number of time steps during a period.
Since :math:`P` and :math:`\omega` are related by :math:`P=2\pi/\omega`,
we get that :math:`\omega\Delta t = 2\pi/N_P`, which shows that
:math:`\omega\Delta t` is directly related to :math:`N_P`.

The plot shows
that at least :math:`N_P\sim 25-30` points per period are necessary for reasonable
accuracy, but this depends on the length of the simulation (:math:`T`) as
the total phase error due to the frequency error grows linearly with time
(see :ref:`vib:exer:phase:err:growth`).

.. _vib:ode1:tildeomega:plot:

.. figure:: fig-vib/discrete_freq.png
   :width: 400

   *Exact discrete frequency and its second-order series expansion*

.. _vib:ode1:analysis:sol:

Exact discrete solution
-----------------------

Perhaps more important than the :math:`\tilde\omega = \omega + {\cal O}(\Delta t^2)`
result found above is the fact that we have an exact discrete solution of
the problem:

.. math::
   :label: vib:ode1:un:exact
        
        u^n = I\cos\left(\tilde\omega n\Delta t\right),\quad
        \tilde\omega = \frac{2}{\Delta t}\sin^{-1}\left(\frac{\omega\Delta t}{2}\right)
        {\thinspace .}
        
        

We can then compute the error mesh function

.. math::
   :label: vib:ode1:en
        
        e^n = {u_{\small\mbox{e}}}(t_n) - u^n =
        I\cos\left(\omega n\Delta t\right) - I\cos\left(\tilde\omega n\Delta t\right){\thinspace .}
        
        

From the formula :math:`\cos 2x - \cos 2y = -2\sin(x-y)\sin(x+y)` we can
rewrite :math:`e^n` so the expression is easier to interpret:

.. math::
   :label: vib:ode1:en2
        
        e^n = -2I\sin\left(t\frac{1}{2}\left( \omega - \tilde\omega\right)\right)
        \sin\left(t\frac{1}{2}\left( \omega + \tilde\omega\right)\right){\thinspace .}
        
        

The error mesh function is ideal for verification purposes
and you are strongly encouraged to make a test based on :eq:`vib:ode1:un:exact`
by doing :ref:`vib:exer:discrete:omega`.

.. _vib:ode1:analysis:conv:

Convergence
-----------

We can use :eq:`vib:ode1:tildeomega:series`, :eq:`vib:ode1:en`, or
:eq:`vib:ode1:en2` to show *convergence* of the
numerical scheme, i.e., :math:`e^n\rightarrow 0` as :math:`\Delta t\rightarrow 0`.
We have that

.. math::
        
        \lim_{\Delta t\rightarrow 0}
        \tilde\omega = \lim_{\Delta t\rightarrow 0}
        \frac{2}{\Delta t}\sin^{-1}\left(\frac{\omega\Delta t}{2}\right)
        = \omega,
        

by L'Hopital's rule or simply asking ``sympy`` or
`WolframAlpha <http://www.wolframalpha.com/input/?i=%282%2Fx%29*asin%28w*x%2F2%29+as+x-%3E0>`__ about the limit:

.. code-block:: python

        >>> import sympy as sym
        >>> dt, w = sym.symbols('x w')
        >>> sym.limit((2/dt)*sym.asin(w*dt/2), dt, 0, dir='+')
        w

Also :eq:`vib:ode1:tildeomega:series` can be used to establish
this result that
:math:`\tilde\omega\rightarrow\omega`. It then follows from the expression(s)
for :math:`e^n` that :math:`e^n\rightarrow 0`.

The global error
----------------

.. index::
   single: error; global

To achieve more analytical insight into the nature of the global error,
we can Taylor expand the error mesh function :eq:`vib:ode1:en`.
Since :math:`\tilde\omega` in :eq:`vib:ode1:tildeomega`
contains :math:`\Delta t` in the denominator we use the series expansion
for :math:`\tilde\omega` inside the cosine function. A relevant ``sympy``
session is

.. code-block:: python

        >>> from sympy import *
        >>> dt, w, t = symbols('dt w t')
        >>> w_tilde_e = 2/dt*asin(w*dt/2)
        >>> w_tilde_series = w_tilde_e.series(dt, 0, 4)
        >>> w_tilde_series
        w + dt**2*w**3/24 + O(dt**4)

Series expansions in ``sympy`` have the inconvenient ``O()`` term that
prevents further calculations with the series. We can use the
``removeO()`` command to get rid of the ``O()`` term:

.. code-block:: python

        >>> w_tilde_series = w_tilde_series.removeO()
        >>> w_tilde_series
        dt**2*w**3/24 + w

Using this ``w_tilde_series`` expression
for :math:`\tilde w` in :eq:`vib:ode1:en`,
dropping :math:`I` (which is a common factor), and performing a series
expansion of the error yields

.. code-block:: python

        >>> error = cos(w*t) - cos(w_tilde_series*t)
        >>> error.series(dt, 0, 6)
        dt**2*t*w**3*sin(t*w)/24 + dt**4*t**2*w**6*cos(t*w)/1152 + O(dt**6)

Since we are mainly interested in the leading-order term in
such expansions (the term with lowest power in :math:`\Delta t` and
goes most slowly to zero), we use the ``.as_leading_term(dt)``
construction to pick out this term:

.. code-block:: python

        >>> error.series(dt, 0, 6).as_leading_term(dt)
        dt**2*t*w**3*sin(t*w)/24

The last result
means that the leading order global (true) error at a point :math:`t`
is proportional to :math:`\omega^3t\Delta t^2`. Now, :math:`t` is related
to :math:`\Delta t` through :math:`t=n\Delta t`. The factor
:math:`\sin(\omega t)` can at most be 1, so we use this value to
bound the leading-order expression to its maximum value

.. math::
         e^n = \frac{1}{24}n\omega^3\Delta t^3{\thinspace .}

This is the dominating term of the error *at a point*.

We are interested in the accumulated global error, which can
be taken as the :math:`\ell^2` norm of :math:`e^n`.
The norm is simply computed by summing contributions from all mesh
points:

.. math::
         ||e^n||_{\ell^2}^2 = \Delta t\sum_{n=0}^{N_t} \frac{1}{24^2}n^2\omega^6\Delta t^6
        =\frac{1}{24^2}\omega^6\Delta t^7 \sum_{n=0}^{N_t} n^2{\thinspace .}

The sum :math:`\sum_{n=0}^{N_t} n^2` is approximately equal to
:math:`\frac{1}{3}N_t^3`. Replacing :math:`N_t` by :math:`T/\Delta t` and taking
the square root gives the expression

.. math::
         ||e^n||_{\ell^2} = \frac{1}{24}\sqrt{\frac{T^3}{3}}\omega^3\Delta t^2{\thinspace .}

This is our expression for the global (or integrated) error.
The main result from this expression is that also the global error
is proportional to :math:`\Delta t^2`.

Stability
---------

Looking at :eq:`vib:ode1:un:exact`, it appears that the numerical
solution has constant and correct amplitude, but an error in the
angular frequency. A constant amplitude is not necessarily the case,
however! To see this, note that if only :math:`\Delta t` is large
enough, the magnitude of the argument to :math:`\sin^{-1}` in
:eq:`vib:ode1:tildeomega` may be larger than 1, i.e.,
:math:`\omega\Delta t/2 > 1`. In this case, :math:`\sin^{-1}(\omega\Delta t/2)`
has a complex value and therefore :math:`\tilde\omega` becomes complex.
Type, for example, ``asin(x)`` in
`wolframalpha.com <http://www.wolframalpha.com>`__ to see basic properties of :math:`\sin^{-1} (x)`).

A complex :math:`\tilde\omega` can be written :math:`\tilde\omega = \tilde\omega_r +
i\tilde\omega_i`. Since :math:`\sin^{-1}(x)` has a *negative* imaginary part for
:math:`x>1`, :math:`\tilde\omega_i < 0`, which means that
:math:`e^{i\tilde\omega t}=e^{-\tilde\omega_i t}e^{i\tilde\omega_r t}`
will lead to exponential growth in time because
:math:`e^{-\tilde\omega_i t}` with :math:`\tilde\omega_i <0` has a positive
exponent.

.. index:: stability criterion


.. admonition:: Stability criterion

   We do not tolerate growth in the amplitude since such growth is not
   present in the exact solution. Therefore, we
   must impose a *stability criterion*  that
   the argument in the inverse sine function leads
   to real and not complex values of :math:`\tilde\omega`. The stability
   criterion reads
   
   .. math::
      :label: _auto8
           
           \frac{\omega\Delta t}{2} \leq 1\quad\Rightarrow\quad
           \Delta t \leq \frac{2}{\omega}
           {\thinspace .}




With :math:`\omega =2\pi`, :math:`\Delta t > \pi^{-1} = 0.3183098861837907` will give
growing solutions. Figure :ref:`vib:ode1:dt:unstable`
displays what happens when :math:`\Delta t =0.3184`,
which is slightly above the critical value: :math:`\Delta t =\pi^{-1} + 9.01\cdot
10^{-5}`.

.. _vib:ode1:dt:unstable:

.. figure:: fig-vib/vib_unstable.png
   :width: 400

   *Growing, unstable solution because of a time step slightly beyond the stability limit*

About the accuracy at the stability limit
-----------------------------------------

An interesting question is whether the stability condition
:math:`\Delta t < 2/\omega` is unfortunate, or more precisely:
would it be meaningful to take larger time steps to speed up computations?
The answer is a clear no. At the stability limit, we have that
:math:`\sin^{-1}\omega\Delta t/2 = \sin^{-1} 1 = \pi/2`, and therefore
:math:`\tilde\omega = \pi/\Delta t`. (Note that the approximate formula
:eq:`vib:ode1:tildeomega:series` is very inaccurate for this
value of :math:`\Delta t` as it predicts :math:`\tilde\omega = 2.34/pi`, which is
a 25 percent reduction.) The corresponding
period of the numerical solution
is :math:`\tilde P=2\pi/\tilde\omega = 2\Delta t`, which means that there is
just one time step :math:`\Delta t` between a peak (maximum)
and a `through <https://simple.wikipedia.org/wiki/Wave_(physics)>`__
(minimum) in the
numerical solution. This is the shortest possible wave that can be
represented in the mesh! In other words, it is not meaningful to
use a larger time step than the stability limit.

Also, the error in angular frequency
when :math:`\Delta t = 2/\omega` is severe: Figure
:ref:`vib:ode1:dt:stablimit` shows a comparison of the numerical and
analytical solution with :math:`\omega = 2\pi` and
:math:`\Delta t = 2/\omega = \pi^{-1}`. Already after one period, the
numerical solution has a through while the exact solution has a peak (!).
The error in frequency when :math:`\Delta t` is at the stability limit
becomes :math:`\omega - \tilde\omega = \omega(1-\pi/2)\approx -0.57\omega`.
The corresponding error in the period is :math:`P - \tilde P \approx 0.36P`.
The error after :math:`m` periods is then :math:`0.36mP`. This error has reached
half a period when :math:`m=1/(2\cdot 0.36)\approx 1.38`, which theoretically
confirms the observations in Figure :ref:`vib:ode1:dt:stablimit`
that the numerical solution is a through ahead of a peak already after
one and a half period. Consequently, :math:`\Delta t` should be chosen much
less than the stability limit to achieve meaningful numerical computations.

.. _vib:ode1:dt:stablimit:

.. figure:: fig-vib/vib_stability_limit.png
   :width: 400

   *Numerical solution with :math:`\Delta t` exactly at the stability limit*


.. admonition:: Summary

   From the accuracy and stability
   analysis we can draw three important conclusions:
   
   1. The key parameter in the formulas is :math:`p=\omega\Delta t`.
      The period of oscillations is :math:`P=2\pi/\omega`, and the
      number of time steps per period is :math:`N_P=P/\Delta t`.
      Therefore, :math:`p=\omega\Delta t = 2\pi N_P`, showing that the
      critical parameter is the number of time steps per period.
      The smallest possible :math:`N_P` is 2, showing that :math:`p\in (0,\pi]`.
   
   2. Provided :math:`p\leq 2`, the amplitude of the numerical solution is
      constant.
   
   3. The ratio of the numerical angular frequency and the exact
      one is
      :math:`\tilde\omega/\omega \approx 1 + \frac{1}{24}p^2`.
      The error :math:`\frac{1}{24}p^2` leads to wrongly displaced peaks of the numerical
      solution, and the error in peak location grows linearly with time
      (see :ref:`vib:exer:phase:err:growth`).




.. _vib:model2x2:

Alternative schemes based on 1st-order equations
================================================

A standard technique for solving second-order ODEs is
to rewrite them as a system of first-order ODEs and then choose a
solution strategy from the
vast collection of methods for first-order ODE systems.
Given the second-order ODE problem

.. math::
         u^{\prime\prime} + \omega^2 u = 0,\quad u(0)=I,\ u^{\prime}(0)=0,

we introduce the auxiliary variable :math:`v=u^{\prime}` and express the ODE problem
in terms of first-order derivatives of :math:`u` and :math:`v`:

.. math::
   :label: vib:model2x2:ueq
        
        u^{\prime} = v,
        
        

.. math::
   :label: vib:model2x2:veq
          
        v' = -\omega^2 u
        
        {\thinspace .}
        

The initial conditions become :math:`u(0)=I` and :math:`v(0)=0`.

The Forward Euler scheme
------------------------

A Forward Euler approximation to our :math:`2\times 2` system of ODEs
:eq:`vib:model2x2:ueq`-:eq:`vib:model2x2:veq`
becomes

.. math::
   :label: _auto9
        
        \lbrack D_t^+ u = v\rbrack^n,
        \lbrack D_t^+ v = -\omega^2 u\rbrack^n,
        
        

or written out,

.. math::
   :label: vib:undamped:FE1
        
        u^{n+1} = u^n + \Delta t v^n,
        
        

.. math::
   :label: vib:undamped:FE2
          
        v^{n+1} = v^n -\Delta t \omega^2 u^n
        
        {\thinspace .}
        

Let us briefly compare this Forward Euler method with the
centered difference scheme for the second-order differential
equation. We have from :eq:`vib:undamped:FE1` and
:eq:`vib:undamped:FE2` applied at levels :math:`n` and :math:`n-1` that

.. math::
         u^{n+1} = u^n + \Delta t v^n = u^n + \Delta t (v^{n-1} -\Delta t \omega^2 u^{n-1}{\thinspace .}

Since from :eq:`vib:undamped:FE1`

.. math::
         v^{n-1} = \frac{1}{\Delta t}(u^{n}-u^{n-1}),

it follows that

.. math::
         u^{n+1} = 2u^n - u^{n-1} -\Delta t^2\omega^2 u^{n-1},

which is very close to the centered difference scheme, but
the last term is evaluated at :math:`t_{n-1}` instead of :math:`t_n`.
Dividing by :math:`\Delta t^2`, the left-hand side is an approximation to
:math:`u^{\prime\prime}` at :math:`t_n`, while the right-hand side is sampled at :math:`t_{n-1}`.
All terms should be sampled at the same mesh point, so using
:math:`\omega^2 u^{n-1}` instead of :math:`\omega^2 u^n` is an inconsistency
in the scheme. This inconsistency turns out to be rather
crucial for the accuracy of
the Forward Euler method applied to vibration problems.

The Backward Euler scheme
-------------------------

A Backward Euler approximation the ODE system is equally easy to
write up in the operator notation:

.. math::
   :label: _auto10
        
        \lbrack D_t^- u = v\rbrack^{n+1},
        
        

.. math::
   :label: _auto11
          
        \lbrack D_t^- v = -\omega u\rbrack^{n+1} {\thinspace .}
        
        

This becomes a coupled system for :math:`u^{n+1}` and :math:`v^{n+1}`:

.. math::
   :label: vib:undamped:BE1
        
        u^{n+1} - \Delta t v^{n+1} = u^{n},
        
        

.. math::
   :label: vib:undamped:BE2
          
        v^{n+1} + \Delta t \omega^2 u^{n+1} = v^{n}
        
        {\thinspace .}
        

We can compare :eq:`vib:undamped:BE1`-:eq:`vib:undamped:BE2` with
the centered scheme :eq:`vib:ode1:step4`
for the second-order differential equation.
To this end, we eliminate :math:`v^{n+1}` in :eq:`vib:undamped:BE1`
using :eq:`vib:undamped:BE2` solved with respect to :math:`v^{n+1}`.
Thereafter, we eliminate :math:`v^n` using :eq:`vib:undamped:BE1`
solved with respect to :math:`v^{n+1}` and replacing :math:`n+1` by :math:`n`.
The resulting equation involving only :math:`u^{n+1}`, :math:`u^n`, and :math:`u^{n-1}`
can be ordered as

.. math::
         \frac{u^{n+1}-2u^n+u^{n-1}}{\Delta t^2} = -\omega^2 u^{n+1},

which has almost the same form as the centered scheme for the
second-order differential equation, but the right-hand side is
evaluated at :math:`u^{n+1}` and not :math:`u^n`. This inconsistent sampling
of terms has a dramatic effect on the numerical solution.

The Crank-Nicolson scheme
-------------------------

The Crank-Nicolson scheme takes this form in the operator notation:

.. math::
   :label: _auto12
        
        \lbrack D_t u = \overline{v}^t\rbrack^{n+\frac{1}{2}},
        
        

.. math::
   :label: _auto13
          
        \lbrack D_t v = -\omega \overline{u}^t\rbrack^{n+\frac{1}{2}}
        {\thinspace .}
        
        

Writing the equations out shows that this is also a coupled system:

.. math::
   :label: _auto14
        
        u^{n+1} - \frac{1}{2}\Delta t v^{n+1} = u^{n} + \frac{1}{2}\Delta t v^{n},
        
        

.. math::
   :label: _auto15
          
        v^{n+1} + \frac{1}{2}\Delta t \omega^2 u^{n+1} = v^{n}
        - \frac{1}{2}\Delta t \omega^2 u^{n}
        {\thinspace .}
        
        

To see the nature of this approximation, and that it is actually
very promising, we write the equations as follows

.. math::
   :label: vib:undamped:CN3a
        
        u^{n+1} - u^n = \frac{1}{2}\Delta t(v^{n+1} + v^n),
        
        

.. math::
   :label: vib:undamped:CN4a
          
        v^{n+1}  = v^n -\frac{1}{2}\Delta t(u^{n+1} + u^n),
        
        

and add the latter at the previous time level as well:

.. math::
   :label: vib:undamped:CN4b1
        
        v^{n}  = v^{n-1} -\frac{1}{2}\Delta t(u^{n} + u^{n-1})
        
        

We can also rewrite :eq:`vib:undamped:CN3a` at the previous time level
as

.. math::
   :label: vib:undamped:CN4b
        
        v^{n+1} + v^n = \frac{2}{\Delta t}(u^{n+1} - u^n){\thinspace .}
        
        

Inserting :eq:`vib:undamped:CN4a` for :math:`v^{n+1}` in
:eq:`vib:undamped:CN3a` and
:eq:`vib:undamped:CN4b1` for :math:`v^{n}` in
:eq:`vib:undamped:CN3a` yields after some reordering:

.. math::
         u^{n+1} - n^n = \frac{1}{2}(-\frac{1}{2}\Delta t\omega^2
        (u^{n+1} + 2u^n + u^{n-1}) + v^ + v^{n-1}){\thinspace .}

Now, :math:`v^n + v^{n-1}` can be eliminated by means of
:eq:`vib:undamped:CN4b`. The result becomes

.. math::
   :label: vib:undamped:CN5
        
        u^{n+1} - 2u^n + u^{n-1} = \Delta t^2\omega^2
        \frac{1}{4}(u^{n+1} + 2u^n + u^{n-1}){\thinspace .}
        
        

We have that

.. math::
         \frac{1}{4}(u^{n+1} + 2u^n + u^{n-1}) \approx u^n + {\mathcal{O}(\Delta t^2)},

meaning that :eq:`vib:undamped:CN5` is an approximation to
the centered scheme :eq:`vib:ode1:step4` for the second-order ODE where
the sampling error in the term :math:`\Delta t^2\omega^2 u^n` is of the same
order as the approximation errors in the finite differences, i.e.,
:math:`{\mathcal{O}(\Delta t^2)}`. The Crank-Nicolson scheme written as
:eq:`vib:undamped:CN5` therefore has consistent sampling of all
terms at the same time point :math:`t_n`. The implication is a much better
method than the Forward and Backward Euler schemes.

.. _vib:model2x2:compare:

Comparison of schemes
---------------------

We can easily compare methods like the ones above (and many more!)
with the aid of the
`Odespy <https://github.com/hplgit/odespy>`__ package. Below is
a sketch of the code.

.. code-block:: python

        import odespy
        import numpy as np
        
        def f(u, t, w=1):
            u, v = u  # u is array of length 2 holding our [u, v]
            return [v, -w**2*u]
        
        def run_solvers_and_plot(solvers, timesteps_per_period=20,
                                 num_periods=1, I=1, w=2*np.pi):
            P = 2*np.pi/w  # duration of one period
            dt = P/timesteps_per_period
            Nt = num_periods*timesteps_per_period
            T = Nt*dt
            t_mesh = np.linspace(0, T, Nt+1)
        
            legends = []
            for solver in solvers:
                solver.set(f_kwargs={'w': w})
                solver.set_initial_condition([I, 0])
                u, t = solver.solve(t_mesh)

There is quite some more code dealing with plots also, and we refer
to the source file `vib_undamped_odespy.py <http://tinyurl.com/nm5587k/vib/vib_undamped_odespy.py>`__
for details. Observe that keyword arguments in ``f(u,t,w=1)`` can
be supplied through a solver parameter ``f_kwargs`` (dictionary of
additional keyword arguments to ``f``).

Specification of the Forward Euler, Backward Euler, and
Crank-Nicolson schemes is done like this:

.. code-block:: python

        solvers = [
            odespy.ForwardEuler(f),
            # Implicit methods must use Newton solver to converge
            odespy.BackwardEuler(f, nonlinear_solver='Newton'),
            odespy.CrankNicolson(f, nonlinear_solver='Newton'),
            ]

.. index:: phase plane plot

The ``vib_undamped_odespy.py``
program makes two plots of the computed solutions with the various
methods in the ``solvers`` list: one plot with :math:`u(t)` versus :math:`t`, and
one *phase plane plot* where :math:`v` is plotted against :math:`u`.
That is, the phase plane plot is the curve :math:`(u(t),v(t))` parameterized
by :math:`t`. Analytically, :math:`u=I\cos(\omega t)` and :math:`v=u^{\prime}=-\omega I\sin(\omega t)`.
The exact curve :math:`(u(t),v(t))` is therefore an ellipse, which often
looks like a circle in a plot if the axes are automatically scaled. The
important feature, however, is that exact curve :math:`(u(t),v(t))` is
closed and repeats itself for every period. Not all numerical schemes
are capable of doing that, meaning that the amplitude instead shrinks or
grows with time.

Figure
:ref:`vib:ode1:1st:odespy:theta:phaseplane` show the results. Note that
Odespy applies the label MidpointImplicit for what we have specified
as ``CrankNicolson`` in the code (``CrankNicolson`` is just a synonym for
class ``MidpointImplicit`` in the Odespy code).
The Forward Euler scheme in Figure
:ref:`vib:ode1:1st:odespy:theta:phaseplane` has a pronounced spiral
curve, pointing to the fact that the amplitude steadily grows, which
is also evident in Figure :ref:`vib:ode1:1st:odespy:theta`.
The Backward Euler scheme has a similar feature, except that the
spriral goes inward and the amplitude is significantly damped.  The
changing amplitude and the sprial form decreases with decreasing time
step.  The Crank-Nicolson scheme looks much more
accurate.  In fact, these plots tell that the Forward and Backward
Euler schemes are not suitable for solving our ODEs with oscillating
solutions.

.. _vib:ode1:1st:odespy:theta:phaseplane:

.. figure:: fig-vib/vib_theta_1_pp.png
   :width: 800

   *Comparison of classical schemes in the phase plane for two time step values*

.. _vib:ode1:1st:odespy:theta:

.. figure:: fig-vib/vib_theta_1_u.png
   :width: 800

   *Comparison of solution curves for classical schemes*

Runge-Kutta methods
-------------------

We may run two popular standard methods for first-order ODEs, the 2nd-
and 4th-order Runge-Kutta methods, to see how they perform. Figures
:ref:`vib:ode1:1st:odespy:RK:phaseplane` and
:ref:`vib:ode1:1st:odespy:RK` show the solutions with larger :math:`\Delta
t` values than what was used in the previous two plots.

.. _vib:ode1:1st:odespy:RK:phaseplane:

.. figure:: fig-vib/vib_RK_1_pp.png
   :width: 800

   *Comparison of Runge-Kutta schemes in the phase plane*

.. _vib:ode1:1st:odespy:RK:

.. figure:: fig-vib/vib_RK_1_u.png
   :width: 800

   *Comparison of Runge-Kutta schemes*

The visual impression is that the
4th-order Runge-Kutta method is very accurate, under all circumstances
in these tests, while the 2nd-order scheme suffers from amplitude errors
unless the time step is very small.

The corresponding results for the Crank-Nicolson scheme are shown in
Figure :ref:`vib:ode1:1st:odespy:CN:long:phaseplane`.
It is clear that the Crank-Nicolson
scheme outperforms the 2nd-order Runge-Kutta method. Both schemes have
the same order of accuracy :math:`{\mathcal{O}(\Delta t^2)}`, but their differences
in the accuracy that matters in a real physical application is very
clearly pronounced in this example.  :ref:`vib:exer:undamped:odespy` invites you to investigate how the amplitude
is computed by a series of famous methods for first-order ODEs.

.. _vib:ode1:1st:odespy:CN:long:phaseplane:

.. figure:: fig-vib/vib_CN_10_pp.png
   :width: 800

   *Long-time behavior of the Crank-Nicolson scheme in the phase plane*

Analysis of the Forward Euler scheme
------------------------------------

We may try to find exact solutions of the discrete
equations :eq:`vib:undamped:FE1`-:eq:`vib:undamped:FE2`
in the Forward Euler method. An "ansatz"
is

.. math::
        
        u^n &= IA^n,\\ 
        v^n &= qIA^n,
        

where :math:`q` and :math:`A` are unknown numbers. We could have used a complex
exponential form :math:`e^{i\tilde\omega n\Delta t}` since we get
oscillatory form, but the oscillations grow in the Forward Euler
method, so the numerical frequency :math:`\tilde\omega` will be complex
anyway (producing an exponentially growing amplitude). Therefore, it is
easier to just work with potentially complex :math:`A` and :math:`q` as introduced
above.

The Forward Euler scheme leads to

.. math::
        
        A &= 1 + \Delta t q,\\ 
        A &= 1 - \Delta t\omega^2 q^{-1}{\thinspace .}
        

We can easily eliminate :math:`A`, get :math:`q^2 + \omega^2=0`, and solve for

.. math::
         q = \pm i\omega,

which gives

.. math::
         A = 1 \pm \Delta t i\omega{\thinspace .}

We shall take the real part of :math:`A^n` as the solution. The two
values of :math:`A` are complex conjugates, and the real part of
:math:`A^n` will be the same for both roots. This is easy to realize if
we rewrite the complex numbers in polar form,
which is also convenient
for further analysis and understanding.
The polar form :math:`re^{i\theta}` of a complex number :math:`x+iy` has
:math:`r=\sqrt{x^2+y^2}` and :math:`\theta = \tan^{-1}(y/x)`.
Hence, the polar form of the two values for :math:`A` become

.. math::
         1 \pm \Delta t i\omega = \sqrt{1+\omega^2\Delta t^2}e^{\pm i\tan^{-1}(\omega\Delta t)}{\thinspace .}

Now it is very easy to compute :math:`A^n`:

.. math::
         (1 \pm \Delta t i\omega)^n = (1+\omega^2\Delta t^2)^{n/2}e^{\pm ni\tan^{-1}(\omega\Delta t)}{\thinspace .}

Since :math:`\cos (\theta n) = \cos (-\theta n)`, the real part of the two
numbers become the same. We therefore continue with the solution that has
the plus sign.

The general solution is :math:`u^n = CA^n`, where
:math:`C` is a constant determined from the initial condition:
:math:`u^0=C=I`. We have :math:`u^n=IA^n` and
:math:`v^n=qIA^n`. The final solutions
are just the real part of the expressions in polar form:

.. math::
   :label: _auto16
        
        u^n  =
        I(1+\omega^2\Delta t^2)^{n/2}\cos (n\tan^{-1}(\omega\Delta t)),
        
        

.. math::
   :label: _auto17
          
        v^n =- \omega
        I(1+\omega^2\Delta t^2)^{n/2}\sin (n\tan^{-1}(\omega\Delta t)){\thinspace .}
        
        

The expression :math:`(1+\omega^2\Delta t^2)^{n/2}` causes growth of
the amplitude, since a number greater than one is raised to a positive
exponent :math:`n/2`. We can develop a series expression to better understand
the formula for the amplitude. Introducing :math:`p=\omega\Delta t` as the
key variable and using ``sympy`` gives

.. code-block:: python

        >>> from sympy import *
        >>> p = symbols('p', real=True)
        >>> n = symbols('n', integer=True, positive=True)
        >>> amplitude = (1 + p**2)**(n/2)
        >>> amplitude.series(p, 0, 4)
        1 + n*p**2/2 + O(p**4)

The amplitude goes like :math:`1 + \frac{1}{2} n\omega^2\Delta t^2`, clearly growing
linearly in time (with :math:`n`).

We can also investigate the error in the angular frequency by a
series expansion:

.. code-block:: python

        >>> n*atan(p).series(p, 0, 4)
        n*(p - p**3/3 + O(p**4))

This means that the solution for :math:`u^n` can be written as

.. math::
         u^n = (1 + \frac{1}{2} n\omega^2\Delta t^2 + \Oof(\Delta t^4))
        \cos\left(\omega t - \frac{1}{3}\omega t\Delta t^2 + {\mathcal{O}(\Delta t^4)}\right)
        {\thinspace .}

The error in the angular frequency is of the same order as in the
scheme :eq:`vib:ode1:step4` for the second-order ODE, but error
in the amplitude is severe.

.. _vib:model1:energy:

Energy considerations
=====================

.. index:: mechanical energy

.. index:: energy principle

The observations of various methods in the previous section can be
better interpreted if we compute a quantity reflecting
the total *energy of the system*. It turns out that this quantity,

.. math::
         E(t) = \frac{1}{2}(u^{\prime})^2 + \frac{1}{2}\omega^2u^2,

is *constant* for all :math:`t`. Checking that :math:`E(t)` really remains constant
brings evidence that the numerical computations are sound.
It turns out that :math:`E` is proportional to the mechanical energy
in the system. Conservation of energy is
much used to check numerical simulations.

Derivation of the energy expression
-----------------------------------

We start out with multiplying

.. math::
         u^{\prime\prime} + \omega^2 u = 0,

by :math:`u^{\prime}` and integrating from :math:`0` to :math:`T`:

.. math::
         \int_0^T u^{\prime\prime}u^{\prime} dt + \int_0^T\omega^2 u u^{\prime} dt = 0{\thinspace .}

Observing that

.. math::
         u^{\prime\prime}u^{\prime} = \frac{d}{dt}\frac{1}{2}(u^{\prime})^2,\quad uu^{\prime} = \frac{d}{dt} {\frac{1}{2}}u^2,

we get

.. math::
        
        \int_0^T (\frac{d}{dt}\frac{1}{2}(u^{\prime})^2 + \frac{d}{dt} \frac{1}{2}\omega^2u^2)dt = E(T) - E(0)=0,
        

where we have introduced

.. math::
   :label: vib:model1:energy:balance1
        
        E(t) = \frac{1}{2}(u^{\prime})^2 + \frac{1}{2}\omega^2u^2{\thinspace .}
        
        

The important result from this derivation is that the total energy
is constant:

.. math::
         E(t) = E(0){\thinspace .}


.. admonition:: :math:`E(t)` is closely related to the system's energy

   The quantity :math:`E(t)` derived above is physically not the mechanical energy of a
   vibrating mechanical system, but the energy per unit mass. To see this,
   we start with Newton's second law :math:`F=ma` (:math:`F` is the sum of forces, :math:`m`
   is the mass of the system, and :math:`a` is the acceleration).
   The displacement :math:`u` is related to :math:`a` through
   :math:`a=u^{\prime\prime}`. With a spring force as the only force we have :math:`F=-ku`, where
   :math:`k` is a spring constant measuring the stiffness of the spring.
   Newton's second law then implies the differential equation
   
   .. math::
            -ku = mu^{\prime\prime}\quad\Rightarrow mu^{\prime\prime} + ku = 0{\thinspace .}
   
   This equation of motion can be turned into an energy balance equation
   by finding the work done by each term during a time interval :math:`[0,T]`.
   To this end, we multiply the equation by :math:`du=u^{\prime}dt` and integrate:
   
   .. math::
            \int_0^T muu^{\prime}dt + \int_0^T kuu^{\prime}dt = 0{\thinspace .}
   
   The result is
   
   .. math::
            \tilde E(t) = E_k(t) + E_p(t) = 0,
   
   where
   
   .. math::
      :label: vib:model1:energy:kinetic
           
           E_k(t) = \frac{1}{2}mv^2,\quad v=u^{\prime},
           
           
   
   is the *kinetic energy* of the system, and
   
   .. math::
      :label: vib:model1:energy:potential
           
           E_p(t) = {\frac{1}{2}}ku^2
           
           
   
   is the *potential energy*. The sum :math:`\tilde E(t)` is the total mechanical energy.
   The derivation demonstrates the famous energy principle that, under
   the right physical circumstances, any
   change in the kinetic energy is due to a change in potential energy
   and vice versa. (This principle breaks down when we introduce damping
   in system, as we do in the section :ref:`vib:model2`.)
   
   The equation :math:`mu^{\prime\prime}+ku=0` can be divided by :math:`m` and written as
   :math:`u^{\prime\prime} + \omega^2u=0` for :math:`\omega=\sqrt{k/m}`. The energy expression
   :math:`E(t)=\frac{1}{2}(u^{\prime})^2 + \frac{1}{2}\omega^2u^2` derived earlier is then
   :math:`\tilde E(t)/m`, i.e., mechanical energy per unit mass.




Energy of the exact solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analytically, we have :math:`u(t)=I\cos\omega t`, if :math:`u(0)=I` and :math:`u^{\prime}(0)=0`,
so we can easily check that the energy evolution and confirm that :math:`E(t)`
is constant:

.. math::
         E(t) = {\frac{1}{2}}I^2 (-\omega\sin\omega t)^2
        + \frac{1}{2}\omega^2 I^2 \cos^2\omega t
        = \frac{1}{2}\omega^2 (\sin^2\omega t + \cos^2\omega t) = \frac{1}{2}\omega^2
        {\thinspace .}
        

An error measure based on energy
--------------------------------

The constant energy is well expressed by its initial value :math:`E(0)`, so that
the error in mechanical energy can be computed as a mesh function by

.. math::
   :label: _auto18
        
        e_E^n = \frac{1}{2}\left(\frac{u^{n+1}-u^{n-1}}{2\Delta t}\right)^2
        + \frac{1}{2}\omega^2 (u^n)^2 - E(0),
        \quad n=1,\ldots,N_t-1,
        
        

where

.. math::
         E(0) = {\frac{1}{2}}V^2 + \frac{1}{2}\omega^2I^2,

if :math:`u(0)=I` and :math:`u^{\prime}(0)=V`. Note that we have used
a centered approximation to :math:`u^{\prime}`: :math:`\boldsymbol{u}^{\prime}(t_n)\approx
[D_{2t}u]^n`.

A useful norm of the mesh function :math:`e_E^n`
for the discrete mechanical energy
can be the maximum absolute value of :math:`e_E^n`:

.. math::
         ||e_E^n||_{\ell^\infty} = \max_{1\leq n <N_t} |e_E^n|{\thinspace .}

Alternatively, we can compute other norms involving integration over
all mesh points, but we are often interested in worst case deviation
of the energy, and then the maximum value is of particular relevance.

A vectorized Python implementation takes the form

.. code-block:: python

        # import numpy as np and compute u, t
        dt = t[1]-t[0]
        E = 0.5*((u[2:] - u[:-2])/(2*dt))**2 + 0.5*w**2*u[1:-1]**2
        E0 = 0.5*V**2 + 0.5**w**2*I**2
        e_E = E - E0
        e_E_norm = np.abs(e_E).max()

The convergence rates of the quantity ``e_E_norm`` can be used for verification.
The value of ``e_E_norm`` is also useful for comparing schemes
through their ability to preserve energy. Below is a table demonstrating
the error in total energy for various schemes. We clearly see that
the Crank-Nicolson and 4th-order Runge-Kutta schemes are superior to
the 2nd-order Runge-Kutta method and better compared to the Forward
and Backward Euler schemes.

=====================  ==========  ================  ========================================  
        Method         :math:`T`   :math:`\Delta t`  :math:`\max \left\vert e_E^n\right\vert`  
=====================  ==========  ================  ========================================  
Forward Euler          :math:`1`   :math:`0.05`      :math:`1.113\cdot 10^{2}`                 
Forward Euler          :math:`1`   :math:`0.025`     :math:`3.312\cdot 10^{1}`                 
Backward Euler         :math:`1`   :math:`0.05`      :math:`1.683\cdot 10^{1}`                 
Backward Euler         :math:`1`   :math:`0.025`     :math:`1.231\cdot 10^{1}`                 
Runge-Kutta 2nd-order  :math:`1`   :math:`0.1`       :math:`8.401`                             
Runge-Kutta 2nd-order  :math:`1`   :math:`0.05`      :math:`9.637\cdot 10^{-1}`                
Crank-Nicolson         :math:`1`   :math:`0.05`      :math:`9.389\cdot 10^{-1}`                
Crank-Nicolson         :math:`1`   :math:`0.025`     :math:`2.411\cdot 10^{-1}`                
Runge-Kutta 4th-order  :math:`1`   :math:`0.1`       :math:`2.387`                             
Runge-Kutta 4th-order  :math:`1`   :math:`0.05`      :math:`6.476\cdot 10^{-1}`                
Crank-Nicolson         :math:`10`  :math:`0.1`       :math:`3.389`                             
Crank-Nicolson         :math:`10`  :math:`0.05`      :math:`9.389\cdot 10^{-1}`                
Runge-Kutta 4th-order  :math:`10`  :math:`0.1`       :math:`3.686`                             
Runge-Kutta 4th-order  :math:`10`  :math:`0.05`      :math:`6.928\cdot 10^{-1}`                
=====================  ==========  ================  ========================================  

[**hpl 6**: The error reductions are not directly in accordance with the order of the schemes, probably caused by :math:`\Delta t` not being in the asympotic regime.]

.. Should build a verification test on the energy error.

.. Link phase plane plot to energy

.. A phase plane plot shows the curve :math:`(u(t), u^{\prime}(t))`.

.. _vib:model2x2:EulerCromer:

The Euler-Cromer method
=======================

While the 4th-order Runge-Kutta method and a
Crank-Nicolson scheme work well for vibration equation modeled as a
first-order ODE system,
both were inferior to the straightforward centered
difference scheme for the second-order equation
:math:`u^{\prime\prime}+\omega^2u=0`. However, there is a similarly successful scheme
available for the first-order system :math:`u^{\prime}=v`, :math:`v'=-\omega^2u`, to be
presented next.

.. index:: forward-backward Euler-Cromer scheme

Forward-backward discretization
-------------------------------

The idea is to apply a Forward Euler discretization to the first
equation and a Backward Euler discretization to the second. In operator
notation this is stated as

.. math::
   :label: _auto19
        
        \lbrack D_t^+u = v\rbrack^n,
        
        

.. math::
   :label: _auto20
          
        \lbrack D_t^-v = -\omega u\rbrack^{n+1}
        {\thinspace .}
        
        

We can write out the formulas and collect the unknowns on the left-hand side:

.. math::
   :label: vib:model2x2:EulerCromer:ueq1
        
        u^{n+1} = u^n + \Delta t v^n,
        
        

.. math::
   :label: vib:model2x2:EulerCromer:veq1
          
        v^{n+1} = v^n -\Delta t \omega^2u^{n+1}
        
        {\thinspace .}
        

We realize that after :math:`u^{n+1}` has been computed from
:eq:`vib:model2x2:EulerCromer:ueq1`, it may be used directly
in
:eq:`vib:model2x2:EulerCromer:veq1` to compute :math:`v^{n+1}`.

In physics, it is more common to update the :math:`v` equation first, with
a forward difference, and thereafter the :math:`u` equation, with a backward
difference that applies the most recently computed :math:`v` value:

.. math::
   :label: vib:model2x2:EulerCromer:veq1b
        
        v^{n+1} = v^n -\Delta t \omega^2u^{n},
        
        

.. math::
   :label: vib:model2x2:EulerCromer:ueq1b
          
        u^{n+1} = u^n + \Delta t v^{n+1}{\thinspace .}
        
        

The advantage of ordering the ODEs as in
:eq:`vib:model2x2:EulerCromer:veq1b`-:eq:`vib:model2x2:EulerCromer:ueq1b`
becomes evident
when consider complicated models. Such models are included if
we write our vibration ODE more generally as

.. math::
         \ddot u + g(u, u^{\prime}, t)=0{\thinspace .}

We can rewrite this second-order ODE as two first-order ODEs,

.. math::
        
        v' &= -g(u,v,t),\\ 
        u^{\prime} &= v{\thinspace .}
        

This rewrite allows the following scheme to be used:

.. math::
        
        v^{n+1} &= v^n -\Delta t\, g(u^n,v^n,t),\\ 
        u^{n+1} &= u^n + \Delta t\, v^{n+1}{\thinspace .}
        

We realize that the first update works well with any :math:`g` since old
values :math:`u^n` and :math:`v^n` are used. Switching the equations would
demand :math:`u^n{+1}` and :math:`v^{n+1}` values in :math:`g`.

.. Despite using a backward difference, there is no need to solve a coupled

.. system for :math:`u^{n+1}` and :math:`v^{n+1}` because the structure of the ODEs

.. allows :eq:`vib:model2x2:EulerCromer:ueq1`

The scheme
:eq:`vib:model2x2:EulerCromer:veq1b`-:eq:`vib:model2x2:EulerCromer:ueq1b`
goes under several names: forward-backward scheme, `semi-implicit Euler method <http://en.wikipedia.org/wiki/Semi-implicit_Euler_method>`__, semi-explicit Euler,
symplectic Euler,
Newton-Stormer-Verlet,
and Euler-Cromer.
We shall stick to the latter name.
Since both time discretizations are based on first-order difference
approximation, one may think that the scheme is only of first-order,
but this is not true: the use of a forward and then a backward
difference make errors cancel so that the overall error in the scheme
is :math:`{\mathcal{O}(\Delta t^2)}`. This is explained below.

.. _vib:model2x2:EulerCromer:equiv:

Equivalence with the scheme for the second-order ODE
----------------------------------------------------

We may eliminate the :math:`v^n` variable from
:eq:`vib:model2x2:EulerCromer:ueq1`-:eq:`vib:model2x2:EulerCromer:veq1`
or
:eq:`vib:model2x2:EulerCromer:veq1b`-:eq:`vib:model2x2:EulerCromer:ueq1b`.
The :math:`v^{n+1}` term in :eq:`vib:model2x2:EulerCromer:veq1b` can
be eliminated from :eq:`vib:model2x2:EulerCromer:ueq1b`:

.. math::
   :label: vib:model2x2:EulerCromer:elim1
        
        u^{n+1} = u^n + \Delta t (v^n - \omega^2\Delta t^2 u^n){\thinspace .}
        
        

The :math:`v^{n}` quantity can be expressed by :math:`u^n` and :math:`u^{n-1}`
using :eq:`vib:model2x2:EulerCromer:ueq1b`:

.. math::
         v^{n} = \frac{u^n - u^{n-1}}{\Delta t},
        

and when this is inserted in :eq:`vib:model2x2:EulerCromer:elim1` we get

.. math::
   :label: _auto21
        
        u^{n+1} = 2u^n - u^{n-1} - \Delta t^2 \omega^2u^{n},
        
        

which is nothing but the centered scheme :eq:`vib:ode1:step4`!
The two seemingly different numerical methods are mathematically
equivalent. Consequently,
the previous analysis of
:eq:`vib:ode1:step4` also applies to the Euler-Cromer
method. In particular, the amplitude is constant, given that the stability
criterion is fulfilled, but there is always an angular frequency error
:eq:`vib:ode1:tildeomega:series`. :ref:`vib:exer:EulerCromer:analysis`
gives guidance on how to derive the exact discrete solution of
the two equations in the Euler-Cromer method.

Although the Euler-Cromer scheme and the method :eq:`vib:ode1:step4` are
equivalent, there could be differences in the way they handle
the initial conditions. Let is look into this topic.
The initial condition :math:`u^{\prime}=0` means :math:`u^{\prime}=v=0`.  From
:eq:`vib:model2x2:EulerCromer:ueq1b` we get :math:`v^1=-\omega^2 u^0`
and :math:`u^1=u^0 - \omega^2\Delta t^2 u^0`. When using
a centered approximation of :math:`u^{\prime}(0)=0` combined with the
discretization :eq:`vib:ode1:step4` of the second-order ODE, we
get :math:`u^1=u^0 - \frac{1}{2}\omega^2\Delta t^2 u^0`. The difference
is :math:`\frac{1}{2}\omega^2\Delta t^2 u^0`, which is of second order in :math:`\Delta t`,
seemingly
consistent with the overall error in the scheme for the differential equation
model.

A different view can also be taken.
If we approximate :math:`u^{\prime}(0)=0` by a backward difference,
:math:`(u^0-u^{-1})/\Delta t =0`, we get :math:`u^{-1}=u^0`, and when combined
with :eq:`vib:ode1:step4`, it results in
:math:`u^1=u^0 - \omega^2\Delta t^2 u^0`. This means that
the Euler-Cromer method based on
:eq:`vib:model2x2:EulerCromer:ueq1b`-:eq:`vib:model2x2:EulerCromer:veq1b`
corresponds to using only a first-order approximation to the initial condition
in the method from the section :ref:`vib:ode1:fdm`.

Correspondingly, using the formulation
:eq:`vib:model2x2:EulerCromer:ueq1`-:eq:`vib:model2x2:EulerCromer:veq1`
with :math:`v^n=0` leads to :math:`u^1=u^0`, which can be interpreted as using
a forward difference approximation for the initial condition :math:`u^{\prime}(0)=0`.
Both Euler-Cromer formulations lead to slightly different values for
:math:`u^1` compared to the method in the section :ref:`vib:ode1:fdm`.
The error is :math:`\frac{1}{2}\omega^2\Delta t^2 u^0` and of the same order
as the overall scheme.

.. _vib:model2x2:EulerCromer:impl:

Implementation          (2)
---------------------------

The function below, found in `vib_EulerCromer.py <http://tinyurl.com/nm5587k/vib/vib_EulerCromer.py>`__ implements the Euler-Cromer scheme
:eq:`vib:model2x2:EulerCromer:veq1b`-:eq:`vib:model2x2:EulerCromer:ueq1b`:

.. code-block:: python

        import numpy as np
        
        def solver(I, w, dt, T):
            """
            Solve v' = - w**2*u, u'=v for t in (0,T], u(0)=I and v(0)=0,
            by an Euler-Cromer method.
            """
            dt = float(dt)
            Nt = int(round(T/dt))
            u = np.zeros(Nt+1)
            v = np.zeros(Nt+1)
            t = np.linspace(0, Nt*dt, Nt+1)
        
            v[0] = 0
            u[0] = I
            for n in range(0, Nt):
                v[n+1] = v[n] - dt*w**2*u[n]
                u[n+1] = u[n] + dt*v[n+1]
            return u, v, t

Since the Euler-Cromer scheme is equivalent to the finite difference
method for the second-order ODE :math:`u^{\prime\prime}+\omega^2u=0` (see the section :ref:`vib:model2x2:EulerCromer:equiv`), the performance of the above
``solver`` function is the same as for the ``solver`` function in the section :ref:`vib:impl1`. The only difference is the formula for the first time
step, as discussed above.  This deviation in the Euler-Cromer scheme
means that the discrete solution listed in the section :ref:`vib:ode1:analysis:sol` is not a solution of the Euler-Cromer
scheme!

To verify the implementation of the Euler-Cromer method we
can adjust ``v[1]`` so that the computer-generated values can be
compared with the formula
:eq:`vib:ode1:un:exact` from in the section :ref:`vib:ode1:analysis:sol`. This
adjustment is done in an alternative solver function, ``solver_ic_fix``
in ``vib_EulerCromer.py``. Since we now have an exact solution of the
discrete equations available, we can write a test function
``test_solver`` for checking the equality of computed values with the
formula :eq:`vib:ode1:un:exact`:

.. code-block:: python

        def test_solver():
            """
            Test solver with fixed initial condition against
            equivalent scheme for the 2nd-order ODE u'' + u = 0.
            """
            I = 1.2; w = 2.0; T = 5
            dt = 2/w  # longest possible time step
            u, v, t = solver_ic_fix(I, w, dt, T)
            from vib_undamped import solver as solver2  # 2nd-order ODE
            u2, t2 = solver2(I, w, dt, T)
            error = np.abs(u - u2).max()
            tol = 1E-14
            assert error < tol

Another function, ``demo``,
visualizes the difference between Euler-Cromer scheme and the scheme
:eq:`vib:ode1:step4`
for the second-oder ODE, arising from the mismatch in the first time level.

[**hpl 7**: Odespy's Euler-Cromer, but it needs more work with the example code.]

.. is anything gained? is v of higher order than D_2t u from the

.. other approach, i.e., if we need v, is this alg better? Probably not

.. since v is related u through a difference

.. make exercises:

.. investigate how important the u^1 wrong formula really is on

.. convergence rate

.. new file: genealizations, systems,

.. new file: apps

.. exercise: damping analysis, see geophysics book first...

The velocity Verlet algorithm
-----------------------------

Another very popular algorithm for vibration problems :math:`u^{\prime\prime}+\omega^2u=0`
can be derived as follows. First, we step :math:`u` forward from :math:`t_n` to
:math:`t_{n+1}` using a three-term Taylor series,

.. math::
         u(t_{n+1}) = u(t_n) + u^{\prime}(t_n)\Delta t + \frac{1}{2}u^{\prime\prime}(t_n)\Delta t^2{\thinspace .}

Using :math:`u^{\prime}=v` and :math:`u^{\prime\prime}=-\omega^2u`, we get the updating formula

.. math::
         u^{n+1} = u^n + v^n\Delta t - \frac{1}{2}\Delta^2\omega^2u^n{\thinspace .}

Second, the first-order equation for :math:`v`,

.. math::
         v'=-\omega^2u,

is discretized by a centered difference
in a Crank-Nicolson fashion at :math:`t_{n+\frac{1}{2}}`:

.. math::
         \frac{v^{n+1}-v^n}{\Delta t} = -\omega^2\frac{1}{2}(u^n + u^{n+1}){\thinspace .}

To summarize, we have the scheme

.. math::
   :label: vib:model2x2:Verlet:dueq
        
        u^{n+1} = u^n + v^n\Delta t - \frac{1}{2}\Delta^2\omega^2u^n
         
        

.. math::
   :label: vib:model2x2:Verlet:dveq
          
        v^{n+1} = v^n -\frac{1}{2}\Delta t\omega^2 (u^n + u^{n+1}),
        
        

known as the *velocity Verlet* algorithm.
Observe that this scheme is explicit since :math:`u^{n+1}` in
:eq:`vib:model2x2:Verlet:dveq` is already computed
from :eq:`vib:model2x2:Verlet:dueq`.

The algorithm can be straightforwardly implemented as shown below (the
code appears in the file `vib_undamped_velocity_Verlet.py <http://tinyurl.com/nm5587k/vib/vib_undamped_velocity_Verlet.py>`__).

.. code-block:: python

        from vib_undamped import convergence_rates, main
        
        def solver(I, w, dt, T, return_v=False):
            """
            Solve u'=v, v'=-w**2*u for t in (0,T], u(0)=I and v(0)=0,
            by the velocity Verlet method with time step dt.
            """
            dt = float(dt)
            Nt = int(round(T/dt))
            u = np.zeros(Nt+1)
            v = np.zeros(Nt+1)
            t = np.linspace(0, Nt*dt, Nt+1)
        
            u[0] = I
            v[0] = 0
            for n in range(Nt):
                u[n+1] = u[n] + v[n]*dt - 0.5*dt**2*w**2*u[n]
                v[n+1] = v[n] - 0.5*dt*w**2*(u[n] + u[n+1])
            if return_v:
                return u, v, t
            else:
                # Return just u and t as in the vib_undamped.py's solver
                return u, t

We provide the option that this ``solver`` function returns the same data
as the ``solver`` function from the section :ref:`vib:impl1:solver` (if ``return_v``
is ``False``), but we may return ``v`` along with ``u`` and ``t``.

The error in the Taylor series expansion behind
:eq:`vib:model2x2:Verlet:dueq` is :math:`{\mathcal{O}(\Delta t^3)}`, while the error
in the central difference for :math:`v` is :math:`{\mathcal{O}(\Delta t^2)}`.  The overall
error is then no better than :math:`{\mathcal{O}(\Delta t^2)}`, which can be verified
empirically using the ``convergence_rates`` function from
:ref:`vib:ode1:verify`:

.. code-block:: python

        >>> import vib_undamped_velocity_Verlet as m
        >>> m.convergence_rates(4, solver_function=m.solver)
        [2.0036366687367346, 2.0009497328124835, 2.000240105995295]

.. The output confirms that the overall convergence rate is 2.

.. !split

.. _vib:model2:

Generalization: damping, nonlinear spring, and external excitation
==================================================================

.. index:: nonlinear restoring force

.. index:: nonlinear spring

.. index:: forced vibrations

We shall now generalize the simple model problem from
the section :ref:`vib:model1` to include a possibly nonlinear damping term :math:`f(u^{\prime})`,
a possibly nonlinear spring (or restoring) force :math:`s(u)`, and
some external excitation :math:`F(t)`:

.. math::
   :label: vib:ode2
        
        mu^{\prime\prime} + f(u^{\prime}) + s(u) = F(t),\quad u(0)=I,\ u^{\prime}(0)=V,\ t\in (0,T]
        {\thinspace .}
        
        

We have also included a possibly nonzero initial value of :math:`u^{\prime}(0)`.
The parameters :math:`m`, :math:`f(u^{\prime})`, :math:`s(u)`, :math:`F(t)`, :math:`I`, :math:`V`, and :math:`T` are
input data.

There are two main types of damping (friction) forces: linear :math:`f(u^{\prime})=bu`, or
quadratic :math:`f(u^{\prime})=bu^{\prime}|u^{\prime}|`. Spring systems often feature linear
damping, while air resistance usually gives rise to quadratic damping.
Spring forces are often linear: :math:`s(u)=cu`, but nonlinear versions
are also common, the most famous is the gravity force on a pendulum
that acts as a spring with :math:`s(u)\sim \sin(u)`.

.. _vib:ode2:fdm:flin:

A centered scheme for linear damping
------------------------------------

Sampling :eq:`vib:ode2` at a mesh point :math:`t_n`, replacing
:math:`u^{\prime\prime}(t_n)` by :math:`[D_tD_tu]^n`, and :math:`u^{\prime}(t_n)` by :math:`[D_{2t}u]^n` results
in the discretization

.. math::
   :label: _auto22
        
        [mD_tD_t u + f(D_{2t}u) + s(u) = F]^n,
        
        

which written out means

.. math::
   :label: vib:ode2:step3b
        
        m\frac{u^{n+1}-2u^n + u^{n-1}}{\Delta t^2}
        + f(\frac{u^{n+1}-u^{n-1}}{2\Delta t}) + s(u^n) = F^n,
        
        

where :math:`F^n` as usual means :math:`F(t)` evaluated at :math:`t=t_n`.
Solving :eq:`vib:ode2:step3b` with respect to the unknown
:math:`u^{n+1}` gives a problem: the :math:`u^{n+1}` inside the :math:`f` function
makes the equation *nonlinear* unless :math:`f(u^{\prime})` is a linear function,
:math:`f(u^{\prime})=bu^{\prime}`. For now we shall assume that :math:`f` is linear in :math:`u^{\prime}`.
Then

.. math::
   :label: vib:ode2:step3b2
        
        m\frac{u^{n+1}-2u^n + u^{n-1}}{\Delta t^2}
        + b\frac{u^{n+1}-u^{n-1}}{2\Delta t} + s(u^n) = F^n,
        
        

which gives an explicit formula for :math:`u` at each
new time level:

.. math::
   :label: vib:ode2:step4
        
        u^{n+1} = (2mu^n + (\frac{b}{2}\Delta t - m)u^{n-1} +
        \Delta t^2(F^n - s(u^n)))(m + \frac{b}{2}\Delta t)^{-1}
        
        {\thinspace .}
        

For the first time step we need to discretize :math:`u^{\prime}(0)=V`
as :math:`[D_{2t}u = V]^0` and combine
with :eq:`vib:ode2:step4` for :math:`n=0`. The discretized initial condition
leads to

.. math::
   :label: vib:ode2:ic:du
        
        u^{-1} = u^{1} - 2\Delta t V,
        
        

which inserted in :eq:`vib:ode2:step4` for :math:`n=0` gives an equation
that can be solved for
:math:`u^1`:

.. math::
   :label: vib:ode2:step4b
        
        u^1 = u^0 + \Delta t\, V
        + \frac{\Delta t^2}{2m}(-bV - s(u^0) + F^0)
        {\thinspace .}
        
        

.. _vib:ode2:fdm:fquad:

A centered scheme for quadratic damping
---------------------------------------

When :math:`f(u^{\prime})=bu^{\prime}|u^{\prime}|`, we get a quadratic equation for :math:`u^{n+1}`
in :eq:`vib:ode2:step3b`. This equation can be straightforwardly
solved by the well-known formula for the roots of a quadratic equation.
However, we can also avoid the nonlinearity by introducing
an approximation with an error of order no higher than what we
already have from replacing derivatives with finite differences.

.. index:: geometric mean

.. index::
   single: averaging; geometric

We start with :eq:`vib:ode2` and only replace
:math:`u^{\prime\prime}` by :math:`D_tD_tu`, resulting in

.. math::
   :label: vib:ode2:quad:idea1
        
        [mD_tD_t u + bu^{\prime}|u^{\prime}| + s(u) = F]^n{\thinspace .}
        
        

Here, :math:`u^{\prime}|u^{\prime}|` is to be computed at time :math:`t_n`. The idea
is now to introduce
a *geometric mean*, defined by

.. math::
         (w^2)^n \approx w^{n-\frac{1}{2}}w^{n+\frac{1}{2}},

for some quantity :math:`w` depending on time. The error in the geometric mean
approximation is :math:`{\mathcal{O}(\Delta t^2)}`, the same as in the
approximation :math:`u^{\prime\prime}\approx D_tD_tu`. With :math:`w=u^{\prime}` it follows
that

.. math::
         [u^{\prime}|u^{\prime}|]^n \approx u^{\prime}(t_{n+\frac{1}{2}})|u^{\prime}(t_{n-\frac{1}{2}})|{\thinspace .}

The next step is to approximate
:math:`u^{\prime}` at :math:`t_{n\pm 1/2}`, and fortunately a centered difference
fits perfectly into the formulas since it involves :math:`u` values at
the mesh points only. With the approximations

.. math::
   :label: vib:ode2:quad:idea2
        
        u^{\prime}(t_{n+1/2})\approx [D_t u]^{n+\frac{1}{2}},\quad u^{\prime}(t_{n-1/2})\approx [D_t u]^{n-\frac{1}{2}},
        
        

we get

.. math::
   :label: _auto23
        
        [u^{\prime}|u^{\prime}|]^n
        \approx [D_tu]^{n+\frac{1}{2}}|[D_tu]^{n-\frac{1}{2}}| = \frac{u^{n+1}-u^n}{\Delta t}
        \frac{|u^n-u^{n-1}|}{\Delta t}
        {\thinspace .}
        
        

The counterpart to :eq:`vib:ode2:step3b` is then

.. math::
   :label: vib:ode2:step3b:quad
        
        m\frac{u^{n+1}-2u^n + u^{n-1}}{\Delta t^2}
        + b\frac{u^{n+1}-u^n}{\Delta t}\frac{|u^n-u^{n-1}|}{\Delta t}
        + s(u^n) = F^n,
        
        

which is linear in the unknown :math:`u^{n+1}`. Therefore, we can easily solve
:eq:`vib:ode2:step3b:quad`
with respect to :math:`u^{n+1}` and achieve the explicit updating formula

.. math::
        
        u^{n+1} =  \left( m + b|u^n-u^{n-1}|\right)^{-1}\times \nonumber
        

.. math::
   :label: vib:ode2:step4:quad
          
         \qquad \left(2m u^n - mu^{n-1} + bu^n|u^n-u^{n-1}| + \Delta t^2 (F^n - s(u^n))
        \right)
        {\thinspace .}
        
        

.. Make exercise to solve complicated u^1 equation with Bisection/Newton

In the derivation of a special equation for the first
time step we run into some trouble: inserting :eq:`vib:ode2:ic:du`
in :eq:`vib:ode2:step4:quad` for :math:`n=0` results in a complicated nonlinear
equation for :math:`u^1`. By thinking differently about the problem we can
easily get away with the nonlinearity again. We have for :math:`n=0` that
:math:`b[u^{\prime}|u^{\prime}|]^0 = bV|V|`. Using this value in :eq:`vib:ode2:quad:idea1`
gives

.. math::
   :label: _auto24
        
        [mD_tD_t u + bV|V| + s(u) = F]^0
        {\thinspace .}
        
        

Writing this equation out and using :eq:`vib:ode2:ic:du` results in the
special equation for the first time step:

.. math::
   :label: vib:ode2:step4b:quad
        
        u^1 = u^0 + \Delta t V + \frac{\Delta t^2}{2m}\left(-bV|V| - s(u^0) + F^0\right)
        {\thinspace .}
        
        

A forward-backward discretization of the quadratic damping term
---------------------------------------------------------------

The previous section first proposed to discretize the quadratic
damping term :math:`|u^{\prime}|u^{\prime}` using centered differences:
:math:`[|D_{2t}|D_{2t}u]^n`. As this gives rise to a nonlinearity in
:math:`u^{n+1}`, it was instead proposed to use a geometric mean combined
with centered differences.  But there are other alternatives. To get
rid of the nonlinearity in :math:`[|D_{2t}|D_{2t}u]^n`, one can think
differently: apply a backward difference to :math:`|u^{\prime}|`, such that
the term involves known values, and apply a forward difference to
:math:`u^{\prime}` to make the term linear in the unknown :math:`u^{n+1}`. With
mathematics,

.. math::
   :label: vib:ode2:nonlin:fbdiff
        
        [\beta |u^{\prime}|u^{\prime}]^n \approx \beta |[D_t^-u]^n|[D_t^+ u]^n =
        \beta\left\vert\frac{u^n-u^{n-1}}{\Delta t}\right\vert
        \frac{u^{n+1}-u^n}{\Delta t}{\thinspace .}
        
        

The forward and backward differences have both an error proportional
to :math:`\Delta t` so one may think the discretization above leads to
a first-order scheme.
However, by looking at the formulas, we realize that the forward-backward
differences in :eq:`vib:ode2:nonlin:fbdiff`
result in exactly the same scheme as in
:eq:`vib:ode2:step3b:quad` where we
used a geometric mean and centered differences and committed errors
of size :math:`{\mathcal{O}(\Delta t^2)}`. Therefore, the forward-backward
differences in :eq:`vib:ode2:nonlin:fbdiff`
act in a symmetric way and actually produce a second-order
accurate discretization of the quadratic damping term.

.. _vib:ode2:solver:

Implementation          (3)
---------------------------

The algorithm arising from the methods in the sections :ref:`vib:ode2:fdm:flin`
and :ref:`vib:ode2:fdm:fquad` is very similar to the undamped case in
the section :ref:`vib:ode1:fdm`. The difference is
basically a question of different formulas for :math:`u^1` and
:math:`u^{n+1}`. This is actually quite remarkable. The equation
:eq:`vib:ode2` is normally impossible to solve by pen and paper, but
possible for some special choices of :math:`F`, :math:`s`, and :math:`f`. On the
contrary, the complexity of the
nonlinear generalized model :eq:`vib:ode2` versus the
simple undamped model is not a big deal when we solve the
problem numerically!

The computational algorithm takes the form

 1. :math:`u^0=I`

 2. compute :math:`u^1` from :eq:`vib:ode2:step4b` if linear
    damping or :eq:`vib:ode2:step4b:quad` if quadratic damping

 3. for :math:`n=1,2,\ldots,N_t-1`:

   1. compute :math:`u^{n+1}` from :eq:`vib:ode2:step4` if linear
      damping or :eq:`vib:ode2:step4:quad` if quadratic damping

Modifying the ``solver`` function for the undamped case is fairly
easy, the big difference being many more terms and if tests on
the type of damping:

.. code-block:: python

        def solver(I, V, m, b, s, F, dt, T, damping='linear'):
            """
            Solve m*u'' + f(u') + s(u) = F(t) for t in (0,T],
            u(0)=I and u'(0)=V,
            by a central finite difference method with time step dt.
            If damping is 'linear', f(u')=b*u, while if damping is
            'quadratic', f(u')=b*u'*abs(u').
            F(t) and s(u) are Python functions.
            """
            dt = float(dt); b = float(b); m = float(m) # avoid integer div.
            Nt = int(round(T/dt))
            u = np.zeros(Nt+1)
            t = np.linspace(0, Nt*dt, Nt+1)
        
            u[0] = I
            if damping == 'linear':
                u[1] = u[0] + dt*V + dt**2/(2*m)*(-b*V - s(u[0]) + F(t[0]))
            elif damping == 'quadratic':
                u[1] = u[0] + dt*V + \ 
                       dt**2/(2*m)*(-b*V*abs(V) - s(u[0]) + F(t[0]))
        
            for n in range(1, Nt):
                if damping == 'linear':
                    u[n+1] = (2*m*u[n] + (b*dt/2 - m)*u[n-1] +
                              dt**2*(F(t[n]) - s(u[n])))/(m + b*dt/2)
                elif damping == 'quadratic':
                    u[n+1] = (2*m*u[n] - m*u[n-1] + b*u[n]*abs(u[n] - u[n-1])
                              + dt**2*(F(t[n]) - s(u[n])))/\ 
                              (m + b*abs(u[n] - u[n-1]))
            return u, t

The complete code resides in the file `vib.py <http://tinyurl.com/nm5587k/vib/vib.py>`__.

.. _vib:ode2:verify:

Verification          (2)
-------------------------

Constant solution
~~~~~~~~~~~~~~~~~

For debugging and initial verification, a constant solution is often
very useful. We choose :math:`{u_{\small\mbox{e}}}(t)=I`, which implies :math:`V=0`.
Inserted in the ODE, we get
:math:`F(t)=s(I)` for any choice of :math:`f`. Since the discrete derivative
of a constant vanishes (in particular, :math:`[D_{2t}I]^n=0`,
:math:`[D_tI]^n=0`, and :math:`[D_tD_t I]^n=0`), the constant solution also fulfills
the discrete equations. The constant should therefore be reproduced
to machine precision. The function ``test_constant`` in ``vib.py``
implements this test.

[**hpl 8**: Add verification tests for constant, linear, quadratic. Check how many bugs that are caught by these tests.]

Linear solution
~~~~~~~~~~~~~~~

Now we choose a linear solution: :math:`{u_{\small\mbox{e}}} = ct + d`. The initial condition
:math:`u(0)=I` implies :math:`d=I`, and :math:`u^{\prime}(0)=V` forces :math:`c` to be :math:`V`.
Inserting :math:`{u_{\small\mbox{e}}}=Vt+I` in the ODE with linear damping results in

.. math::
         0 + bV + s(Vt+I) = F(t),

while quadratic damping requires the source term

.. math::
         0 + b|V|V + s(Vt+I) = F(t){\thinspace .}

Since the finite difference approximations used to compute :math:`u^{\prime}` all
are exact for a linear function, it turns out that the linear :math:`{u_{\small\mbox{e}}}`
is also a solution of the discrete equations.
:ref:`vib:exer:verify:gen:linear` asks you to carry out
all the details.

Quadratic solution
~~~~~~~~~~~~~~~~~~

Choosing :math:`{u_{\small\mbox{e}}} = bt^2 + Vt + I`, with :math:`b` arbitrary,
fulfills the initial conditions and
fits the ODE if :math:`F` is adjusted properly. The solution also solves
the discrete equations with linear damping. However, this quadratic
polynomial in :math:`t` does not fulfill the discrete equations in case
of quadratic damping, because the geometric mean used in the approximation
of this term introduces an error.
Doing :ref:`vib:exer:verify:gen:linear` will reveal
the details. One can fit :math:`F^n` in the discrete equations such that
the quadratic polynomial is reproduced by the numerical method (to
machine precision).

.. More: classes, cases with pendulum approx u vs sin(u),

.. making UI via parampool

.. _vib:ode2:viz:

Visualization
-------------

The functions for visualizations differ significantly from
those in the undamped case in the ``vib_undamped.py`` program because,
in the present general case, we do not have an exact solution to
include in the plots. Moreover, we have no good estimate of
the periods of the oscillations as there will be one period
determined by the system parameters, essentially the
approximate frequency :math:`\sqrt{s'(0)/m}` for linear :math:`s` and small damping,
and one period dictated by :math:`F(t)` in case the excitation is periodic.
This is, however,
nothing that the program can depend on or make use of.
Therefore, the user has to specify :math:`T` and the window width
to get a plot that moves with the graph and shows
the most recent parts of it in long time simulations.

The ``vib.py`` code
contains several functions for analyzing the time series signal
and for visualizing the solutions.

.. _vib:ode2:ui:

User interface
--------------

.. index:: ArgumentParser (Python class)

.. index:: argparse (Python module)

The ``main`` function is changed substantially from
the ``vib_undamped.py`` code, since we need to
specify the new data :math:`c`, :math:`s(u)`, and :math:`F(t)`.  In addition, we must
set :math:`T` and the plot window width (instead of the number of periods we
want to simulate as in ``vib_undamped.py``). To figure out whether we
can use one plot for the whole time series or if we should follow the
most recent part of :math:`u`, we can use the ``plot_empricial_freq_and_amplitude``
function's estimate of the number of local maxima. This number is now
returned from the function and used in ``main`` to decide on the
visualization technique.

.. code-block:: python

        def main():
            import argparse
            parser = argparse.ArgumentParser()
            parser.add_argument('--I', type=float, default=1.0)
            parser.add_argument('--V', type=float, default=0.0)
            parser.add_argument('--m', type=float, default=1.0)
            parser.add_argument('--c', type=float, default=0.0)
            parser.add_argument('--s', type=str, default='u')
            parser.add_argument('--F', type=str, default='0')
            parser.add_argument('--dt', type=float, default=0.05)
            parser.add_argument('--T', type=float, default=140)
            parser.add_argument('--damping', type=str, default='linear')
            parser.add_argument('--window_width', type=float, default=30)
            parser.add_argument('--savefig', action='store_true')
            a = parser.parse_args()
            from scitools.std import StringFunction
            s = StringFunction(a.s, independent_variable='u')
            F = StringFunction(a.F, independent_variable='t')
            I, V, m, c, dt, T, window_width, savefig, damping = \ 
               a.I, a.V, a.m, a.c, a.dt, a.T, a.window_width, a.savefig, \ 
               a.damping
        
            u, t = solver(I, V, m, c, s, F, dt, T)
            num_periods = empirical_freq_and_amplitude(u, t)
            if num_periods <= 15:
                figure()
                visualize(u, t)
            else:
                visualize_front(u, t, window_width, savefig)
            show()

The program ``vib.py`` contains
the above code snippets and can solve the model problem
:eq:`vib:ode2`. As a demo of ``vib.py``, we consider the case
:math:`I=1`, :math:`V=0`, :math:`m=1`, :math:`c=0.03`, :math:`s(u)=\sin(u)`, :math:`F(t)=3\cos(4t)`,
:math:`\Delta t = 0.05`, and :math:`T=140`. The relevant command to run is

.. code-block:: text

        Terminal> python vib.py --s 'sin(u)' --F '3*cos(4*t)' --c 0.03

This results in a `moving window following the function <http://tinyurl.com/opdfafk/pub/mov-vib/vib_generalized_dt0.05/index.html>`__ on the screen.
Figure :ref:`vib:ode2:fig:demo` shows a part of the time series.

.. _vib:ode2:fig:demo:

.. figure:: fig-vib/vib_gen_demo.png
   :width: 600

   *Damped oscillator excited by a sinusoidal function*

The Euler-Cromer scheme for the generalized model
-------------------------------------------------

The ideas of the Euler-Cromer method from the section :ref:`vib:model2x2:EulerCromer`
carry over to the generalized model. We write :eq:`vib:ode2`
as two equations for :math:`u` and :math:`v=u^{\prime}`. The first equation is taken as the
one with :math:`v'` on the left-hand side:

.. math::
   :label: vib:ode2:EulerCromer:veq
        
        v' = \frac{1}{m}(F(t)-s(u)-f(v)),
        
        

.. math::
   :label: vib:ode2:EulerCromer:ueq
          
        u^{\prime} = v{\thinspace .}
        
        

The idea is to step :eq:`vib:ode2:EulerCromer:veq` forward using
a standard Forward Euler method, while we update :math:`u` from
:eq:`vib:ode2:EulerCromer:ueq` with a Backward Euler method,
utilizing the recent, computed :math:`v^{n+1}` value. In detail,

.. math::
   :label: vib:ode2:EulerCromer:dveq0a
        
        \frac{v^{n+1}-v^n}{\Delta t} = \frac{1}{m}(F(t_n)-s(u^n)-f(v^n)),
        
        

.. math::
   :label: vib:ode2:EulerCromer:dueq0a
          
        \frac{u^{n+1}-u^n}{\Delta t} = v^{n+1},
        
        

resulting in the explicit scheme

.. math::
   :label: vib:ode2:EulerCromer:dveq
        
        v^{n+1} = v^n + \Delta t\frac{1}{m}(F(t_n)-s(u^n)-f(v^n)),
        
        

.. math::
   :label: vib:ode2:EulerCromer:dueq0
          
        u^{n+1} = u^n + \Delta t\,v^{n+1}{\thinspace .}
        
        

We immediately note one very favorable feature of this scheme: all the
nonlinearities in :math:`s(u)` and :math:`f(v)` are evaluated at a previous time
level. This makes the Euler-Cromer method easier to apply and
hence much more convenient than the centered scheme for the second-order
ODE :eq:`vib:ode2`.

The initial conditions are trivially set as

.. math::
   :label: _auto25
        
        v^0 = V,
        
        

.. math::
   :label: _auto26
          
        u^0 = I{\thinspace .}
        
        

[**hpl 9**: odespy for the generalized problem]

Exercises and Problems
======================

.. --- begin exercise ---

.. _vib:exer:undamped:verify:linquad:

Problem 1: Use linear/quadratic functions for verification
----------------------------------------------------------

Consider the ODE problem

.. math::
         u^{\prime\prime} + \omega^2u=f(t), \quad u(0)=I,\ u^{\prime}(0)=V,\ t\in(0,T]{\thinspace .}

Discretize this equation according to
:math:`[D_tD_t u + \omega^2 u = f]^n`.

**a)**
Derive the equation for the
first time step (:math:`u^1`).

**b)**
For verification purposes,
we use the method of manufactured solutions (MMS) with the
choice of :math:`{u_{\small\mbox{e}}}(x,t)= ct+d`.
Find restrictions on :math:`c` and :math:`d` from
the initial conditions. Compute the corresponding source term :math:`f` by term.
Show that :math:`[D_tD_t t]^n=0` and use the fact
that the :math:`D_tD_t` operator is linear,
:math:`[D_tD_t (ct+d)]^n = c[D_tD_t t]^n + [D_tD_t d]^n = 0`, to show that
:math:`{u_{\small\mbox{e}}}` is also a perfect solution of the discrete equations.

**c)**
Use ``sympy`` to do the symbolic calculations above. Here is a
sketch of the program ``vib_undamped_verify_mms.py``:

.. code-block:: python

        import sympy as sym
        V, t, I, w, dt = sym.symbols('V t I w dt')  # global symbols
        f = None  # global variable for the source term in the ODE
        
        def ode_source_term(u):
            """Return the terms in the ODE that the source term
            must balance, here u'' + w**2*u.
            u is symbolic Python function of t."""
            return sym.diff(u(t), t, t) + w**2*u(t)
        
        def residual_discrete_eq(u):
            """Return the residual of the discrete eq. with u inserted."""
            R = ...
            return sym.simplify(R)
        
        def residual_discrete_eq_step1(u):
            """Return the residual of the discrete eq. at the first
            step with u inserted."""
            R = ...
            return sym.simplify(R)
        
        def DtDt(u, dt):
            """Return 2nd-order finite difference for u_tt.
            u is a symbolic Python function of t.
            """
            return ...
        
        def main(u):
            """
            Given some chosen solution u (as a function of t, implemented
            as a Python function), use the method of manufactured solutions
            to compute the source term f, and check if u also solves
            the discrete equations.
            """
            print '=== Testing exact solution: %s ===' % u
            print "Initial conditions u(0)=%s, u'(0)=%s:" % \ 
                  (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))
        
            # Method of manufactured solution requires fitting f
            global f  # source term in the ODE
            f = sym.simplify(ode_lhs(u))
        
            # Residual in discrete equations (should be 0)
            print 'residual step1:', residual_discrete_eq_step1(u)
            print 'residual:', residual_discrete_eq(u)
        
        def linear():
            main(lambda t: V*t + I)
        
        if __name__ == '__main__':
            linear()

Fill in the various functions such that the calls in the ``main``
function works.

**d)**
The purpose now is to choose a quadratic function
:math:`{u_{\small\mbox{e}}} = bt^2 + ct + d` as exact solution. Extend the ``sympy``
code above with a function ``quadratic`` for fitting ``f`` and checking
if the discrete equations are fulfilled. (The function is very similar
to ``linear``.)

.. Check with hand calculations that the ``sympy`` implementation

.. is correct.

**e)**
Will a polynomial of degree three fulfill the discrete equations?

**f)**
Implement a ``solver`` function for computing the numerical
solution of this problem.

**g)**
Write a nose test for checking that the quadratic solution
is computed to correctly (too machine precision, but the
round-off errors accumulate and increase with :math:`T`) by the ``solver``
function.

Filename: ``vib_undamped_verify_mms``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:phase:err:growth:

Exercise 2: Show linear growth of the phase with time
-----------------------------------------------------

Consider an exact solution :math:`I\cos (\omega t)` and an
approximation :math:`I\cos(\tilde\omega t)`.
Define the phase error as time lag between the peak :math:`I`
in the exact solution and the corresponding peak in the approximation
after :math:`m` periods of oscillations. Show that this phase error
is linear in :math:`m`.
Filename: ``vib_phase_error_growth``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:w:adjust:

Exercise 3: Improve the accuracy by adjusting the frequency
-----------------------------------------------------------

According to :eq:`vib:ode1:tildeomega:series`, the numerical
frequency deviates from the exact frequency by a (dominating) amount
:math:`\omega^3\Delta t^2/24 >0`. Replace the ``w`` parameter in the algorithm
in the ``solver`` function in ``vib_undamped.py`` by ``w*(1 -
(1./24)*w**2*dt**2`` and test how this adjustment in the numerical
algorithm improves the accuracy (use :math:`\Delta t =0.1` and simulate
for 80 periods, with and without adjustment of :math:`\omega`).
Filename: ``vib_adjust_w``.

.. How does this go if

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:undamped:adaptive:

Exercise 4: See if adaptive methods improve the phase error
-----------------------------------------------------------

Adaptive methods for solving ODEs aim at adjusting :math:`\Delta t` such
that the error is within a user-prescribed tolerance. Implement the
equation :math:`u^{\prime\prime}+u=0` in the `Odespy <https://github.com/hplgit/odespy>`__
software. Use the example `on adaptive
schemes <http://hplgit.github.io/decay-book/doc/pub/book/sphinx/._book006.html#example-adaptive-runge-kutta-methods>`__
in [Ref1]_.  Run the scheme with a very low
tolerance (say :math:`10^{-14}`) and for a long time, check the number of
time points in the solver's mesh (``len(solver.t_all)``), and compare
the phase error with that produced by the simple finite difference
method from the section :ref:`vib:ode1:fdm` with the same number of (equally
spaced) mesh points. The question is whether it pays off to use an
adaptive solver or if equally many points with a simple method gives
about the same accuracy.
Filename: ``vib_undamped_adaptive``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:step4b:alt:

Exercise 5: Use a Taylor polynomial to compute :math:`u^1`
----------------------------------------------------------

As an alternative to the derivation of :eq:`vib:ode1:step4b` for
computing :math:`u^1`, one can use a Taylor polynomial with three terms
for :math:`u^1`:

.. math::
         u(t_1) \approx u(0) + u^{\prime}(0)\Delta t + {\frac{1}{2}}u^{\prime\prime}(0)\Delta t^2

With :math:`u^{\prime\prime}=-\omega^2 u` and :math:`u^{\prime}(0)=0`, show that this method also leads to
:eq:`vib:ode1:step4b`. Generalize the condition on :math:`u^{\prime}(0)` to
be :math:`u^{\prime}(0)=V` and compute :math:`u^1` in this case with both methods.
Filename: ``vib_first_step``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:wdt:limit:

Exercise 6: Find the minimal resolution of an oscillatory function
------------------------------------------------------------------

.. Short: Find the largest relevant value of :math:`\omega\Delta t`

Sketch the function on a given mesh which has the highest possible
frequency. That is, this oscillatory "cos-like" function has its
maxima and minima at every two grid points.  Find an expression for
the frequency of this function, and use the result to find the largest
relevant value of :math:`\omega\Delta t` when :math:`\omega` is the frequency
of an oscillating function and :math:`\Delta t` is the mesh spacing.
Filename: ``vib_largest_wdt``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:fd:exp:plot:

Exercise 7: Visualize the accuracy of finite differences for a cosine function
------------------------------------------------------------------------------

.. Short: Visualize the accuracy of finite differences

We introduce the error fraction

.. math::
         E = \frac{[D_tD_t u]^n}{u^{\prime\prime}(t_n)} 

to measure the error in the finite difference approximation :math:`D_tD_tu` to
:math:`u^{\prime\prime}`.
Compute :math:`E`
for the specific choice of a cosine/sine function of the
form :math:`u=\exp{(i\omega t)}` and show that

.. math::
         E = \left(\frac{2}{\omega\Delta t}\right)^2
        \sin^2(\frac{\omega\Delta t}{2})
        {\thinspace .}
        

Plot :math:`E` as a function of :math:`p=\omega\Delta t`. The relevant
values of :math:`p` are :math:`[0,\pi]` (see :ref:`vib:exer:wdt:limit`
for why :math:`p>\pi` does not make sense).
The deviation of the curve from unity visualizes the error in the
approximation. Also expand :math:`E` as a Taylor polynomial in :math:`p` up to
fourth degree (use, e.g., ``sympy``).
Filename: ``vib_plot_fd_exp_error``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:energy:convrate:

Exercise 8: Verify convergence rates of the error in energy
-----------------------------------------------------------

We consider the ODE problem :math:`u^{\prime\prime} + \omega^2u=0`, :math:`u(0)=I`, :math:`u^{\prime}(0)=V`,
for :math:`t\in (0,T]`. The total energy of the solution
:math:`E(t)=\frac{1}{2}(u^{\prime})^2 + \frac{1}{2}\omega^2 u^2` should stay
constant.
The error in energy can be computed as explained in
the section :ref:`vib:model1:energy`.

Make a nose test in a file ``test_error_conv.py``, where code from
``vib_undamped.py`` is imported, but the ``convergence_rates`` and
``test_convergence_rates`` functions are copied and modified to also
incorporate computations of the error in energy and the convergence
rate of this error. The expected rate is 2.
Filename: ``test_error_conv``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:verify:gen:linear:

Exercise 9: Use linear/quadratic functions for verification
-----------------------------------------------------------

This exercise is a generalization of :ref:`vib:exer:undamped:verify:linquad` to the extended model problem
:eq:`vib:ode2` where the damping term is either linear or quadratic.
Solve the various subproblems and see how the results and problem
settings change with the generalized ODE in case of linear or
quadratic damping. By modifying the code from :ref:`vib:exer:undamped:verify:linquad`, ``sympy`` will do most
of the work required to analyze the generalized problem.
Filename: ``vib_verify_mms``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:discrete:omega:

Exercise 10: Use an exact discrete solution for verification
------------------------------------------------------------

Write a nose test function in a separate file
that employs the exact discrete solution
:eq:`vib:ode1:un:exact` to verify the implementation of the
``solver`` function in the file ``vib_undamped.py``.
Filename: ``test_vib_undamped_exact_discrete_sol``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:conv:rate:

Exercise 11: Use analytical solution for convergence rate tests
---------------------------------------------------------------

The purpose of this exercise is to perform convergence tests of the
problem :eq:`vib:ode2` when :math:`s(u)=\omega^2u` and :math:`F(t)=A\sin\phi t`.
Find the complete analytical solution to the problem in this case
(most textbooks on mechanics or ordinary differential equations list
the various elements you need to write down the exact solution).
Modify the ``convergence_rate`` function from the ``vib_undamped.py``
program to perform experiments with the extended model.  Verify that
the error is of order :math:`\Delta t^2`.
Filename: ``vib_conv_rate``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:undamped:odespy:

Exercise 12: Investigate the amplitude errors of many solvers
-------------------------------------------------------------

Use the program ``vib_undamped_odespy.py`` from the section :ref:`vib:model2x2:compare` and the amplitude estimation from the
``amplitudes`` function in the ``vib_undamped.py`` file (see the section :ref:`vib:ode1:empirical`) to investigate how well famous methods for
1st-order ODEs can preserve the amplitude of :math:`u` in undamped
oscillations.  Test, for example, the 3rd- and 4th-order Runge-Kutta
methods (``RK3``, ``RK4``), the Crank-Nicolson method (``CrankNicolson``),
the 2nd- and 3rd-order Adams-Bashforth methods (``AdamsBashforth2``,
``AdamsBashforth3``), and a 2nd-order Backwards scheme
(``Backward2Step``).  The relevant governing equations are listed in
the beginning of the section :ref:`vib:model2x2`.
Filename: ``vib_amplitude_errors``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:memsave:

Exercise 13: Minimize memory usage of a vibration solver
--------------------------------------------------------

The program `vib.py <http://tinyurl.com/nm5587k/vib/vib.py>`__
store the complete solution :math:`u^0,u^1,\ldots,u^{N_t}` in memory, which is
convenient for later plotting.
Make a memory minimizing version of this program where only the last three
:math:`u^{n+1}`, :math:`u^n`, and :math:`u^{n-1}` values are stored in memory.
Write each computed :math:`(t_{n+1}, u^{n+1})` pair to file.
Visualize the data in the file (a cool solution is to
read one line at a time and
plot the :math:`u` value using the line-by-line plotter in the
``visualize_front_ascii`` function - this technique makes it trivial
to visualize very long time simulations).
Filename: ``vib_memsave``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:gen:class:

Exercise 14: Implement the solver via classes
---------------------------------------------

Reimplement the ``vib.py``
program
using a class ``Problem`` to hold all the physical parameters of the problem,
a class ``Solver`` to hold the numerical parameters and compute the
solution, and a class ``Visualizer`` to display the solution.

.. --- begin hint in exercise ---

**Hint.**
Use the ideas and examples
for an `ODE model <http://hplgit.github.io/decay-book/doc/pub/book/sphinx/._book009.html#classes-for-problem-and-solution-method>`__ in [Ref1]_.
More specifically, make a superclass ``Problem`` for holding the scalar
physical parameters of a problem and let subclasses implement the
:math:`s(u)` and :math:`F(t)` functions as methods.
Try to call up as much existing functionality in ``vib.py`` as possible.

.. --- end hint in exercise ---

Filename: ``vib_class``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:DtDt:asDtpDtm:

Exercise 15: Interpret :math:`[D_tD_t u]^n` as a forward-backward difference
----------------------------------------------------------------------------

Show that the difference :math:`[D_t D_tu]^n` is equal to :math:`[D_t^+D_t^-u]^n`
and :math:`D_t^-D_t^+u]^n`. That is, instead of applying a centered difference
twice one can alternatively apply a mixture forward and backward
differences.
Filename: ``vib_DtDt_fw_bw``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:quad:damping:bw:

Exercise 16: Use a backward difference for the damping term
-----------------------------------------------------------

As an alternative to discretizing the damping terms :math:`\beta u^{\prime}` and
:math:`\beta |u^{\prime}|u^{\prime}` by centered differences, we may apply
backward differences:

.. math::
        
        [u^{\prime}]^n &\approx [D_t^-u]^n,\\ 
        & [|u^{\prime}|u^{\prime}]^n &\approx [|D_t^-u|D_t^-u]^n
        = |[D_t^-u]^n|[D_t^-u]^n{\thinspace .}
        

The advantage of the backward difference is that the damping term is
evaluated using known values :math:`u^n` and :math:`u^{n-1}` only.
Extend the `vib.py <http://tinyurl.com/nm5587k/vib/vib.py>`__ code with a scheme based
on using backward differences in the damping terms. Add statements
to compare the original approach with centered difference and the
new idea launched in this exercise. Perform numerical experiments
to investigate how much accuracy that is lost by using the backward
differences.
Filename: ``vib_gen_bwdamping``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:EulerCromer:analysis:

Exercise 17: Analysis of the Euler-Cromer scheme
------------------------------------------------

The Euler-Cromer scheme for the model problem
:math:`u^{\prime\prime} + \omega^2 u =0`, :math:`u(0)=I`, :math:`u^{\prime}(0)=0`, is given in
:eq:`vib:model2x2:EulerCromer:ueq1b`-:eq:`vib:model2x2:EulerCromer:veq1b`.
Find the exact discrete solutions of this scheme and show that the solution
for :math:`u^n` coincides with that found in the section :ref:`vib:ode1:analysis`.

.. --- begin hint in exercise ---

**Hint.**
Use an "ansatz" :math:`u^n=I\exp{(i\tilde\omega\Delta t\,n)}` and
:math:`v^n=qu^n`, where :math:`\tilde\omega` and :math:`q` are unknown parameters. The
following formula is handy:

.. math::
         \boldsymbol{e}^{i\tilde\omega\Delta t} + e^{i\tilde\omega(-\Delta t)} - 2
        = 2\left(\cosh(i\tilde\omega\Delta t) -1 \right)
        =-4\sin^2(\frac{\tilde\omega\Delta t}{2}){\thinspace .}

.. --- end hint in exercise ---

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

.. mech systems: horizontal, vertical/hanging

.. box with mu*M*g*v/|v| friction force, treat nonlinearity with geometric mean

.. pendulum

.. elastic pendulum

.. bouncing ball (just move text from exercise)

.. bumpy road

.. moored ship

.. electrical circuits, see ode2.p.tex

.. 0D blood flow?

.. waves: 1D blood flow

.. general particle laws and velocity verlet, make exercises

.. see `<http://en.wikipedia.org/wiki/Velocity_Verlet>`_

.. --- end exercise ---

.. _vib:app:

Applications of vibration models
================================

The following text derives some of the most well-known physical problems
that lead to
second-order ODE models of the type addressed in this document.
We consider a simple spring-mass system; thereafter extended with
nonlinear spring, damping, and external excitation; a spring-mass system
with sliding friction; a simple and a physical (classical) pendulum;
and an elastic pendulum.

.. _vib:app:mass_spring:

Oscillating mass attached to a spring
-------------------------------------

.. _vib:app:mass_spring:fig:

.. figure:: fig-vib/oscillator_spring.png
   :width: 500

   *Simple oscillating mass*

The most fundamental mechanical vibration system is depicted in Figure
:ref:`vib:app:mass_spring:fig`. A body with mass :math:`m` is attached to a
spring and can move horizontally without friction (in the wheels). The
position of the body is given by the vector :math:`\boldsymbol{r}(t) = u(t)\boldsymbol{i}`, where
:math:`\boldsymbol{i}` is a unit vector in :math:`x` direction.
There is
only one force acting on the body: a spring force :math:`\boldsymbol{F}_s =-ku\boldsymbol{i}`, where
:math:`k` is a constant. The point :math:`x=0`, where :math:`u=0`, must therefore
correspond to the body's position
where the spring is neither extended nor compressed, so the force
vanishes.

The basic physical principle that governs the motion of the body is
Newton's second law of motion: :math:`\boldsymbol{F}=m\boldsymbol{a}`, where
:math:`\boldsymbol{F}` is the sum of forces on the body, :math:`m` is its mass, and :math:`\boldsymbol{a}=\ddot\boldsymbol{r}`
is the acceleration. We use the dot for differentiation with respect
to time, which is
usual in mechanics. Newton's second law simplifies here
to :math:`-\boldsymbol{F}_s=m\ddot u\boldsymbol{i}`, which translates to

.. math::
         -ku = m\ddot u{\thinspace .}

Two initial conditions are needed: :math:`u(0)=I`, :math:`\dot u(0)=V`.
The ODE problem is normally written as

.. math::
   :label: vib:app:mass_spring:eqx
        
        m\ddot u + ku = 0,\quad u(0)=I,\ \dot u(0)=V{\thinspace .}
        
        

It is
not uncommon to divide by :math:`m`
and introduce the frequency :math:`\omega = \sqrt{k/m}`:

.. math::
   :label: vib:app:mass_spring:equ
        
        \ddot u + \omega^2 u = 0,\quad u(0)=I,\  \dot u(0)=V{\thinspace .}
        
        

This is the model problem in the first part of this chapter, with the
small difference that we write the time derivative of :math:`u` with a dot
above, while we used :math:`u^{\prime}` and :math:`u^{\prime\prime}` in previous
parts of the document.

.. index:: DOF (degree of freedom)

Since only one scalar mathematical quantity, :math:`u(t)`, describes the
complete motion, we say that the mechanical system has one degree of freedom
(DOF).

Scaling          (1)
~~~~~~~~~~~~~~~~~~~~

For numerical simulations it is very convenient to scale :eq:`vib:app:mass_spring:equ` and thereby get rid of the problem of finding relevant values
for all the parameters :math:`m`, :math:`k`, :math:`I`, and :math:`V`.
Since the amplitude of the oscillations are dictated by :math:`I` and :math:`V`
(or more precisely, :math:`V/\omega`), we scale :math:`u` by :math:`I` (or :math:`V/omega` if
:math:`I=0`):

.. math::
         \bar u = \frac{u}{I},\quad \bar t = \frac{t}{t_c}{\thinspace .}

The time scale :math:`t_c` is normally chosen as the inverse period :math:`2\pi/\omega` or
angular frequency :math:`1/\omega`, most often as :math:`t_c=1/\omega`.
Inserting the dimensionless quantities :math:`\bar u` and :math:`\bar t` in
:eq:`vib:app:mass_spring:equ` results in the scaled problem

.. math::
         \frac{d^2\bar u}{d\bar t^2} + \bar u = 0,\quad \bar u(0)=1,\ \frac{\bar u}{\bar t}(0)=\beta = \frac{V}{I\omega},

where :math:`\beta` is a dimensionless number. Any motion that starts from rest
(:math:`V=0`) is free of parameters in the scaled model!

The physics
~~~~~~~~~~~

The typical physics of the system in Figure :ref:`vib:app:mass_spring:fig` can
be described as follows.
Initially, we displace the body to some position :math:`I`, say at rest
(:math:`V=0`). After releasing the body, the spring, which is extended, will
act with a force :math:`-kI\boldsymbol{i}` and pull the body to the left. This force
causes an acceleration and therefore increases velocity. The body passes
the point :math:`x=0`, where :math:`u=0`,
and the spring will then be compressed and act with a
force :math:`kx\boldsymbol{i}` against the motion and cause retardation. At some point,
the motion stops and the velocity is zero, before the spring force
:math:`kx\boldsymbol{i}` accelerates the body in positive direction. The result is that
the body accelerates back and forth. As long as there is no friction
forces to damp the motion, the oscillations will continue forever.

.. _vib:app:mass_gen:

General mechanical vibrating system
-----------------------------------

.. _vib:app:mass_gen:fig:

.. figure:: fig-vib/oscillator_general.png
   :width: 500

   *General oscillating system*

The mechanical system in Figure :ref:`vib:app:mass_spring:fig` can easily be
extended to the more general system in Figure :ref:`vib:app:mass_gen:fig`,
where the body is attached to a spring and a dashpot, and also subject
to an environmental force :math:`F(t)\boldsymbol{i}`. The system has still only one
degree of freedom since the body can only move back and forth parallel to
the :math:`x` axis. The spring force was linear, :math:`\boldsymbol{F}_s=-ku\boldsymbol{i}`,
in the section :ref:`vib:app:mass_spring`, but in more general cases it can
depend nonlinearly on the position. We therefore set :math:`\boldsymbol{F}_s=s(u)\boldsymbol{i}`.
The dashpot, which acts
as a damper, results in a force :math:`\boldsymbol{F}_d` that depends on the body's
velocity :math:`\dot u` and that always acts against the motion.
The mathematical model of the force is written :math:`\boldsymbol{F}_d =f(\dot u)\boldsymbol{i}`.
A positive :math:`\dot u` must result in a force acting in the positive :math:`x`
direction.
Finally, we have the external environmental force :math:`\boldsymbol{F}_e = F(t)\boldsymbol{i}`.

Newton's second law of motion now involves three forces:

.. math::
         F(t)\boldsymbol{i} + f(\dot u)ii - s(u)\boldsymbol{i} = m\ddot u \boldsymbol{i}{\thinspace .}

The common mathematical form of the ODE problem is

.. math::
   :label: vib:app:mass_gen:equ
        
        m\ddot u + f(\dot u) + s(u) = F(t),\quad u(0)=I,\ \dot u(0)=V{\thinspace .}
        
        

This is the generalized problem treated in the last part of the
present chapter, but with prime denoting the derivative instead of the dot.

The most common models for the spring and dashpot are linear: :math:`f(\dot u)
=b\dot u` with a constant :math:`b\geq 0`, and :math:`s(u)=ku` for a constant :math:`k`.

Scaling          (2)
~~~~~~~~~~~~~~~~~~~~

A specific scaling requires specific choices of :math:`f`, :math:`s`, and :math:`F`.
Suppose we have

.. math::
         f(\dot u) = b|\dot u|\dot u,\quad s(u)=ku,\quad F(t)=A\sin(\phi t){\thinspace .}

We introduce dimensionless variables as usual, :math:`\bar u = u/u_c` and
:math:`\bar t = t/t_c`. The scale :math:`u_c` depends both on the initial conditions
and :math:`F`, but as time grows, the effect of the initial conditions die out
and :math:`F` will drive the motion. Inserting :math:`\bar u` and :math:`\bar t` in the
ODE gives

.. math::
         m\frac{u_c}{t_c^2}\frac{d^2\bar u}{d\bar t^2}
        + b\frac{u_c^2}{t_c^2}\left\vert\frac{d\bar u}{d\bar t}\right\vert
        \frac{d\bar u}{d\bar t} + ku_c\bar u = A\sin(\phi t_c\bar t){\thinspace .}

We divide by :math:`u_c/t_c^2` and demand the coefficients of the
:math:`\bar u` and the forcing term from :math:`F(t)` to have unit coefficients.
This leads to the scales

.. math::
         t_c = \sqrt{\frac{m}{k}},\quad u_c = \frac{A}{k}{\thinspace .}

The scaled ODE becomes

.. math::
   :label: vib:app:mass_gen:scaled
        
        \frac{d^2\bar u}{d\bar t^2}
        + 2\beta\left\vert\frac{d\bar u}{d\bar t}\right\vert
        \frac{d\bar u}{d\bar t} + \bar u = \sin(\gamma\bar t),
        
        

where there are two dimensionless numbers:

.. math::
         \beta = \frac{Ab}{2mk},\quad\gamma =\phi\sqrt{\frac{m}{k}}{\thinspace .}

The :math:`\beta` number measures the size of the damping term (relative to unity)
and is assumed to be small, basically because :math:`b` is small. The :math:`\phi`
number is the ratio of the time scale of free vibrations and the time scale
of the forcing.
The scaled initial conditions have two other dimensionless numbers
as values:

.. math::
         \bar u(0) = \frac{Ik}{A},\quad \frac{d\bar u}{d\bar t}=\frac{t_c}{u_c}V = \frac{V}{A}\sqrt{mk}{\thinspace .}

.. _vib:app:mass_sliding:

A sliding mass attached to a spring
-----------------------------------

Consider a variant of the oscillating body in the section :ref:`vib:app:mass_spring`
and Figure :ref:`vib:app:mass_spring:fig`: the body rests on a flat
surface, and there is sliding friction between the body and the surface.
Figure :ref:`vib:app:mass_sliding:fig` depicts the problem.

.. _vib:app:mass_sliding:fig:

.. figure:: fig-vib/oscillator_sliding.png
   :width: 500

   *Sketch of a body sliding on a surface*

The body is attached to a spring with spring force :math:`-s(u)\boldsymbol{i}`.
The friction force is proportional to the normal force on the surface,
:math:`-mg\boldsymbol{j}`, and given by :math:`-f(\dot u)\boldsymbol{i}`, where

.. math::
         f(\dot u) = \left\lbrace\begin{array}{ll}
        -\mu mg,& \dot u < 0,\\ 
        \mu mg, & \dot u > 0,\\ 
        0,      & \dot u=0
        \end{array}\right.

Here, :math:`\mu` is a friction coefficient. With the signum function

.. math::
         \mbox{sign(x)} = \left\lbrace\begin{array}{ll}
        -1,& x < 0,\\ 
        1, & x > 0,\\ 
        0, & x=0
        \end{array}\right.

we can simply write :math:`f(\dot u) = \mu mg\,\hbox{sign}(\dot u)`
(the sign function is implemented by ``numpy.sign``).

The equation of motion becomes

.. math::
   :label: vib:app:mass_sliding:equ
        
        m\ddot u + \mu mg\hbox{sign}(\dot u) + s(u) = 0,\quad u(0)=I,\ \dot u(0)=V{\thinspace .}
        
        

.. _vib:app:washmach:

A jumping washing machine
-------------------------

A washing machine is placed on four springs with efficient dampers.
If the machine contains just a few clothes, the circular motion of
the machine induces a sinusoidal external force and the machine will
jump up and down if the frequency of the external force is close to
the natural frequency of the machine and its spring-damper system.

[**hpl 10**: Not finished. This is a good example on resonance.]

.. _vib:app:pendulum:

Motion of a pendulum
--------------------

A classical problem in mechanics is the motion of a pendulum. We first
consider a `simple pendulum <https://en.wikipedia.org/wiki/Pendulum>`__:
a small body of mass :math:`m` is attached to a massless wire and can oscillate back and forth
in the gravity field. Figure :ref:`vib:app:pendulum:fig_problem` shows
a sketch of the problem.

.. _vib:app:pendulum:fig_problem:

.. figure:: fig-vib/pendulum_problem.png
   :width: 300

   *Sketch of a simple pendulum*

The motion is governed by Newton's 2nd law, so we need to find expressions
for the forces and the acceleration. Three forces on the body are
considered: an unknown force :math:`S` from the wire, the gravity force :math:`mg`,
and an air resistance force, :math:`\frac{1}{2}C_D\varrho A |v|v`,
hereafter called the drag force,
directed against the velocity of the body. Here, :math:`C_D` is a drag coefficient,
:math:`\varrho` is the density of air, :math:`A` is the cross section area of the body,
and :math:`v` is the velocity.

We introduce a coordinate system with polar coordinates and unit
vectors :math:`{\boldsymbol{i}_r}` and :math:`\boldsymbol{i}_{\theta}` as shown in Figure :ref:`vib:app:pendulum:fig_forces`.
The position of the center of mass of the body is

.. math::
         \boldsymbol{r}(t) = x_0\boldsymbol{i} + y_0\boldsymbol{j} + L{\boldsymbol{i}_r},

where :math:`\boldsymbol{i}` and :math:`\boldsymbol{j}` are unit vectors in the corresponding Cartesian
coordinate system in the :math:`x` and :math:`y` directions, respectively. We have
that :math:`{\boldsymbol{i}_r} = \cos\theta\boldsymbol{i} +\sin\theta\boldsymbol{j}`.

.. _vib:app:pendulum:fig_forces:

.. figure:: fig-vib/pendulum_forces.png
   :width: 400

   *Forces acting on a simple pendulum*

The forces are now expressed as follows.

 * Wire force: :math:`-S{\boldsymbol{i}_r}`

 * Gravity force: :math:`-mg\boldsymbol{j} = mg(-\sin\theta\boldsymbol{i}_{\theta} + \cos\theta{\boldsymbol{i}_r})`

 * Drag force: :math:`-\frac{1}{2}C_D\varrho A |v|v\boldsymbol{i}_{\theta}`

Since a positive velocity means movement in the direction of :math:`\boldsymbol{i}_{\theta}`,
the drag force must be directed along :math:`-\boldsymbol{i}_{\theta}`.

The velocity of the body is found from :math:`\boldsymbol{r}`:

.. math::
         \boldsymbol{v}(t) = \dot\boldsymbol{r} (t) = \frac{d}{d\theta}(x_0\boldsymbol{i} + y_0\boldsymbol{j} + L{\boldsymbol{i}_r})\frac{d\theta}{dt} = L\dot\theta\boldsymbol{i}_{\theta},

since :math:`\frac{d}{d\theta}{\boldsymbol{i}_r} = \boldsymbol{i}_{\theta}`. It follows that :math:`v=|\boldsymbol{v}|=L\dot\theta`.
The acceleration is

.. math::
         \boldsymbol{a}(t) = \dot\boldsymbol{v}(r) = \frac{d}{dt}(L\dot\theta\boldsymbol{i}_{\theta})
        = L\ddot\theta\boldsymbol{i}_{\theta} + L\dot\theta\frac{d\boldsymbol{i}_{\theta}}{d\theta}\dot\theta =
        = L\ddot\theta\boldsymbol{i}_{\theta} - L\dot\theta^2{\boldsymbol{i}_r},

since :math:`\frac{d}{d\theta}\boldsymbol{i}_{\theta} = -{\boldsymbol{i}_r}`.

Newton's 2nd law of motion becomes

.. math::
         -S{\boldsymbol{i}_r} + mg(-\sin\theta\boldsymbol{i}_{\theta} + \cos\theta{\boldsymbol{i}_r}) -
        \frac{1}{2}C_D\varrho AL^2|\dot\theta|\dot\theta\boldsymbol{i}_{\theta}
        = mL\ddot\theta\dot\theta\boldsymbol{i}_{\theta} - L\dot\theta^2{\boldsymbol{i}_r},

leading to two component equations

.. math::
   :label: vib:app:pendulum:ir
        
        -S + mg\cos\theta = -L\dot\theta^2,
        
        

.. math::
   :label: vib:app:pendulum:ith
          
        -mg\sin\theta - \frac{1}{2}C_D\varrho AL^2|\dot\theta|\dot\theta
        = mL\ddot\theta{\thinspace .}
        
        

From :eq:`vib:app:pendulum:ir` we get an expression for
:math:`S=mg\cos\theta + L\dot\theta^2`, and from :eq:`vib:app:pendulum:ith`
we get a differential equation for the angle :math:`\theta(t)`. This latter
equation is ordered as

.. math::
   :label: vib:app:pendulum:thetaeq
        
        m\ddot\theta + + \frac{1}{2}C_D\varrho AL|\dot\theta|\dot\theta
        + \frac{mg}{L}\sin\theta = 0{\thinspace .}
        
        

Two initial conditions are needed: :math:`\theta=\Theta` and :math:`\dot\theta = \Omega`.
Normally, the pendulum motion is started from rest, which means :math:`\Omega =0`.

Equation :eq:`vib:app:pendulum:thetaeq` fits the general model
used in :eq:`vib:ode2` in the section :ref:`vib:model2` if we define
:math:`u=\theta`, :math:`f(u^{\prime}) = \frac{1}{2}C_D\varrho AL|\dot\theta|\dot\theta`,
:math:`s(u) = L^{-1}mg\sin u`, and :math:`F=0`.
If the body is a sphere with radius :math:`R`, we can take :math:`C_D=0.4` and :math:`A=\pi R^2`.

The motion of a compound or physical pendulum where the wire is a rod with
mass, can be modeled very similarly. The governing equation is
:math:`I\boldsymbol{a} = \boldsymbol{T}` where :math:`I` is the moment of inertia of the entire body about
the point :math:`(x_0,y_0)`, and :math:`\boldsymbol{T}` is the sum of moments of the forces
with respect to :math:`(x_0,y_0)`. The vector equation reads

.. math::
         \boldsymbol{r}\times(-S{\boldsymbol{i}_r} + mg(-\sin\theta\boldsymbol{i}_{\theta} + \cos\theta{\boldsymbol{i}_r}) -
        \frac{1}{2}C_D\varrho AL^2|\dot\theta|\dot\theta\boldsymbol{i}_{\theta})
        = I(L\ddot\theta\dot\theta\boldsymbol{i}_{\theta} - L\dot\theta^2{\boldsymbol{i}_r}){\thinspace .}

The component equation in :math:`\boldsymbol{i}_{\theta}` direction gives the equation of motion
for :math:`\theta(t)`:

.. math::
   :label: vib:app:pendulum:thetaeq_physical
        
        I\ddot\theta + \frac{1}{2}C_D\varrho AL^3|\dot\theta|\dot\theta
        + mgL\sin\theta = 0{\thinspace .}
        
        

[**hpl 11**: Scale the equations to arrive at the model problem with :math:`\sin\theta` spring.]

.. _vib:app:pendulum_elastic:

Motion of an elastic pendulum
-----------------------------

Consider a pendulum as in Figure :ref:`vib:app:pendulum:fig_problem`, but
this time the wire is elastic. The length of the wire when it is not
stretched is :math:`L_0`, while :math:`L(t)` is the stretched
length at time :math:`t` during the motion.

Stretching the elastic wire a distance :math:`\Delta L`
gives rise to a spring force :math:`k\Delta L` in the opposite direction of the
stretching. Let :math:`\boldsymbol{n}` be a unit normal vector along the wire
from the point :math:`\boldsymbol{r}_0=(x_0,y_0)` and in the direction of :math:`\boldsymbol{i}_{\theta}`,
see Figure :ref:`vib:app:pendulum:fig_forces` for definition of
:math:`(x_0,y_0)` and :math:`\boldsymbol{i}_{\theta}`. Obviously, we have :math:`\boldsymbol{n}=\boldsymbol{i}_{\theta}`, but in
this modeling of an elastic pendulum we do not need polar coordinates.
Instead, it is more straightforward to develop the equation in
Cartesian coordinates.

A mathematical expression for :math:`\boldsymbol{n}` is

.. math::
         \boldsymbol{n} = \frac{\boldsymbol{r}-\boldsymbol{r}_0}{L(t)},

where :math:`L(t)=||\boldsymbol{r}-\boldsymbol{r}_0||` is the current length of the elastic wire.
The position vector :math:`\boldsymbol{r}` in Cartesian coordinates reads
:math:`\boldsymbol{r}(t) = x(t)\boldsymbol{i} + y(t)\boldsymbol{j}`, where :math:`\boldsymbol{i}` and :math:`\boldsymbol{j}` are unit vectors
in the :math:`x` and :math:`y` directions, respectively.
It is convenient to introduce the Cartesian components :math:`n_x` and :math:`n_y`
of the normal vector:

.. math::
         \boldsymbol{n} = \frac{\boldsymbol{r}-\boldsymbol{r}_0}{L(t)} = \frac{x(t)-x_0}{L(t)}\boldsymbol{i} + \frac{y(t)-y_0}{L(t)}\boldsymbol{j} = n_x\boldsymbol{i} + n_y\boldsymbol{j}{\thinspace .}

The stretch :math:`\Delta L` in the wire is

.. math::
         \Delta t = L(t) - L_0{\thinspace .}

The force in the wire is then :math:`-S\boldsymbol{n}=-k\Delta L\boldsymbol{n}`.

The other forces are the gravity and the air resistance, just
as in Figure :ref:`vib:app:pendulum:fig_forces`. The main difference
is that we have a *model* for :math:`S` in terms of the motion (as soon as
we have expressed :math:`\Delta L` by :math:`\boldsymbol{r}`). For simplicity, we drop
the air resistance term (but :ref:`vib:exer:pendulum_elastic_drag`
asks you to include it).

Newton's second law of motion applied to the body now results in

.. math::
   :label: vib:app:pendulum_elastic:eq1
        
        m\ddot\boldsymbol{r} = -k(L-L_0)\boldsymbol{n} - mg\boldsymbol{j}
        
        

The two components of
:eq:`vib:app:pendulum_elastic:eq1` are

.. math::
   :label: _auto27
        
        \ddot x = -\frac{k}{m}(L-L_0)n_x,
        
        

.. math::
   :label: vib:app:pendulum_elastic:eq2a
          
         
        

.. math::
   :label: vib:app:pendulum_elastic:eq2b
          
        \ddot y = - \frac{k}{m}(L-L_0)n_y - g
        {\thinspace .}
        

Remarks about an elastic vs a non-elastic pendulum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the derivation of the ODEs for an elastic pendulum is more
straightforward than for a classical, non-elastic pendulum,
since we avoid the details
with polar coordinates, but instead work with Newton's second law
directly in Cartesian coordinates. The reason why we can do this is that
the elastic pendulum undergoes a general two-dimensional motion where
all the forces are known or expressed as functions of :math:`x(t)` and :math:`y(t)`,
such that we get two ordinary differential equations.
The motion of the non-elastic pendulum, on the other hand, is constrained:
the body has to move along a circular path, and the force :math:`S` in the
wire is unknown.

The non-elastic pendulum therefore leads to
a *differential-algebraic* equation, i.e., ODEs for :math:`x(t)` and :math:`y(t)`
combined with an extra constraint :math:`(x-x_0)^2 + (y-y_0)^2 = L^2`
ensuring that the motion takes place along a circular path.
The extra constraint (equation) is compensated by an extra unknown force
:math:`-S\boldsymbol{n}`. Differential-algebraic equations are normally hard
to solve, especially with pen and paper.
Fortunately, for the non-elastic pendulum we can do a
trick: in polar coordinates the unknown force :math:`S` appears only in the
radial component of Newton's second law, while the unknown
degree of freedom for describing the motion, the angle :math:`\theta(t)`,
is completely governed by the asimuthal component. This allows us to
decouple the unknowns :math:`S` and :math:`\theta`. But this is a kind of trick and
not a widely applicable method. With an elastic pendulum we use straightforward
reasoning with Newton's 2nd law and arrive at a standard ODE problem that
(after scaling) is easy solve on a computer.

Initial conditions
~~~~~~~~~~~~~~~~~~

What is the initial position of the body? We imagine that first the
pendulum hangs in equilibrium in its vertical position, and then it is
displaced an angle :math:`\Theta`. The equilibrium position is governed
by the ODEs with the accelerations set to zero.
The :math:`x` component leads to :math:`x(t)=x_0`, while the :math:`y` component gives

.. math::
         0 = - \frac{k}{m}(L-L_0)n_y - g = \frac{k}{m}(L(0)-L_0) - g\quad\Rightarrow\quad
        L(0) = L_0 + mg/k,

since :math:`n_y=-11` in this position. The corresponding :math:`y` value is then
from :math:`n_y=-1`:

.. math::
         y(t) = y_0 - L(0) = y_0 - (L_0 + mg/k){\thinspace .}

Let us now choose :math:`(x_0,y_0)` such that the body is at the origin
in the equilibrium position:

.. math::
         x_0 =0,\quad y_0 = L_0 + mg/k{\thinspace .}

Displacing the body an angle :math:`\Theta` to the right leads to the
initial position

.. math::
         x(0)=(L_0+mg/k)\sin\Theta,\quad y(0)=(L_0+mg/k)(1-\cos\Theta){\thinspace .}

The initial velocities can be set to zero: :math:`x'(0)=y'(0)=0`.

The complete ODE problem
~~~~~~~~~~~~~~~~~~~~~~~~

We can summarize all the equations as follows:

.. math::
        
        \ddot x &= -\frac{k}{m}(L-L_0)n_x,
        \\ 
        \ddot y &= -\frac{k}{m}(L-L_0)n_y - g,
        \\ 
        L &= \sqrt{(x-x_0)^2 + (y-y_0)^2},
        \\ 
        n_x &= \frac{x-x_0}{L},
        \\ 
        n_y &= \frac{y-y_0}{L},
        \\ 
        x(0) &= (L_0+mg/k)\sin\Theta,
        \\ 
        x'(0) &= 0,
        \\ 
        y(0) & =(L_0+mg/k)(1-\cos\Theta),
        \\ 
        y'(0) &= 0{\thinspace .}
        

We insert :math:`n_x` and :math:`n_y`  in the ODEs:

.. math::
   :label: vib:app:pendulum_elastic:x
        
        \ddot x = -\frac{k}{m}\left(1 -\frac{L_0}{L}\right)(x-x_0),
        
        

.. math::
   :label: vib:app:pendulum_elastic:y
          
        \ddot y = -\frac{k}{m}\left(1 -\frac{L_0}{L}\right)(y-y_0) - g,
        
        

.. math::
   :label: vib:app:pendulum_elastic:L
          
        L = \sqrt{(x-x_0)^2 + (y-y_0)^2},
        
        

.. math::
   :label: vib:app:pendulum_elastic:x0
          
        x(0) = (L_0+mg/k)\sin\Theta,
        
        

.. math::
   :label: vib:app:pendulum_elastic:vx0
          
        x'(0) = 0,
        
        

.. math::
   :label: vib:app:pendulum_elastic:y0
          
        y(0)  =(L_0+mg/k)(1-\cos\Theta),
        
        

.. math::
   :label: vib:app:pendulum_elastic:vy0
          
        y'(0) = 0{\thinspace .}
        
        

Scaling          (3)
~~~~~~~~~~~~~~~~~~~~

The elastic pendulum model can be used to study both an elastic pendulum
and a classic, non-elastic pendulum. The latter problem is obtained
by letting :math:`k\rightarrow\infty`. Unfortunately,
a serious problem with the ODEs
:eq:`vib:app:pendulum_elastic:x`-:eq:`vib:app:pendulum_elastic:y` is that for large :math:`k`, we have a very large factor :math:`k/m` multiplied by a
very small number :math:`1-L_0/L`, since for large :math:`k`, :math:`L\approx L_0` (very
small deformations of the wire). The product is subject to
significant round-off errors for many relevant physical values of
the parameters. To circumvent the problem, we introduce a scaling. This
will also remove physical parameters from the problem such that we end
up with only one dimensionless parameter,
closely related to the elasticity of the wire. Simulations can then be
done by setting just this dimensionless parameter.

The characteristic length can be taken such that in equilibrium, the
scaled length is unity, i.e., the characteristic length is :math:`L_0+mg/k`:

.. math::
         \bar x = \frac{x}{L_0+mg/k},\quad \bar y = \frac{y}{L_0+mg/k}{\thinspace .}

We must then also work with the scaled length :math:`\bar L = L/(L_0+mg/k)`.

Introducing :math:`\bar t=t/t_c`, where :math:`t_c` is a characteristic time we
have to decide upon later, one gets

.. math::
        
        \frac{d^2\bar x}{d\bar t^2} &=
        -t_c^2\frac{k}{m}\left(1 -\frac{L_0}{L_0+mg/k}\frac{1}{\bar L}\right)\bar x,\\ 
        \frac{d^2\bar y}{d\bar t^2} &=
        -t_c^2\frac{k}{m}\left(1 -\frac{L_0}{L_0+mg/k}\frac{1}{\bar L}\right)(\bar y-1)
        -t_c^2\frac{g}{L_0 + mg/k},\\ 
        \bar L &= \sqrt{\bar x^2 + (\bar y-1)^2},\\ 
        \bar x(0) &= \sin\Theta,\\ 
        \bar x'(0) &= 0,\\ 
        \bar y(0) & = 1 - \cos\Theta,\\ 
        \bar y'(0) &= 0{\thinspace .}
        

For a non-elastic pendulum with small angles, we know that the
frequency of the oscillations are :math:`\omega = \sqrt{L/g}`. It is therefore
natural to choose a similar expression here, either the length in
the equilibrium position,

.. math::
         t_c^2 = \frac{L_0+mg/k}{g}{\thinspace .}

or simply the unstretched length,

.. math::
         t_c^2 = \frac{L_0}{g}{\thinspace .}

These quantities are not very different (since the elastic model
is valid only for quite small elongations), so we take the latter as it is
the simplest one.

The ODEs become

.. math::
        
        \frac{d^2\bar x}{d\bar t^2} &=
        -\frac{L_0k}{mg}\left(1 -\frac{L_0}{L_0+mg/k}\frac{1}{\bar L}\right)\bar x,\\ 
        \frac{d^2\bar y}{d\bar t^2} &=
        -\frac{L_0k}{mg}\left(1 -\frac{L_0}{L_0+mg/k}\frac{1}{\bar L}\right)(\bar y-1)
        -\frac{L_0}{L_0 + mg/k},\\ 
        \bar L &= \sqrt{\bar x^2 + (\bar y-1)^2}{\thinspace .}
        

We can now identify a dimensionless number

.. math::
         \beta = \frac{L_0}{L_0 + mg/k} = \frac{1}{1+\frac{mg}{L_0k}},

which is the ratio of the unstretched length and the
stretched length in equilibrium. The non-elastic pendulum will have
:math:`\beta =1` (:math:`k\rightarrow\infty`).
With :math:`\beta` the ODEs read

.. math::
   :label: vib:app:pendulum_elastic:x:s
        
        \frac{d^2\bar x}{d\bar t^2} =
        -\frac{\beta}{1-\beta}\left(1- \frac{\beta}{\bar L}\right)\bar x,
        
        

.. math::
   :label: vib:app:pendulum_elastic:y:s
          
        \frac{d^2\bar y}{d\bar t^2} =
        -\frac{\beta}{1-\beta}\left(1- \frac{\beta}{\bar L}\right)(\bar y-1)
        -\beta,
        
        

.. math::
   :label: vib:app:pendulum_elastic:L:s
          
        \bar L = \sqrt{\bar x^2 + (\bar y-1)^2},
        
        

.. math::
   :label: vib:app:pendulum_elastic:x0:s
          
        \bar x(0) = (1+\epsilon)\sin\Theta,
        
        

.. math::
   :label: vib:app:pendulum_elastic:vx0:s
          
        \frac{d\bar x}{d\bar t}(0) = 0,
        
        

.. math::
   :label: vib:app:pendulum_elastic:y0:s
          
        \bar y(0) = 1 - (1+\epsilon)\cos\Theta,
        
        

.. math::
   :label: vib:app:pendulum_elastic:vy0:s
          
        \frac{d\bar y}{d\bar t}(0) = 0,
        
        

We have here added a parameter :math:`\epsilon`, which is an additional
downward stretch of the wire at :math:`t=0`. This parameter makes it possible
to do a desired test: vertical oscillations of the pendulum. Without
:math:`\epsilon`, starting the motion from :math:`(0,0)` with zero velocity will
result in :math:`x=y=0` for all times (also a good test!), but with
an initial stretch so the body's position is :math:`(0,\epsilon)`, we
will have oscillatory vertical motion with amplitude :math:`\epsilon` (see
:ref:`vib:exer:pendulum_elastic`).

Remark on the non-elastic limit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We immediately see that as :math:`k\rightarrow\infty` (i.e., we obtain a non-elastic
pendulum), :math:`\beta\rightarrow 1`, :math:`\bar L\rightarrow 1`, and we have
very small values :math:`1-\beta\bar L^{-1}` divided by very small values
:math:`1-\beta` in the ODEs. However, it turns out that we can set :math:`\beta`
very close to one and obtain a path of the body that within the visual
accuracy of a plot does not show any elastic oscillations.
(Should the division of very small values become a problem, one can
study the limit by L'Hospital's rule:

.. math::
         \lim_{\beta\rightarrow 1}\frac{1 - \beta \bar L^{-1}}{1-\beta}
        = \frac{1}{\bar L},

and use the limit :math:`\bar L^{-1}` in the ODEs for :math:`\beta` values very
close to 1.)

.. _vib:app:bouncing_ball:

Bouncing ball
-------------

A bouncing ball is a body in free vertically fall until it impacts the
ground. During the impact, some kinetic energy is lost, and a new
motion upwards with reduced velocity starts.  At some point the
velocity close to the ground is so small that the ball is considered
to be finally at rest.

The motion of the ball falling in air is governed by Newton's second
law :math:`F=ma`, where :math:`a` is the acceleration of the body, :math:`m` is the mass,
and :math:`F` is the sum of all forces. Here, we neglect the air resistance
so that gravity :math:`-mg` is the only force. The height of the ball is
denoted by :math:`h` and :math:`v` is the velocity. The relations between :math:`h`, :math:`v`, and
:math:`a`,

.. math::
         h'(t)= v(t),\quad v'(t) = a(t),

combined with Newton's second law gives the ODE model

.. math::
   :label: vib:app:bouncing:ball:h2eq
        
        h^{\prime\prime}(t) = -g,
        
        

or expressed alternatively as a system of first-order equations:

.. math::
   :label: vib:app:bouncing:ball:veq
        
        v'(t) = -g,
         
        

.. math::
   :label: vib:app:bouncing:ball:heq
          
        h'(t) = v(t){\thinspace .}
        
        

These equations govern the motion as long as the ball is away from
the ground by a small distance :math:`\epsilon_h > 0`. When :math:`h<\epsilon_h`,
we have two cases.

1. The ball impacts the ground, recognized by a sufficiently large negative
   velocity (:math:`v<-\epsilon_v`). The velocity then changes sign and is
   reduced by a factor :math:`C_R`, known as the `coefficient of restitution <http://en.wikipedia.org/wiki/Coefficient_of_restitution>`__.
   For plotting purposes, one may set :math:`h=0`.

2. The motion stops, recognized by a sufficiently small velocity
   (:math:`|v|<\epsilon_v`) close to the ground.

Electric circuits
-----------------

Although the term "mechanical vibrations" is used in the present
document, we must mention that the same type of equations arise
when modeling electric circuits.
The current :math:`I(t)` in a
circuit with an inductor with inductance :math:`L`, a capacitor with
capacitance :math:`C`, and overall resistance :math:`R`, is governed by

.. math::
   :label: _auto28
        
        \ddot I + \frac{R}{L}\dot I + \frac{1}{LC}I = \dot V(t),
        
        

where :math:`V(t)` is the voltage source powering the circuit.
This equation has the same form as the general model considered in
Section ref:ref:`vib:model2` if we set :math:`u=I`, :math:`f(u^{\prime}=bu^{\prime}`
and define :math:`b=R/L`, :math:`s(u) = L^{-1}C^{-1}u`, and :math:`F(t)=\dot V(t)`.

Exercises
=========

.. --- begin exercise ---

.. _vib:exer:resonance:

Exercise 18: Simulate resonance
-------------------------------

.. index:: resonance

We consider the scaled ODE model
:eq:`vib:app:mass_gen:scaled` from the section :ref:`vib:app:mass_gen`.
After scaling, the amplitude of :math:`u` will have a size about unity
as time grows and the effect of the initial conditions die out due
to damping. However, as :math:`\gamma\rightarrow 1`, the amplitude of :math:`u`
increases, especially if :math:`\beta` is small. This effect is called
*resonance*. The purpose of this exercise is to explore resonance.

**a)**
Figure out how the ``solver`` function in ``vib.py`` can be called
for the scaled ODE :eq:`vib:app:mass_gen:scaled`.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

**b)**
Run :math:`\gamma =5, 1.5, 1.1, 1` for :math:`\beta=0.005, 0.05, 0.2`.
For each :math:`\beta` value, present an image with plots of :math:`u(t)` for
the four :math:`\gamma` values.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

Filename: ``resonance``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:sliding_box:

Exercise 19: Simulate oscillations of a sliding box
---------------------------------------------------

Consider a sliding box on a flat surface as modeled in the section :ref:`vib:app:mass_sliding`. As spring force we choose the nonlinear
formula

.. math::
         s(u) = \frac{k}{\alpha}\tanh(\alpha u) = ku + \frac{1}{3}\alpha^2 ku^3 + \frac{2}{15}\alpha^4 k u^5 + {\mathcal{O}(u^6)}{\thinspace .}

**a)**
Plot :math:`g(u)=\alpha^{-1}\tanh(\alpha u)` for various values of :math:`\alpha`.
Assume :math:`u\in [-1,1]`.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

**b)**
Scale the equations using :math:`I` as scale for :math:`u` and :math:`m/k` as
time scale.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

**c)**
Implement the scaled model in b). Run it for some values of
the dimensionless parameters.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

Filename: ``sliding_box``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:bouncing:ball:

Exercise 20: Simulate a bouncing ball
-------------------------------------

The section :ref:`vib:app:bouncing_ball` presents a model for a bouncing
ball.
Choose one of the two ODE formulation, :eq:`vib:app:bouncing:ball:h2eq` or
:eq:`vib:app:bouncing:ball:veq`-:eq:`vib:app:bouncing:ball:heq`,
and simulate the motion of a bouncing ball. Plot :math:`h(t)`. Think about how to
plot :math:`v(t)`.

.. --- begin hint in exercise ---

**Hint.**
A naive implementation may get stuck in repeated impacts for large time
step sizes. To avoid this situation, one can introduce a state
variable that holds the mode of the motion: free fall, impact, or rest.
Two consecutive impacts imply that the motion has stopped.

.. --- end hint in exercise ---

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

Filename: ``bouncing_ball``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:pendulum_elastic:

Exercise 21: Simulate an elastic pendulum
-----------------------------------------

The section :ref:`vib:app:pendulum_elastic` describes a model for an elastic
pendulum, resulting in a system of two ODEs. The purpose of this
exercise is to implement the scaled model, test the software, and
generalize the model.

**a)**
Write a function ``simulate``
that can simulate an elastic pendulum using the scaled model.
The function should have the following arguments:

.. code-block:: python

        def simulate(
            beta=0.9,                 # dimensionless parameter
            Theta=30,                 # initial angle in degrees
            epsilon=0,                # initial stretch of wire
            num_periods=6,            # simulate for num_periods
            time_steps_per_period=60, # time step resolution
            plot=True,                # make plots or not
            ):

To set the total simulation time and the time step, we
use our knowledge of the scaled, classical, non-elastic pendulum:
:math:`u^{\prime\prime} + u = 0`, with solution
:math:`u = \Theta\cos \bar t`.
The period of these oscillations is :math:`P=2\pi`
and the frequency is unity. The time
for simulation is taken as ``num_periods`` times :math:`P`. The time step
is set as :math:`P` divided by ``time_steps_per_period``.

The ``simulate`` function should return the arrays of
:math:`x`, :math:`y`, :math:`\theta`, and :math:`t`, where :math:`\theta = \tan^{-1}(x/(1-y))` is
the angular displacement of the elastic pendulum corresponding to the
position :math:`(x,y)`.

If ``plot`` is ``True``, make a plot of :math:`\bar y(\bar t)`
versus :math:`\bar x(\bar t)`, i.e., the physical motion
of the mass at :math:`(\bar x,\bar y)`. Use the equal aspect ratio on the axis
such that we get a physically correct picture of the motion. Also
make a plot of :math:`\theta(\bar t)`, where :math:`\theta` is measured in degrees.
If :math:`\Theta < 10` degrees, add a plot that compares the solutions of
the scaled, classical, non-elastic pendulum and the elastic pendulum
(:math:`\theta(t)`).

Although the mathematics here employs a bar over scaled quantities, the
code should feature plain names ``x`` for :math:`\bar x`, ``y`` for :math:`\bar y`, and
``t`` for :math:`\bar t` (rather than ``x_bar``, etc.). These variable names make
the code easier to read and compare with the mathematics.

.. --- begin hint in exercise ---

**Hint 1.**
Equal aspect ratio is set by ``plt.gca().set_aspect('equal')`` in
Matplotlib (``import matplotlib.pyplot as plt``)
and by ``plot(..., daspect=[1,1,1], daspectmode='equal')``
in SciTools (``import scitools.std as plt``).

.. --- end hint in exercise ---

.. --- begin hint in exercise ---

**Hint 2.**
If you want to use Odespy to solve the equations, order the ODEs
like :math:`\dot \bar x, \bar x, \dot\bar y,\bar y` such that the Euler-Cromer
scheme can (also) be used (``odespy.EulerCromer``).

.. --- end hint in exercise ---

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

**b)**
Write a test function for testing that :math:`\Theta=0` and :math:`\epsilon=0`
gives :math:`x=y=0` for all times.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

**c)**
Write another test function for checking that the pure vertical
motion of the elastic pendulum is correct.
Start with simplifying the ODEs for pure vertical motion and show that
:math:`\bar y(\bar t)` fulfills a vibration equation with
frequency :math:`\sqrt{\beta/(1-\beta)}`. Set up the exact solution.

Write a test function that
uses this special case to verify the ``simulate`` function. There will
be numerical approximation errors present in the results from
``simulate`` so you have to believe in correct results and set a
(low) tolerance that corresponds to the computed maximum error.
Use a small :math:`\Delta t` to obtain a small numerical approximation error.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

**d)**
Make a function ``demo(beta, Theta)`` for simulating an elastic pendulum with a
given :math:`\beta` parameter and initial angle :math:`\Theta`. Use 600 time steps
per period to get every accurate results, and simulate for 3 periods.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

Filename: ``elastic_pendulum``.

.. --- end exercise ---

.. --- begin exercise ---

.. _vib:exer:pendulum_elastic_drag:

Exercise 22: Simulate an elastic pendulum with air resistance
-------------------------------------------------------------

This is a continuation :ref:`vib:exer:pendulum_elastic_drag`.
Air resistance on the body with mass :math:`m` can be modeled by the
force :math:`-\frac{1}{2}\varrho C_D A|\boldsymbol{v}|\boldsymbol{v}`,
where :math:`C_D` is a drag coefficient (0.2 for a sphere), :math:`\varrho`
is the density of air (1.2 :math:`\hbox{kg }\,{\hbox{m}}^{-3}`), :math:`A` is the
cross section area (:math:`A=\pi R^2` for a sphere, where :math:`R` is the radius),
and :math:`\boldsymbol{v}` is the velocity of the body.
Include air resistance in the original model, scale the model,
write a function ``simulate_drag`` that is a copy of the ``simulate``
function from :ref:`vib:exer:pendulum_elastic_drag`, but with the
new ODEs included, and show plots of how air resistance
influences the motion.

.. removed !bsol ... !esol environment (because of the command-line option --without_solutions)

Filename: ``elastic_pendulum_drag``.

.. Closing remarks for this Exercise

Remarks
~~~~~~~

Test functions are challenging to construct for the problem with
air resistance. You can reuse the tests from
:ref:`vib:exer:pendulum_elastic_drag` for ``simulate_drag``,
but these tests does not verify the new terms arising from air
resistance.

.. in vb_odespy examples: add 20 RK4 1000 to show RK4 in the long run

.. mu'' + bu' + ku = F(t)

.. set up analytical solution for reference

.. compare for F = sin qt, demonstrate resonance by having

.. F = sin q t and q = piecewise constant in time with four

.. levels: 0.1, 0.75 1, 1.25, 2 of the resonance frequency,

.. make each platou act for a while to see the effect

.. mu'' + bu' + s(u) = F(t) as exercise, pendulum

.. mu'' + f(x) + s() = F(t) via odespy RK4

.. odespy: ForwardBackward on a 2n system? Need special formula for first

.. step to ensure dt^2 accuracy there.

.. apps: planet around a star, box horizontal and vertical, bumpy,

.. jumping washing machine, pendulum, moored ship, look to Irgens

.. --- end exercise ---

References
==========

.. [Ref1]
   **H. P. Langtangen**. *Finite Difference Computing with Exponential Decay Models*,
   2015,
   `http://tinyurl.com/nclmcng/web <http://tinyurl.com/nclmcng/web>`_.

