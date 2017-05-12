# -*- coding: utf-8 -*-
"""
Class implementation for solving of the wave equation 
u_tt = (c**2*u_x)_x + f(x,t) with t in [0,T] and x in (0,L). 
We have u=U_0 or du/dn=0 on x=0, and u=u_L or du/dn=0 on x = L.  
For simplicity, we use a constant c here and compare with a
known exact solution.
"""
import time, glob, shutil, os
import numpy as np

class Parameters(object):
    def __init__(self):
        """
        Subclasses must initialize self.prm with
        parameters and default values, self.type with
        the corresponding types, and self.help with
        the corresponding descriptions of parameters.
        self.type and self.help are optional, but
        self.prms must be complete and contain all parameters.
        """
        pass

    def ok(self):
        """Check if attr. prm, type, and help are defined."""
        if hasattr(self, 'prm') and \
           isinstance(self.prm, dict) and \
           hasattr(self, 'type') and \
           isinstance(self.type, dict) and \
           hasattr(self, 'help') and \
           isinstance(self.help, dict):
            return True
        else:
            raise ValueError(
                'The constructor in class %s does not '\
                'initialize the\ndictionaries '\
                'self.prm, self.type, self.help!' %
                self.__class__.__name__)

    def _illegal_parameter(self, name):
        """Raise exception about illegal parameter name."""
        raise ValueError(
            'parameter "%s" is not registered.\nLegal '\
            'parameters are\n%s' %
            (name, ' '.join(list(self.prm.keys()))))

    def set(self, **parameters):
        """Set one or more parameters."""
        for name in parameters:
            if name in self.prm:
                self.prm[name] = parameters[name]
            else:
                self._illegal_parameter(name)

    def get(self, name):
        """Get one or more parameter values."""
        if isinstance(name, (list,tuple)):   # get many?
            for n in name:
                if n not in self.prm:
                    self._illegal_parameter(name)
            return [self.prm[n] for n in name]
        else:
            if name not in self.prm:
                self._illegal_parameter(name)
            return self.prm[name]

    def __getitem__(self, name):
        """Allow obj[name] indexing to look up a parameter."""
        return self.get(name)

    def __setitem__(self, name, value):
        """
        Allow obj[name] = value syntax to assign a parameter's value.
        """
        return self.set(name=value)

    def define_command_line_options(self, parser=None):
        self.ok()
        if parser is None:
            import argparse
            parser = argparse.ArgumentParser()

        for name in self.prm:
            tp = self.type[name] if name in self.type else str
            help = self.help[name] if name in self.help else None
            parser.add_argument(
                '--' + name, default=self.get(name), metavar=name,
                type=tp, help=help)

        return parser

    def init_from_command_line(self, args):
        for name in self.prm:
            self.prm[name] = getattr(args, name)


class Problem(Parameters):
    """
    Physical parameters for the wave equation 
    u_tt = (c**2*u_x)_x + f(x,t) with t in [0,T] and
    x in (0,L). The problem definition is implied by
    the method of manufactured solution, choosing 
    u(x,t)=x(L-x)(1+t/2) as our solution. This solution
    should be exactly reproduced when c is const.
    """
            
    def __init__(self):
        self.prm  = dict(L=2.5, c=1.5, T=18)
        self.type = dict(L=float, c=float, T=float)
        self.help = dict(L='1D domain',
                         c='coefficient (wave velocity) in PDE',
                         T='end time of simulation')
    def u_exact(self, x, t):
        L = self['L']
        return x*(L-x)*(1+0.5*t)   
    def I(self, x):
        return self.u_exact(x, 0)
    def V(self, x):
        return 0.5*self.u_exact(x, 0)
    def f(self, x, t):
        c = self['c']
        return 2*(1+0.5*t)*c**2
    def U_0(self, t):
        return self.u_exact(0, t)
    U_L = None
        
        
class Solver(Parameters):
    """
    Numerical parameters for solving the wave equation 
    u_tt = (c**2*u_x)_x + f(x,t) with t in [0,T] and
    x in (0,L). The problem definition is implied by
    the method of manufactured solution, choosing 
    u(x,t)=x(L-x)(1+t/2) as our solution. This solution
    should be exactly reproduced, provided c is const.
    We simulate in [0, L/2] and apply a symmetry condition
    at the end x=L/2.    
    """ 
    
    def __init__(self, problem):       
        self.problem = problem
        self.prm  = dict(C = 0.75, Nx=3, stability_safety_factor=1.0)
        self.type = dict(C=float, Nx=int, stability_safety_factor=float)
        self.help = dict(C='Courant number',
                         Nx='No of spatial mesh points',
                         stability_safety_factor='stability factor')
                         
        from UniformFDMesh import Mesh, Function
        # introduce some local help variables to ease reading
        L_end = self.problem['L']
        dx = (L_end/2)/float(self['Nx'])
        t_interval = self.problem['T']
        dt = dx*self['stability_safety_factor']*self['C']/ \
                                    float(self.problem['c'])
        self.m = Mesh(L=[0,L_end/2], 
                      d=[dx], 
                      Nt = int(round(t_interval/float(dt))),
                      T=t_interval)
        # The mesh function f will, after solving, contain
        # the solution for the whole domain and all time steps.
        self.f = Function(self.m, num_comp=1, space_only=False)

    def solve(self, user_action=None, version='scalar'):
        # ...use local variables to ease reading
        L, c, T = self.problem['L c T'.split()]
        L = L/2     # compute with half the domain only (symmetry)
        C, Nx, stability_safety_factor = self[
                            'C Nx stability_safety_factor'.split()]
        dx = self.m.d[0]
        I = self.problem.I
        V = self.problem.V
        f = self.problem.f
        U_0 = self.problem.U_0
        U_L = self.problem.U_L
        Nt = self.m.Nt
        t = np.linspace(0, T, Nt+1)      # Mesh points in time
        x = np.linspace(0, L, Nx+1)      # Mesh points in space
        
        # Make sure dx and dt are compatible with x and t
        dx = x[1] - x[0]
        dt = t[1] - t[0]
    
        # Treat c(x) as array
        if isinstance(c, (float,int)):
            c = np.zeros(x.shape) + c
        elif callable(c):
            # Call c(x) and fill array c
            c_ = np.zeros(x.shape)
            for i in range(Nx+1):
                c_[i] = c(x[i])
            c = c_
    
        q = c**2
        C2 = (dt/dx)**2; dt2 = dt*dt    # Help variables in the scheme
       
        # Wrap user-given f, I, V, U_0, U_L if None or 0
        if f is None or f == 0:
            f = (lambda x, t: 0) if version == 'scalar' else \
                lambda x, t: np.zeros(x.shape)
        if I is None or I == 0:
            I = (lambda x: 0) if version == 'scalar' else \
                lambda x: np.zeros(x.shape)
        if V is None or V == 0:
            V = (lambda x: 0) if version == 'scalar' else \
                lambda x: np.zeros(x.shape)
        if U_0 is not None:
            if isinstance(U_0, (float,int)) and U_0 == 0:
                U_0 = lambda t: 0
        if U_L is not None:
            if isinstance(U_L, (float,int)) and U_L == 0:
                U_L = lambda t: 0
    
        # Make hash of all input data
        import hashlib, inspect
        data = inspect.getsource(I) + '_' + inspect.getsource(V) + \
               '_' + inspect.getsource(f) + '_' + str(c) + '_' + \
               ('None' if U_0 is None else inspect.getsource(U_0)) + \
               ('None' if U_L is None else inspect.getsource(U_L)) + \
               '_' + str(L) + str(dt) + '_' + str(C) + '_' + str(T) + \
               '_' + str(stability_safety_factor)
        hashed_input = hashlib.sha1(data).hexdigest()
        if os.path.isfile('.' + hashed_input + '_archive.npz'):
            # Simulation is already run
            return -1, hashed_input
                
        # use local variables to make code closer to mathematical
        # notation in computational scheme
        u_1 = self.f.u[0,:]
        u   = self.f.u[1,:]
        
        import time;  t0 = time.clock()  # CPU time measurement
    
        Ix = range(0, Nx+1)
        It = range(0, Nt+1)
    
        # Load initial condition into u_1
        for i in range(0,Nx+1):
            u_1[i] = I(x[i])
    
        if user_action is not None:
            user_action(u_1, x, t, 0)
        
        # Special formula for the first step
        for i in Ix[1:-1]:
            u[i] = u_1[i] + dt*V(x[i]) + \
            0.5*C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i]) - \
                    0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
            0.5*dt2*f(x[i], t[0])
    
        i = Ix[0]
        if U_0 is None:
            # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
            # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
            ip1 = i+1
            im1 = ip1  # i-1 -> i+1
            u[i] = u_1[i] + dt*V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_0(dt)
    
        i = Ix[-1]
        if U_L is None:
            im1 = i-1
            ip1 = im1  # i+1 -> i-1
            u[i] = u_1[i] + dt*V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_L(dt)
    
        if user_action is not None:
            user_action(u, x, t, 1)
            
        for n in It[1:-1]:
            # u corresponds to u^{n+1} in the mathematical scheme            
            u_2 = self.f.u[n-1,:]
            u_1 = self.f.u[n,:]
            u   = self.f.u[n+1,:]
            
            # Update all inner points
            if version == 'scalar':
                for i in Ix[1:-1]:
                    u[i] = - u_2[i] + 2*u_1[i] + \
                        C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i])  - \
                            0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
                    dt2*f(x[i], t[n])
    
            elif version == 'vectorized':
                u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
                C2*(0.5*(q[1:-1] + q[2:])*(u_1[2:] - u_1[1:-1]) -
                    0.5*(q[1:-1] + q[:-2])*(u_1[1:-1] - u_1[:-2])) + \
                dt2*f(x[1:-1], t[n])
            else:
                raise ValueError('version=%s' % version)
    
            # Insert boundary conditions
            i = Ix[0]
            if U_0 is None:
                # Set boundary values
                # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
                # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
                ip1 = i+1
                im1 = ip1
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
                dt2*f(x[i], t[n])
            else:
                u[i] = U_0(t[n+1])
    
            i = Ix[-1]
            if U_L is None:
                im1 = i-1
                ip1 = im1
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
                dt2*f(x[i], t[n])
            else:
                u[i] = U_L(t[n+1])
    
            if user_action is not None:
                if user_action(u, x, t, n+1):
                    break
                
        cpu_time = time.clock() - t0
        return cpu_time, hashed_input
        
    def assert_no_error(self):
        """Run through mesh and check error"""   
        Nx = self['Nx']
        Nt = self.m.Nt
        L, T = self.problem['L T'.split()]
        L = L/2     # only half the domain used (symmetry)
        x = np.linspace(0, L, Nx+1)     # Mesh points in space        
        t = np.linspace(0, T, Nt+1)     # Mesh points in time
        
        for n in range(len(t)):
            u_e = self.problem.u_exact(x, t[n])
            diff = np.abs(self.f.u[n,:] - u_e).max()
            print 'diff:', diff
            tol = 1E-13
            assert diff < tol
                
def test_quadratic_with_classes():
    """
    Check the scalar and vectorized versions for a quadratic 
    u(x,t)=x(L-x)(1+t/2) that is exactly reproduced,
    provided c(x) is constant. We simulate in [0, L/2] and 
    apply a symmetry condition at the end x=L/2.
    """

    problem = Problem()
    solver = Solver(problem)

    # Read input from the command line
    parser = problem.define_command_line_options()
    parser = solver. define_command_line_options(parser)
    args = parser.parse_args()
    problem.init_from_command_line(args)
    solver. init_from_command_line(args)

    print parser.parse_args()              # parameters ok?

    solver.solve()   
    print 'Check error.........................'
    solver.assert_no_error()


if __name__ == '__main__':
    test_quadratic_with_classes()
