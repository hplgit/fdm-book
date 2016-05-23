import numpy as np

def solver(I, V, g, dt, T):
    """
    Solve v' = - g(u,v), u'=v for t in (0,T], u(0)=I and v(0)=V,
    by the PEFRL method.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    v = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    xi = 0.1720865590295143
    lambda_ = -0.09156203075515678
    chi = -0.1616217622107222

    v[0] = V
    u[0] = I
    # Compare with eq 22 in http://arxiv.org/pdf/cond-mat/0110585.pdf
    for n in range(0, Nt):
        u_ = u[n] + xi*dt*v[n]
        v_ = v[n] + 0.5*(1-2*lambda_)*dt*g(u_, v[n])
        u_ = u_ + chi*dt*v_
        v_ = v+ + lambda_*dt*g(u_, v_)
        u_ = u_ + 0.5*(1-2*lambda_)*dt*v_
        v_ = v_ + lambda*g(u_, v_)
        u_ = u_ + chi*dt*v_
        v[n+1] = v_ + 0.5*(1-2*lambda_)*dt*g(u_, v_)
        u[n+1] = u_ + xi*dt*v[n+1]
        print 'v[%d]=%g, u[%d]=%g' % (n+1,v[n+1],n+1,u[n+1])
    return u, v, t
