import random, numpy as np
import matplotlib.pyplot as plt
random.seed(10)
np.random.seed(10)

def random_walk1D(x0, N, p, random=random):
    """1D random walk with 1 particle and N moves."""
    # random is argument so we can use np.random instead
    # and use it for testing equivalence with random_walk1D_vec

    # Store position in step k in position[k]
    position = np.zeros(N+1)
    position[0] = x0
    current_pos = x0
    for k in range(N):
        r = random.uniform(0, 1)
        if r <= p:
            current_pos -= 1
        else:
            current_pos += 1
        position[k+1] = current_pos
    return position

def random_walk1D_vec(x0, N, p):
    """Vectorized version of random_walk1D."""
    # Store position in step k in position[k]
    position = np.zeros(N+1)
    position[0] = x0
    r = np.random.uniform(0, 1, size=N)
    steps = np.where(r <= p, -1, 1)
    position[1:] = x0 + np.cumsum(steps)
    return position

def test_random_walk1D():
    # For fixed seed, check that scalar and vectorized versions
    # produce the same result
    x0 = 2;  N = 4;  p = 0.6
    np.random.seed(10)
    scalar_computed = random_walk1D(x0, N, p, random=np.random)
    np.random.seed(10)
    vectorized_computed = random_walk1D_vec(x0, N, p)
    assert (scalar_computed == vectorized_computed).all()
    # We do not need tolerance since all arithmetics is w/int
    #diff = np.abs(scalar_computed - vectorized_computed).max()
    # Positions: [2, 3, 2, 3, 4]
    #tol = 1E-14
    #assert diff < 1E-14

def demo_random_walk1D(N=50000):
    np.random.seed(10)
    pos = random_walk1D_vec(x0=0, N=N, p=0.5)
    plt.figure()
    plt.plot(pos)
    plt.savefig('tmp1.pdf');  plt.savefig('tmp1.png')
    plt.figure()
    plt.plot(pos*pos)
    plt.savefig('tmp2.pdf');  plt.savefig('tmp2.png')
    plt.show()

def demo_fig_random_walk1D(N=200):
    """Make ensamble of positions (to illustrate E[] operator)."""
    np.random.seed(10)
    num_plots = 4
    for n in range(num_plots):
        plt.subplot(num_plots, 1, n+1)
        pos = random_walk1D_vec(x0=0, N=N, p=0.5)
        plt.plot(pos)
        plt.axis([0, N, -15, 20])
    plt.savefig('tmp.pdf');  plt.savefig('tmp.png')
    plt.show()

def demo_random_walk1D_timing():
    import time
    x0 = 0
    N = 10000000
    p = 0.5

    t0 = time.clock()
    np.random.seed(10)
    pos = random_walk1D(x0, N, p, random=np.random)
    t1 = time.clock()
    cpu_scalar = t1 - t0
    print 'CPU scalar: %.1f' % cpu_scalar
    np.random.seed(10)
    pos = random_walk1D_vec(x0, N, p)
    t2 = time.clock()
    cpu_vec = t2 - t1
    print 'CPU vectorized: %.1f' % cpu_vec
    print 'CPU scalar/vectorized: %.1f' % (cpu_scalar/cpu_vec)

def random_walks1D(x0, N, p, num_walks=1, num_times=1,
                   random=random):
    """Simulate num_walks random walks from x0 with N steps."""
    position = np.zeros(N+1)    # Accumulated positions
    position[0] = x0*num_walks
    position2 = np.zeros(N+1)   # Accumulated positions**2
    position2[0] = x0**2*num_walks
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]
    #print 'save hist:', post_hist_times

    for n in range(num_walks):
        num_times_counter = 0
        current_pos = x0
        for k in range(N):
            if k in pos_hist_times:
                #print 'save, k:', k, num_times_counter, n
                pos_hist[n,num_times_counter] = current_pos
                num_times_counter += 1
            # current_pos corresponds to step k+1
            r = random.uniform(0, 1)
            if r <= p:
                current_pos -= 1
            else:
                current_pos += 1
            position [k+1] += current_pos
            position2[k+1] += current_pos**2
    return position, position2, pos_hist, np.array(pos_hist_times)

def random_walks1D_vec1(x0, N, p, num_walks=1, num_times=1):
    """Vectorized version of random_walks1D."""
    position  = np.zeros(N+1)    # Accumulated positions
    position2 = np.zeros(N+1)    # Accumulated positions**2
    walk = np.zeros(N+1)         # Positions of current walk
    walk[0] = x0
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]

    for n in range(num_walks):
        r = np.random.uniform(0, 1, size=N)
        steps = np.where(r <= p, -1, 1)
        walk[1:] = x0 + np.cumsum(steps)  # Positions of this walk
        position  += walk
        position2 += walk**2
        pos_hist[n,:] = walk[pos_hist_times]
    return position, position2, pos_hist, np.array(pos_hist_times)

def random_walks1D_vec2(x0, N, p, num_walks=1, num_times=1):
    """Vectorized version of random_walks1D; no loops."""
    position  = np.zeros(N+1)    # Accumulated positions
    position2 = np.zeros(N+1)    # Accumulated positions**2
    walks = np.zeros((num_walks, N+1))  # Positions of each walk
    walks[:,0] = x0
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]

    r = np.random.uniform(0, 1, size=N*num_walks)
    steps = np.where(r <= p, -1, 1).reshape(num_walks, N)
    walks[:,1:] = x0 + np.cumsum(steps, axis=1)
    position  = np.sum(walks,    axis=0)
    position2 = np.sum(walks**2, axis=0)
    pos_hist[:,:] = walks[:,pos_hist_times]
    return position, position2, pos_hist, np.array(pos_hist_times)

def test_random_walks1D():
    # For fixed seed, check that scalar and vectorized versions
    # produce the same result
    x0 = 0;  N = 4;  p = 0.5

    # First, check that random_walks1D for 1 walk reproduces
    # the walk in random_walk1D
    num_walks = 1
    np.random.seed(10)
    computed = random_walks1D(
        x0, N, p, num_walks, random=np.random)
    np.random.seed(10)
    expected = random_walk1D(
        x0, N, p, random=np.random)
    assert (computed[0] == expected).all()

    # Same for vectorized versions
    np.random.seed(10)
    computed = random_walks1D_vec1(x0, N, p, num_walks)
    np.random.seed(10)
    expected = random_walk1D_vec(x0, N, p)
    assert (computed[0] == expected).all()
    np.random.seed(10)
    computed = random_walks1D_vec2(x0, N, p, num_walks)
    np.random.seed(10)
    expected = random_walk1D_vec(x0, N, p)
    assert (computed[0] == expected).all()

    # Second, check multiple walks: scalar == vectorized
    num_walks = 3
    num_times = N
    np.random.seed(10)
    serial_computed = random_walks1D(
        x0, N, p, num_walks, num_times, random=np.random)
    np.random.seed(10)
    vectorized1_computed = random_walks1D_vec1(
        x0, N, p, num_walks, num_times)
    np.random.seed(10)
    vectorized2_computed = random_walks1D_vec2(
        x0, N, p, num_walks, num_times)
    # positions: [0, 1, 0, 1, 2]
    # Can test without tolerance since everything is +/- 1
    return_values = ['pos', 'pos2', 'pos_hist', 'pos_hist_times']
    for s, v, r in zip(serial_computed,
                       vectorized1_computed,
                       return_values):
        msg = '%s: %s (serial) vs %s (vectorized)' % (r, s, v)
        assert (s == v).all(), msg
    for s, v, r in zip(serial_computed,
                       vectorized2_computed,
                       return_values):
        msg = '%s: %s (serial) vs %s (vectorized)' % (r, s, v)
        assert (s == v).all(), msg

def demo_random_walks1D(N=1000, num_walks=10000, EX_minmax=None):
    import time
    t0 = time.clock()
    pos, pos2, hist, hist_times = random_walks1D_vec1(
        x0=0, N=N, p=0.5, num_walks=num_walks, num_times=10,)
    t1 = time.clock()
    print 'histogram times:', hist_times
    print 'random walk: %.1fs' % (t1-t0)
    E_X = pos/float(num_walks)
    Var_X = pos2/float(num_walks) - E_X**2
    if N <= 50:
        print pos

    plt.figure()
    plt.plot(E_X)
    if EX_minmax is not None:
        plt.axis([0, N, EX_minmax[0], EX_minmax[1]])
    plt.title('Expected position (%d walks)' % num_walks)
    plt.savefig('tmp1.png');  plt.savefig('tmp1.pdf')
    plt.figure()
    plt.plot(Var_X)
    plt.title('Variance of position (%d walks)' % num_walks)
    plt.savefig('tmp2.png');  plt.savefig('tmp2.pdf')

    # Compare histogram and diffusion equation solution
    plt.figure()
    a = 0.5
    exact = lambda x, t: 1./np.sqrt(4*np.pi*t*a)*np.exp(-x**2/(4.0*t*a))
    hist_time_index = -2
    n, bins, patches = plt.hist(hist[:,hist_time_index], bins=30, normed=True)
    x = np.linspace(bins[0], bins[-1], 301)
    t = hist_times[hist_time_index]
    plt.plot(x, exact(x, t), 'r--')
    plt.title('Histogram of positions (%d walks)' % num_walks)
    plt.savefig('tmp3.png');  plt.savefig('tmp3.pdf')
    plt.show()

def demo_fig_random_walks1D():
    """Make figures with statistics and dependence on no of walks."""
    import shutil, os
    N = 1000
    num_walks = [100, 10000, 100000, 1000000]
    for n in num_walks:
        np.random.seed(10)  # Use same seq. for all experiments
        # Tailor y axis for E[X]
        if n == 100:
            demo_random_walks1D(N=N, num_walks=n, EX_minmax=None)
        else:
            demo_random_walks1D(N=N, num_walks=n, EX_minmax=[-0.1, 0.5])
        d = 'tmp_%d' % n
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.mkdir(d)
        for p in 1,2,3:
            os.rename('tmp%d.png' % p, os.path.join(d, 'tmp%d.png' % p))
            os.rename('tmp%d.pdf' % p, os.path.join(d, 'tmp%d.pdf' % p))
    # E[X] tmp1 figs
    plots = ['EX', 'VarX', 'HistX']
    for j, plot in enumerate(plots):
        for ext in 'png', 'pdf':
            files = [os.path.join('tmp_%d' % n, 'tmp%d.%s' % (j+1, ext))
                     for n in num_walks]
            cmd = 'doconce combine_images -%d %s rw1D_%s_%s.%s' % \
                  (3 if len(num_walks) == 3 else 2,
                   ' '.join(files), plot,
                   '_'.join([str(n) for n in num_walks]), ext)
            print cmd
            os.system(cmd)

def demo_random_walks1D_timing():
    import time
    x0 = 0
    N = 1000
    num_walks = 50000
    p = 0.5

    t0 = time.clock()
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walks1D(
        x0, N, p, num_walks, num_times=4,
        random=np.random)
    t1 = time.clock()
    cpu_scalar = t1 - t0
    print 'CPU scalar: %.1f' % cpu_scalar
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walks1D_vec1(
        x0, N, p, num_walks, num_times=4)
    t2 = time.clock()
    cpu_vec1 = t2 - t1
    print 'CPU vectorized1: %.1f' % cpu_vec1
    print 'CPU scalar/vectorized1: %.1f' % (cpu_scalar/cpu_vec1)
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walks1D_vec2(
        x0, N, p, num_walks, num_times=4)
    t3 = time.clock()
    cpu_vec2 = t3 - t2
    print 'CPU vectorized2: %.1f' % cpu_vec2
    print 'CPU scalar/vectorized2: %.1f' % (cpu_scalar/cpu_vec2)

def random_walk2D(x0, N, p, random=random):
    """2D random walk with 1 particle and N moves: N, E, W, S."""
    # Store position in step k in position[k]
    d = len(x0)
    position = np.zeros((N+1, d))
    position[0,:] = x0
    current_pos = np.array(x0, dtype=float)
    for k in range(N):
        r = random.uniform(0, 1)
        if r <= 0.25:
            current_pos += np.array([0, 1])   # Move north
        elif 0.25 < r <= 0.5:
            current_pos += np.array([1, 0])   # Move east
        elif 0.5 < r <= 0.75:
            current_pos += np.array([0, -1])  # Move south
        else:
            current_pos += np.array([-1, 0])  # Move west
        position[k+1,:] = current_pos
    return position

def demo_random_walk2D():
    x0 = (0,0)
    N = 200
    p = 0.5
    np.random.seed(10)
    pos = random_walk2D(x0, N, p, random=np.random)
    #print pos, pos.shape
    plt.plot(pos[:,0], pos[:,1])
    plt.savefig('tmp1.png');  plt.savefig('tmp1.pdf')
    plt.show()

def random_walkdD(x0, N, p, random=random):
    """Any-D (diagonal) random walk with 1 particle and N moves."""
    # Store position in step k in position[k]
    d = len(x0)
    position = np.zeros((N+1, d))
    position[0,:] = x0
    current_pos = np.array(x0, dtype=float)
    for k in range(N):
        for i in range(d):
            r = random.uniform(0, 1)
            if r <= p:
                current_pos[i] -= 1
            else:
                current_pos[i] += 1
        position[k+1,:] = current_pos
    return position

def random_walkdD_vec(x0, N, p):
    """Vectorized version of random_walkdD."""
    d = len(x0)
    # Store position in step k in position[k]
    position = np.zeros((N+1,d))
    position[0] = np.array(x0, dtype=float)
    r = np.random.uniform(0, 1, size=N*d)
    steps = np.where(r <= p, -1, 1).reshape(N,d)
    position[1:,:] = x0 + np.cumsum(steps, axis=0)
    return position

def demo_random_walkdD():
    x0 = (0,0)
    N = 200
    p = 0.5
    np.random.seed(10)
    pos = random_walkdD(x0, N, p, random=np.random)
    #print pos, pos.shape
    plt.plot(pos[:,0], pos[:,1])
    plt.savefig('tmp1.png');  plt.savefig('tmp1.pdf')
    plt.show()

def demo_random_walkdD_timing():
    import time
    x0 = (0,0)
    N = 4000000
    p = 0.5

    t0 = time.clock()
    np.random.seed(10)
    pos = random_walkdD(x0, N, p, random=np.random)
    t1 = time.clock()
    cpu_scalar = t1 - t0
    print 'CPU scalar: %.1f' % cpu_scalar
    np.random.seed(10)
    pos = random_walkdD_vec(x0, N, p)
    t2 = time.clock()
    cpu_vec = t2 - t1
    print 'CPU vectorized: %.1f' % cpu_vec
    print 'CPU scalar/vectorized: %.1f' % (cpu_scalar/cpu_vec)

def demo_fig_random_walkdD():
    x0 = (0,0)
    N = 5000
    p = 0.5
    n = 2   # nxn subplots
    f, axarr = plt.subplots(n, n, sharex=True, sharey=True)
    for i in range(n):
        for j in range(n):
            seed = 3*i+8*j
            np.random.seed(seed)
            pos = random_walkdD(x0, N, p, random=np.random)
            axarr[i,j].plot(pos[:,0], pos[:,1])
            #axarr[i,j].set_axis([-200,100,-200,50])
            #axarr[i,j].set_title('Seed=%d' % seed)
    plt.savefig('tmp1.png');  plt.savefig('tmp1.pdf')
    plt.show()

def test_ramdom_walkdD():
    x0 = (0,0)
    N = 7
    p = 0.5
    np.random.seed(10)
    scalar_computed = random_walkdD(x0, N, p, random=np.random)
    np.random.seed(10)
    vectorized_computed = random_walkdD_vec(x0, N, p)
    assert (scalar_computed == vectorized_computed).all()

def random_walksdD(x0, N, p, num_walks=1, num_times=1,
                   random=random):
    """Simulate num_walks random walks from x0 with N steps."""
    d = len(x0)
    position  = np.zeros((N+1, d))   # Accumulated positions
    position2 = np.zeros((N+1, d))   # Accumulated positions**2
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times, d))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]

    for n in range(num_walks):
        num_times_counter = 0
        current_pos = np.array(x0, dtype=float)
        for k in range(N):
            if k in pos_hist_times:
                pos_hist[n,num_times_counter,:] = current_pos
                num_times_counter += 1
            # current_pos corresponds to step k+1
            for i in range(d):
                r = random.uniform(0, 1)
                if r <= p:
                    current_pos[i] -= 1
                else:
                    current_pos[i] += 1
            position [k+1,:] += current_pos
            position2[k+1,:] += current_pos**2
    return position, position2, pos_hist, np.array(pos_hist_times)

def random_walksdD_vec(x0, N, p, num_walks=1, num_times=1):
    """Vectorized version of random_walks1D; no loops."""
    d = len(x0)
    position  = np.zeros((N+1, d))  # Accumulated positions
    position2 = np.zeros((N+1, d))  # Accumulated positions**2
    walks = np.zeros((num_walks, N+1, d))  # Positions of each walk
    walks[:,0,:] = x0
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times, d))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]

    r = np.random.uniform(0, 1, size=N*num_walks*d)
    steps = np.where(r <= p, -1, 1).reshape(num_walks, N, d)
    walks[:,1:,:] = x0 + np.cumsum(steps, axis=1)
    position  = np.sum(walks,    axis=0)
    position2 = np.sum(walks**2, axis=0)
    pos_hist[:,:,:] = walks[:,pos_hist_times,:]
    return position, position2, pos_hist, np.array(pos_hist_times)

def test_random_walksdD():
    # For fixed seed, check that scalar and vectorized versions
    # produce the same result
    x0 = (0,0);  N = 4;  p = 0.5

    # First, check that random_walksdD for 1 walk reproduces
    # the walk in random_walkdD
    num_walks = 1
    np.random.seed(10)
    computed = random_walksdD(
        x0, N, p, num_walks, random=np.random)
    np.random.seed(10)
    expected = random_walkdD(
        x0, N, p, random=np.random)
    assert (computed[0] == expected).all()

    # Same for vectorized version
    np.random.seed(10)
    computed = random_walksdD_vec(x0, N, p, num_walks)
    np.random.seed(10)
    expected = random_walkdD_vec(x0, N, p)
    assert (computed[0] == expected).all()

    # Second, check multiple walks: scalar == vectorized
    num_walks = 3
    num_times = N
    np.random.seed(10)
    serial_computed = random_walksdD(
        x0, N, p, num_walks, num_times, random=np.random)
    np.random.seed(10)
    vectorized_computed = random_walksdD_vec(
        x0, N, p, num_walks, num_times)
    return_values = ['pos', 'pos2', 'pos_hist', 'pos_hist_times']
    for s, v, r in zip(serial_computed,
                       vectorized_computed,
                       return_values):
        msg = '%s: %s\n%s (serial)\nvs\n%s\n%s (vectorized)' % \
              (r, s.shape, s, v.shape, v)
        assert (s == v).all(), msg

def demo_random_walksdD():
    x0 = (0,0)
    N = 1000
    num_walks = 1000
    p = 0.5
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walksdD(
        x0, N, p, num_walks, num_times=4,
        random=np.random)
    print pos_hist_times
    plt.figure()
    plt.plot(pos[:,0], pos[:,1])
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walksdD_vec(
        x0, N, p, num_walks, num_times=4)
    plt.figure()
    plt.plot(pos[:,0], pos[:,1])

    #plt.savefig('tmp1.png');  plt.savefig('tmp1.pdf')
    plt.show()

def demo_random_walksdD_timing():
    import time
    # Scalar version is almost independent of d, while vec is sensitive to d
    x0 = (0,0,0)
    #x0 = (0,)
    #x0 = (0,0)
    N = 1000
    num_walks = 10000
    p = 0.5

    t0 = time.clock()
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walksdD(
        x0, N, p, num_walks, num_times=4,
        random=np.random)
    t1 = time.clock()
    cpu_scalar = t1 - t0
    print 'CPU scalar: %.1f' % cpu_scalar
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walksdD_vec(
        x0, N, p, num_walks, num_times=4)
    t2 = time.clock()
    cpu_vec = t2 - t1
    print 'CPU vectorized: %.1f' % cpu_vec
    print 'CPU scalar/vectorized: %.1f' % (cpu_scalar/cpu_vec)

# Reversed loops in 1D walks (all functions named walks1D2)

def random_walks1D2(x0, N, p, num_walks=1, num_times=1,
                    random=random):
    """Simulate num_walks random walks from x0 with N steps."""
    position = np.zeros(N+1)    # Accumulated positions
    position[0] = x0*num_walks
    position2 = np.zeros(N+1)   # Accumulated positions**2
    position2[0] = x0**2*num_walks
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]

    current_pos = x0 + np.zeros(num_walks)
    num_times_counter = -1

    for k in range(N):
        if k in pos_hist_times:
	    num_times_counter += 1
	    store_hist = True  # Store histogram data for this k
	else:
	    store_hist = False

        for n in range(num_walks):
            # current_pos corresponds to step k+1
            r = random.uniform(0, 1)
            if r <= p:
                current_pos[n] -= 1
            else:
                current_pos[n] += 1
            position [k+1] += current_pos[n]
            position2[k+1] += current_pos[n]**2
            if store_hist:
                pos_hist[n,num_times_counter] = current_pos[n]
    return position, position2, pos_hist, np.array(pos_hist_times)

def random_walks1D2_vec1(x0, N, p, num_walks=1, num_times=1):
    """Vectorized version of random_walks1D2."""
    position  = np.zeros(N+1)    # Accumulated positions
    position2 = np.zeros(N+1)    # Accumulated positions**2
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]

    current_pos = np.zeros(num_walks)
    current_pos[0] = x0
    num_times_counter = -1

    for k in range(N):
        if k in pos_hist_times:
	    num_times_counter += 1
	    store_hist = True  # Store histogram data for this k
	else:
	    store_hist = False

        # Move all walks one step
        r = np.random.uniform(0, 1, size=num_walks)
        steps = np.where(r <= p, -1, 1)
        current_pos += steps
        position[k+1]  = np.sum(current_pos)
        position2[k+1] = np.sum(current_pos**2)
        if store_hist:
            pos_hist[:,num_times_counter] = current_pos
    return position, position2, pos_hist, np.array(pos_hist_times)

def test_random_walks1D2():
    x0 = 0;  N = 4;  p = 0.5
    num_walks = 3
    num_times = N
    np.random.seed(10)
    serial_computed = random_walks1D2(
        x0, N, p, num_walks, num_times, random=np.random)
    np.random.seed(10)
    vectorized_computed = random_walks1D2_vec1(
        x0, N, p, num_walks, num_times)
    # Can test without tolerance since everything is +/- 1
    return_values = ['pos', 'pos2', 'pos_hist', 'pos_hist_times']
    for s, v, r in zip(serial_computed,
                       vectorized_computed,
                       return_values):
        msg = '%s: %s (serial) vs %s (vectorized)' % (r, s, v)
        assert (s == v).all(), msg

def demo_random_walks1D2_timing():
    """Timing of random 1D walks with reversed loops."""
    import time
    x0 = 0
    N = 1000
    num_walks = 50000
    p = 0.5

    t0 = time.clock()
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walks1D2(
        x0, N, p, num_walks, num_times=4,
        random=np.random)
    t1 = time.clock()
    cpu_scalar = t1 - t0
    print 'CPU scalar: %.1f' % cpu_scalar
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walks1D2_vec1(
        x0, N, p, num_walks, num_times=4)
    t2 = time.clock()
    cpu_vec1 = t2 - t1
    print 'CPU vectorized1: %.1f' % cpu_vec1
    print 'CPU scalar/vectorized1: %.1f' % (cpu_scalar/cpu_vec1)
    # Compare with the other version without loops too
    np.random.seed(10)
    pos, pos2, pos_hist, pos_hist_times = random_walks1D_vec2(
        x0, N, p, num_walks, num_times=4)
    t3 = time.clock()
    cpu_vec2 = t3 - t2
    print 'CPU vectorized2: %.1f' % cpu_vec2
    print 'CPU scalar/vectorized2: %.1f' % (cpu_scalar/cpu_vec2)

if __name__ == '__main__':
    #demo_fig_random_walk1D(N=200)
    #demo_random_walks1D(N=1000, num_walks=1000000)
    #demo_fig_random_walks1D()
    #demo_random_walk1D()
    #demo_fig_random_walkdD()
    #demo_random_walkdD()
    #test_ramdom_walkdD()
    #test_random_walks1D()
    #demo_random_walksdD_timing()
    #demo_random_walks1D_timing()
    demo_random_walks1D2_timing()
    print '----'
    demo_random_walks1D_timing()
