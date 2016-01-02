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
    print random
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
    current_pos = x0
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

def random_walks1D(x0, N, p, num_walks=1, num_times=1,
                   random=random):
    """Simulate num_walks random walks from x0 with N steps."""
    # random is argument so we can use np.random instead
    # and use it for testing equivalence with random_walk1D_vec

    # Store position in step k in position[k]
    position = np.zeros(N+1)
    position[0] = x0*(N+1)
    position2 = np.zeros(N+1)
    position2[0] = x0**2*(N+1)
    # Histogram at num_times selected time points
    pos_hist = []
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]
    print 'save at times', pos_hist_times

    for n in range(num_walks):
        current_pos = x0
        for k in range(N):
            r = random.uniform(0, 1)
            if r <= p:
                current_pos -= 1
            else:
                current_pos += 1
            position[k+1] += current_pos
            if k in pos_hist_times:
                pos_hist.append(position[k])
            position2[k+1] += current_pos**2
            #print k+1, current_pos, freq[x2i(position[k+1], N)]
    pos_hist = np.array(pos_hist).reshape(num_walks, num_times)
    return position, position2, pos_hist, np.array(pos_hist_times)

def random_walks1D_vec(x0, N, p, num_walks=1, num_times=1):
    """Vectorized version of random_walks1D."""
    # Store position in step k in position[k]
    position = np.zeros(N+1)
    position[0] = x0*num_walks
    position2 = np.zeros(N+1)
    position2[0] = x0**2*num_walks
    # Histogram at num_times selected time points
    pos_hist = np.zeros((num_walks, num_times))
    pos_hist_times = [(N//num_times)*i for i in range(num_times)]

    for n in range(num_walks):
        current_pos = x0
        r = np.random.uniform(0, 1, size=N)
        steps = np.where(r <= p, -1, 1)
        walk = x0 + np.cumsum(steps)  # Positions of this walk
        position[1:] += walk
        position2[1:] += walk**2
        pos_hist[n,:] = position[pos_hist_times]
    return position, position2, pos_hist, np.array(pos_hist_times)

def test_random_walks1D():
    # For fixed seed, check that scalar and vectorized versions
    # produce the same result
    x0 = 0;  N = 4;  p = 0.5

    # First check that random_walks1D for 1 walk reproduces
    # the walk in random_walk1D
    num_walks = 1
    np.random.seed(10)
    computed = random_walks1D(
        x0, N, p, num_walks, random=np.random)
    np.random.seed(10)
    expected = random_walk1D(
        x0, N, p, random=np.random)
    assert (computed[0] == expected).all()
    print 'positions 1 walk serial:', computed[0]

    # Same for vectorized version
    np.random.seed(10)
    computed = random_walks1D_vec(x0, N, p, num_walks)
    np.random.seed(10)
    expected = random_walk1D_vec(x0, N, p)
    assert (computed[0] == expected).all()

    print 'positions 1 walk vec:', computed[0]
    # Test multiple walks
    num_walks = 3
    num_times = N
    np.random.seed(10)
    serial_computed = random_walks1D(
        x0, N, p, num_walks, num_times, random=np.random)
    np.random.seed(10)
    vectorized_computed = random_walks1D_vec(
        x0, N, p, num_walks, num_times)
    # positions: [0`, 1, 0, 1, 2]
    # Can test without tolerance since everything is +/- 1
    return_values = ['pos', 'pos2', 'pos_hist', 'pos_hist_times']
    for s, v, r in zip(serial_computed,
                       vectorized_computed,
                       return_values):
        msg = '%s: %s (serial) vs %s (vectorized)' % (r, s, v)
        assert (s == v).all(), msg

def demo_random_walks1D(N=1000, num_walks=10000):
    import time
    pos, pos2, freq, hist, hist_times = random_walks1D_vec(
        x0=0, N=N, p=0.5, num_walks=num_walks)
    t1 = time.clock()
    print 'random walk: %.1fs' % (t1-t0)
    E_X = pos/float(num_walks)
    Var_X = pos2/float(num_walks) - E_X**2
    t2 = time.clock()
    if N <= 50:
        print pos
    print 'E/Var: %.1fs' % (t2-t1)


    plt.figure()
    """
    plt.plot(E_X)
    plt.title('Expected position')
    """
    plt.plot(E_X)
    plt.title('Average position over the walks')
    plt.savefig('tmp1.png');  plt.savefig('tmp1.png')
    t3 = time.clock()
    print 'plotting 1: %.1fs' % (t3-t2)
    plt.figure()
    """
    plt.plot(Var_X)
    plt.title('Variance of position')
    """
    plt.plot(Var_X)
    plt.title('Average position squared over the walks')
    plt.savefig('tmp2.png');  plt.savefig('tmp2.png')
    t4 = time.clock()
    print 'plotting 2: %.1fs' % (t4-t3)
    plt.figure()
    #plt.bar(x, freq)
    # WRONG: strip away zeros in freq or track min/max pos ever
    # debug that freq is correct, unit test for that? two walks are enough,
    # with N=4 :-)[[[
    #nonzeros = freq[-10,:] > 0
    #plt.bar(range(len(freq[-10,nonzeros])), freq[-10,nonzeros])
    plt.hist(hist[:,-3], bins=50) #, normed=True)
    plt.savefig('tmp3.png');  plt.savefig('tmp3.png')
    t5 = time.clock()
    print 'plotting 3: %.1fs' % (t5-t4)
    plt.show()

if __name__ == '__main__':
    demo_fig_random_walk1D(N=200)
    #demo_random_walks1D()
    #demo_random_walk1D()
