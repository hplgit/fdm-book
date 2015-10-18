def rub_solver1(N, version='scalar'):
    I = lambda x: sin(2*x*pi/L)

    class Action:
        def __init__(self):
            self.solutions = []
            self.time_level_counter = 0

        def __call__(self, u, x, t, n):
            if self.time_level_counter % N == 0:
                self.solutions.append(u.copy())
            self.time_level_counter += 1

    action = Action()
    n = 100; T = 6; L = 10
    dt, x, cpu = solver(I, None, None, 1.0, 0, 0,
                        L, n, 1, T,
                        user_action=action, version=version)
    print 'CPU time:', cpu
    print 'Max value in final u:', arrmax(action.solutions[-1])

