import numpy as np
from scipy import optimize
from multiprocessing import Pool


def square(x, a=1):
    return [np.sum(x**2 + a), 2*x]


def minimize(args):
    print(args)
    f, x, a = args
    res = optimize.minimize(f, x, method='BFGS', jac=True, args=[a])
    return res.x


if __name__ == "__main__":
    # your a values
    a = np.arange(1, 17)

    # initial guess for all the x values
    x = np.empty(len(a))
    x[:] = 25

    args = [(square, a[i], x[i]) for i in range(len(a))]
    p = Pool(4)
    a = np.array(p.map(minimize, args))
    a.reshape((4, 4))
    print(a)
