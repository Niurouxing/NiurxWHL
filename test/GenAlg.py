import numpy as np
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.core.callback import Callback
from pymoo.util.display.column import Column
from pymoo.util.display.output import Output
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
import mimo as m
from pymoo.visualization.scatter import Scatter
import multiprocessing
from multiprocessing import Pool
import time
import os


hyperparameters = {
    "ModType": 8,
    "SNRdB": 20,
    "TxAntNum": 64,
    "RxAntNum": 128,
    "samplesPreIter": 1000,
    "errorBitsTarget": 10000,
    "pop_size": 10,
    "max_gen": 1000,
    "beta": 0.9,
    "loop": 7,
    "NSAIter": 40,
}


def scale_coordinates(matrix):
    """
    Scales the coordinates in the given matrix uniformly to have the same average values
    for the horizontal and vertical coordinates.
    """
    # Compute the means of the horizontal and vertical coordinates
    x_mean = np.mean(matrix[:, 0])
    y_mean = np.mean(matrix[:, 1])

    # Compute the scaling factors for the horizontal and vertical coordinates
    x_scale = y_mean / x_mean
    y_scale = 1.0

    # Scale the coordinates in the matrix
    scaled_matrix = np.copy(matrix)
    scaled_matrix[:, 0] *= x_scale
    scaled_matrix[:, 1] *= y_scale

    return scaled_matrix


def closest_to_origin(matrix):
    """
    Finds the point closest to the origin among all the scaled coordinates in the given matrix.
    """
    # Scale the coordinates uniformly
    scaled_matrix = scale_coordinates(matrix)

    # Compute the distances from the origin to each point in the scaled matrix
    distances = np.sqrt(np.sum(np.square(scaled_matrix), axis=1))

    # Find the index of the point with the smallest distance to the origin
    min_index = np.argmin(distances)

    # Return the coordinates of the closest point
    return min_index


def work(alphaVec, accuaVec):
 
    errorBits, errorFrames = m.EPAwNSADet(
        hyperparameters["TxAntNum"],
        hyperparameters["RxAntNum"],
        hyperparameters["ModType"],
        hyperparameters["SNRdB"],
        hyperparameters["samplesPreIter"],
        hyperparameters["beta"],
        hyperparameters["NSAIter"],
        hyperparameters["loop"],
        alphaVec,
        accuaVec,
    )
    return errorBits, errorFrames


class myProblem(Problem):
    def __init__(self):
        super().__init__(
            n_var=2 * hyperparameters["TxAntNum"] + hyperparameters["NSAIter"],
            n_obj=1,
            n_constr=0,
        )

    def _evaluate(self, x, out, *args, **kwargs):
        global hyperparameters

        pop_size = x.shape[0]
        f = np.zeros((pop_size, 1))

        alphaVec = x[:, 0 : 2 * hyperparameters["TxAntNum"]]
        accuaVec = x[:, 2 * hyperparameters["TxAntNum"] :]

        p = Pool(processes=10)
        res_l = list()

        for i in range(pop_size):
            res = p.apply_async(work, args=(alphaVec[i], accuaVec[i]))
            res_l.append(res)

        p.close()
        p.join()

        for i in range(pop_size):
            errorBits, errorFrames = res_l[i].get()
            f[i] = errorBits 

        out["F"] = f


class MyCallback(Callback):
    def __init__(self) -> None:
        super().__init__()
        self.current = np.zeros((1, 4))

    def notify(self, algorithm):
        currentArg = algorithm.pop.get("X")
        currentPer = algorithm.pop.get("F")

        # find smallest errorBits in current population
        minIndex = np.argmin(currentPer)
        minErrorBits = currentPer[minIndex]
        minArg = currentArg[minIndex]
        print("minErrorBits:", minErrorBits)
        print("minArg:", minArg)
        print("currentPer:", currentPer)


if __name__ == "__main__":
    myPro = myProblem()
    print("The problem has been initialized!")
    pop_size = hyperparameters["pop_size"]
    n_gen = hyperparameters["max_gen"]

    good_init1 = [0.5] * (2 * hyperparameters["TxAntNum"])
    good_init2 = [1.0] * (hyperparameters["NSAIter"])
    good_init = good_init1 + good_init2
    # print("good_init:",good_init)
    all_pop = np.tile(good_init, (pop_size, 1))

    # add some random noise to the initial population, +-0.1 for each element
    noise = np.random.uniform(
        -0.1,
        0.1,
        (pop_size, 2 * hyperparameters["TxAntNum"] + hyperparameters["NSAIter"]),
    )
    all_pop = all_pop + noise

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=all_pop,
        # sampling=IntegerRandomSampling(),
        crossover=SBX(prob=0.9, eta=15, vtype=float),
        mutation=PM(
            eta=20,
            vtype=float,
        ),
        eliminate_duplicates=True,
    )

    callback = MyCallback()
    res = minimize(
        myPro, algorithm, callback=callback, seed=1, verbose=True
    )
