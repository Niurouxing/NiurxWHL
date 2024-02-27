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
from pymoo.termination import get_termination

hyperparameters = {
    "ModType": 8,
    "SNRdB": 20 ,
    "TxAntNum": 64,
    "RxAntNum": 128,
    "samplesPreIter": 20000,
    "errorBitsTarget": 10000,
    "pop_size": 32,
    "max_gen": 1000,
    "beta": 0.9,
    "loop": 8,
    "NSAIter": 10
    ,
}


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
            n_var=2 * hyperparameters["NSAIter"],
            n_obj=1,
            n_constr=0,
            xl=-100,
            xu=100,
        )
        # 创建一个类成员进程池
        self.pool = Pool(processes=os.cpu_count())

    def _evaluate(self, x, out, *args, **kwargs):
        pop_size = x.shape[0]
        f = np.zeros((pop_size, 1))

        alphaVec = x[:, 0 : hyperparameters["NSAIter"]]
        accuaVec = x[:, hyperparameters["NSAIter"] :]

        # 使用已创建的进程池
        results = self.pool.starmap(work, zip(alphaVec, accuaVec))

        for i in range(pop_size):
            errorBits, errorFrames = results[i]
            f[i] = (
                errorBits
                / hyperparameters["samplesPreIter"]
                / hyperparameters["ModType"]
                / hyperparameters["TxAntNum"]
            )

        out["F"] = f

    def __del__(self):
        # 在对象被销毁时关闭进程池
        self.pool.close()
        self.pool.join()


class MyCallback(Callback):
    def __init__(self) -> None:
        super().__init__()
        self.log = 0

    def notify(self, algorithm):
        self.log += 1
        if self.log % 10 == 0:
            currentArg = algorithm.pop.get("X")
            currentPer = algorithm.pop.get("F")

            minIndex = np.argmin(currentPer)
            minErrorBits = currentPer[minIndex]
            print("minErrorBits:", minErrorBits)

            filename = f'pop_{hyperparameters["ModType"]}_{hyperparameters["TxAntNum"]}_{hyperparameters["RxAntNum"]}_{hyperparameters["SNRdB"]}_{hyperparameters["beta"]}_{hyperparameters["loop"]}_{hyperparameters["NSAIter"]}.txt'

            # 将数组转换为逗号分隔的字符串保存
            with open(filename, "w") as f:
                for individual in currentArg:
                    line = ",".join(map(str, individual))
                    f.write(f"{line}\n")

        # if the first time, just report the best performance
        if self.log == 1:
 
            currentPer = algorithm.pop.get("F")

            minIndex = np.argmin(currentPer)
            minErrorBits = currentPer[minIndex]
       
            print("minErrorBits:", minErrorBits)
 

def main():
    filename = f'pop_{hyperparameters["ModType"]}_{hyperparameters["TxAntNum"]}_{hyperparameters["RxAntNum"]}_{hyperparameters["SNRdB"]}_{hyperparameters["beta"]}_{hyperparameters["loop"]}_{hyperparameters["NSAIter"]}.txt'

    if os.path.exists(filename):
        print(
            f"Found existing population file: {filename}. Loading as initial population."
        )
        with open(filename, "r") as f:
            lines = f.readlines()
            all_pop = np.array(
                [list(map(float, line.strip().split(","))) for line in lines]
            )
    else:
        print("No existing population file found. Generating a new initial population.")
        pop_size = hyperparameters["pop_size"]

        good_init1 = [0.5] * (hyperparameters["NSAIter"])
        good_init2 = [0] * (hyperparameters["NSAIter"])
        good_init = good_init1 + good_init2
        all_pop = np.tile(good_init, (pop_size, 1))

        # add some random noise to the initial population, +-0.1 for each element
        noise = np.random.uniform(
            -0.3,
            0.3,
            (pop_size, 2 *hyperparameters["NSAIter"]),
        )
        all_pop = all_pop + noise

    myPro = myProblem()
    print("The problem has been initialized!")
    pop_size = hyperparameters["pop_size"]
    n_gen = hyperparameters["max_gen"]

    termination = get_termination("n_gen", hyperparameters["max_gen"])
    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=all_pop,
        # sampling=IntegerRandomSampling(),
        termination=termination,
        crossover=SBX(prob=0.9, eta=15, vtype=float),
        mutation=PM(
            eta=15,
            vtype=float,
        ),
        eliminate_duplicates=True,
    )

    callback = MyCallback()
 
    res = minimize(myPro, algorithm, callback=callback, seed=1, verbose=True)
    print("normal exit")

if __name__ == "__main__":
    while True:
        main()

