"""Conversion examples python side."""
import numpy as np

import mimo as mimo

TxAntNum = 4
RxAntNum = 8
SNRdB = 10.0
ModType = 4
sample = 100000

for i in range(100):

    err = mimo.Detection(TxAntNum, RxAntNum, SNRdB, ModType, sample)

    print(err)

    a=1

print(err)


a = 1
