"""Conversion examples python side."""
import numpy as np

import mimo as mimo

TxAntNum = 4
RxAntNum = 4
SNRdB = 20
ModType = 4
sample = 10000

err = mimo.Detection(TxAntNum, RxAntNum, SNRdB, ModType, sample)

print(err)


a = 1
