import os
os.add_dll_directory("C:/Program Files/mingw64/bin")
os.add_dll_directory("C:/Users/Niurouxing/Desktop/NiurxWHL/external/openblas/bin")

import mimo as m
import multiprocessing
import time
from tqdm import tqdm

 


TxAntNum = 64
RxAntNum = 128
ModType = 8


startSNR=10
endSNR=50
stepSNR=1
SNR=range(startSNR,endSNR,stepSNR)

target = 20000000
samplesPreIter=10


def worker(queue,samples,snr):
    while True:
        errorBits,errorFrames = m.idd(TxAntNum,RxAntNum,ModType,snr,samples)
        queue.put((errorBits,errorFrames))
 
if __name__ == '__main__':
    BERall = []
    FERall = []
 

    for snr in SNR:
        num_processes = 4
        queue = multiprocessing.Queue()
        processes = [multiprocessing.Process(target=worker, args=(queue,samplesPreIter,snr)) for i in range(num_processes)]
        for process in processes:
            process.start()

        errorBits = 0
        errorFrames=0
        samples=0

        with tqdm(total=target) as pbar:
            while errorFrames < target:
                (newErrorBits,newErrorFrames) = queue.get()
        
                errorBits += newErrorBits
                errorFrames += newErrorFrames
                samples += samplesPreIter

                pbar.set_postfix(snr=snr,误码率=errorBits/samples/TxAntNum/ModType,误帧率=errorFrames/samples,总错比特=str(errorBits) ,总错帧=str(errorFrames),总样本=str(samples))
                
                pbar.update(newErrorFrames)
         
        for process in processes:
            process.terminate()

        BERall.append(errorBits/samples/TxAntNum/ModType)
        FERall.append(errorFrames/samples)