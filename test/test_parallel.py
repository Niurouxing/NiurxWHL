import os
import multiprocessing
from tqdm import tqdm
import numpy as np
import mimo as m

TxAntNum = 64
RxAntNum = 128
ModType = 8

startSNR = 10
endSNR = 25
stepSNR = 3
SNR = np.arange(startSNR, endSNR, stepSNR)

target = 2000

def worker(queue, snr, stop_signal):
    m.detectionInit(TxAntNum, RxAntNum, ModType, snr, False)
    m.MMSEInit(False)
    while not stop_signal.is_set():
        m.execute()
        errorBits, errorFrames = m.report()
        queue.put((errorBits, errorFrames))

def run_experiment(snr):
    num_processes = os.cpu_count()
    queue = multiprocessing.Queue()
    stop_signal = multiprocessing.Event()

    processes = [multiprocessing.Process(target=worker, args=(queue, snr, stop_signal)) for _ in range(num_processes)]
    for process in processes:
        process.start()

    errorBits = 0
    errorFrames = 0
    samples = 0

    try:
        with tqdm(total=target) as pbar:
            while errorFrames < target:
                newErrorBits, newErrorFrames = queue.get()
                errorBits += newErrorBits
                errorFrames += newErrorFrames
                samples += 1
                pbar.set_postfix(snr=snr, 误码率=errorBits/samples/TxAntNum/ModType, 误帧率=errorFrames/samples, 总错比特=str(errorBits), 总错帧=str(errorFrames), 总样本=str(samples))
                pbar.update(newErrorFrames)
    finally:
        stop_signal.set()
        for process in processes:
            process.join()

    return errorBits/samples/TxAntNum/ModType, errorFrames/samples

if __name__ == '__main__':
    BERall = []
    FERall = []

    for snr in SNR:
        ber, fer = run_experiment(snr)
        BERall.append(ber)
        FERall.append(fer)