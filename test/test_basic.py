import mimo as m
import multiprocessing
import time

TxAntNum = 32
RxAntNum = 64
SNR = 10
ModType = 4
sample = 10000



def run_detection():
    errorBits, errorFrames = m.detection(TxAntNum, RxAntNum, ModType, SNR, sample)

# Create a list to store the processes
processes = []

start = time.time()

# Create and start a process for each iteration
for _ in range(multiprocessing.cpu_count()):
    process = multiprocessing.Process(target=run_detection)
    process.start()
    processes.append(process)

# Wait for all processes to finish
for process in processes:
    process.join()

end = time.time()
print("Time taken: ", end - start)
 
 
