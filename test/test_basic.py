import mimo as m
 
from multiprocessing import Pool

TxAntNum = 4
RxAntNum = 6
ModType = 4
SNRdB = 12
sample = 1000
beta=0.9
loop=9
NSAIter=18

alphaVec = [0.5] * (2 * TxAntNum)
accuaVec = [1.0] * (NSAIter)
 
errorBits,errorFrames = m.EPAwNSADet(TxAntNum, RxAntNum, ModType, SNRdB, sample,beta,NSAIter,loop,alphaVec,accuaVec)

print ("Error Bits: ", errorBits)
print ("Error Frames: ", errorFrames)

# just to wrap the function into a worker function
def work(a,b):
    errorBits,errorFrames = m.EPAwNSADet(TxAntNum, RxAntNum, ModType, SNRdB, sample,beta,NSAIter,loop,a,b)
    return errorBits,errorFrames


# multi-processing test


p = Pool(processes=10)

results = p.map(work, [(alphaVec,accuaVec)]*10)

p.close()
p.join()

for errorBits,errorFrames in results:
    print ("Error Bits: ", errorBits)
    print ("Error Frames: ", errorFrames)


