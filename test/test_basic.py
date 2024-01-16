import mimo as m
 
 
TxAntNum = 32
RxAntNum = 64
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



                                       