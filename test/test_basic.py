import mimo as m
 
 
TxAntNum = 4
RxAntNum = 8
ModType = 6
SNRdB = 10
sample = 100

for i in range(10):

    errorBits,errorFrames = m.idd(TxAntNum, RxAntNum, ModType, SNRdB, sample)

    print ("Error Bits: ", errorBits)
    print ("Error Frames: ", errorFrames)



                                       