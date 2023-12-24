import mimo as m



TxAntNum = 8
RxAntNum = 8
ModType = 4
SNR=10
samples=10000


m.detectionInit(TxAntNum,RxAntNum,ModType,SNR,True)

m.ExBsPInit(3,4)

for i in range(samples):
    m.execute()
    newErrorBits,newErrorFrames = m.report()

    print(newErrorBits,newErrorFrames)





