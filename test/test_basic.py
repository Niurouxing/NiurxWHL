import mimo as m



TxAntNum = 4
RxAntNum = 8
ModType = 4
SNR=15
samples=10000


errorBits,errorFrames = m.det(TxAntNum,RxAntNum,ModType,SNR,samples)

print("误码率=",errorBits)
print("误帧率=",errorFrames)

