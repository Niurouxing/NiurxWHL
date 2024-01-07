import os
os.add_dll_directory("C:/Program Files/mingw64/bin")
os.add_dll_directory("C:/Users/Niurouxing/Desktop/NiurxWHL/external/openblas/bin")
import mimo as m
 
 
TxAntNum = 64
RxAntNum = 128
ModType = 8
SNRdB = 10
sample = 100

for i in range(10):

    errorBits,errorFrames = m.idd(TxAntNum, RxAntNum, ModType, SNRdB, sample)

    print ("Error Bits: ", errorBits)
    print ("Error Frames: ", errorFrames)



                                       