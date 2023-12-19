import mimo as m

errorBits,errorFrames=m.detection(4,8,4,30,10000)

print("errorBits: ",errorBits)
print("errorFrames: ",errorFrames)