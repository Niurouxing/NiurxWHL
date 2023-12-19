import mimo as m

errorBits,errorFrames=m.detection(4,8,4,10,100)

print("errorBits: ",errorBits)
print("errorFrames: ",errorFrames)