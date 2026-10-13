import struct, numpy as np
def load(fn):
    b=open(fn,'rb').read(); z=b.index(b'\0'); p=z+1+4+8+8+8+2
    rc,=struct.unpack('>i',b[p:p+4]); p+=4
    return rc, np.frombuffer(b[p:],dtype='>f4').astype(float).reshape(180,180)
def binof(x): return int((x+180)//2) % 180
def center(i): return -179+2*i
