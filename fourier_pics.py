from z3 import *

def getPath (o,dx,dy,xy):
    """the specified path taken by a long elephant starting at o"""
    # o is a complex number
    assert (dx in [1,-1])
    assert (dy in [1j,-1j])
    yield o
    kx = o+(dx if xy else dy)
    for i in range(4):
        yield kx
        kx+=dx+dy
def crange(c):
    for x in range(int(c.real)):
        for y in range(int(c.imag)):
            yield x+y*1j

l = []#list of all relevant paths
from collections import defaultdict
d=defaultdict(int) # the kernel we care about

sz=8
box = sz+sz*1j
#ker = [i+j*1j for i in range(sz) for j in range(1)]
#ker = [0j,1j,1+0j,1+1j]
ker=[0j]

def ky(c):
    return (c.real,c.imag)
def norm(p):
    """take a complex number to a canonical version (by applying rotations and reflections within the box"""
    off = (box-1-1j)/2
    q = p-(box-1-1j)/2 # center it
    return min([c*r+off for c in [q,q.conjugate()] for r in [1,-1,1j,-1j] ],key=ky)
    


counts = defaultdict(int)
dmul={} # the first 
for origin in crange(box):
    d[origin]=1
    counts[norm(origin)]+=1
    for dx in [1,-1]:
        for dy in [1j,-1j]:
            for xy in [True,False]:
                l.append(list(getPath(origin,dx,dy,xy)))

#the first successful kernel I found 
g1 = """  2  4  6  8  8  6  4  2
  4  8 11 16 16 11  8  4
  6 11 18 24 24 18 11  6
  8 16 24 30 30 24 16  8
  8 16 24 30 30 24 16  8
  6 11 18 24 24 18 11  6
  4  8 11 16 16 11  8  4
  2  4  6  8  8  6  4  2"""

g2 = """  2  1  1  1  1  1  1  2
  1  3  2  2  2  2  3  1
  1  2  4  4  4  4  2  1
  1  2  4  7  7  4  2  1
  1  2  4  7  7  4  2  1
  1  2  4  4  4  4  2  1
  1  3  2  2  2  2  3  1
  2  1  1  1  1  1  1  2"""


def getd(s):
    """turn a string into a dictionary"""
    d=defaultdict(int)
    ls = s.split("\n")
    rs = [[int(float(x)*2) for x in  l.split()] for l in ls]
    for i,ll in enumerate(rs):
        for j,x in enumerate(ll):
            if x!=0: d[i+1j*j]=x
    return d

touches0 = set(d)
vs = {}
for pt in touches0:
    vs[pt] = Bool(str(pt)) #be careful here: str(-0j)!=str(0j)!=str(0)

dg1 = getd(g1)
dg2 = getd(g2)
dres=defaultdict(int)

def slv(a,b=0,count=1):
    o = Optimize()
    global d2
    d2=defaultdict(int)
    for x in l:
        if all(y in touches0 for y in x):
            o.add(Or(*map(vs.get,x)))
            for y in x: d2[y]+=1
    z = 0
    tot = 0
    for pt in touches0:
        k = dg1[pt]*a +dg2[pt]*b # +(b if pt in [str(x) for x in ker] else 0) + c
        z = z + If(vs[pt],k,0)
        tot +=k
        d2[pt]=k
    print("kernel")
    for x in range(sz):
        for y in range(sz):
            p = x+y*1j
            print(str(d2[p]).rjust(3),end="")
        print()
    # o.add(Not(vs[str(0j)])) # <- 8
    #for pt in touches0:
    #    z = z + If(vs[pt],1,0)
    mn = o.minimize(z)
    print("checking")
    ii=0
    while ii<count and o.check() == sat:
        ii+=1
        val = mn.value()
        if ii<2: print (val)
        iv = int(str(val))
        g = math.gcd(iv,tot)
        if ii<2: print("implicit density", val/tot,f"{iv//g}/{tot//g}", iv/tot)
        m = o.model()
        for x in range(sz):
            for y in range(sz):
                p = x+y*1j
                if m[vs[p]]:
                    if ii<2:print("#",end="")
                    dres[norm(p)]+=1
                #elif p=="0j": print("*",end="")
                elif ii<2:print(".",end="")
            if ii<2: print()
    return (iv/tot)

r = slv(1,b=0,count=1)

if False:# __name__=="__main__":
    r=0
    while r<1/4:
        dres=defaultdict(int)
        r = slv(1,count=10)
        mx = 0
        mxvs = []
        print("dres:")
        for x in range(sz):
            for y in range(sz):
                p = x+y*1j
                v=dres[p]
                if v/counts[norm(p)]>mx:
                    mx=v/counts[norm(p)]
                    mxvs=[p]
                elif mx==v/counts[norm(p)]:
                    mxvs.append(p)
                print(str(dres[p]).rjust(3),end="")
            print()
        print(len(mxvs))
        for p in crange(box):
            if norm(p) in mxvs:
                d[p]+=1 # 8//counts[norm(p)]#not the right place
                

