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

sz=5
box = sz+sz*1j
#ker = [i+j*1j for i in range(sz) for j in range(1)]
#ker = [0j,1j,1+0j,1+1j]
ker=[0j]

# gen
for origin in crange(box):
    d[origin]=1
    for dx in [1,-1]:
        for dy in [1j,-1j]:
            for xy in [True,False]:
                l.append(list(getPath(origin,dx,dy,xy)))


touches0 = set(d)
vs = {}
for pt in touches0:
    vs[pt] = Bool(str(pt)) #be careful here: str(-0j)!=str(0j)!=str(0)




def slv(a,count=1):
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
        k = d[pt]*a # +(b if pt in [str(x) for x in ker] else 0) + c
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
    while o.check() == sat and ii<count:
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
                    dres[p]+=1
                #elif p=="0j": print("*",end="")
                elif ii<2:print(".",end="")
            if ii<2: print()
    return (iv/tot)

r=0
while r<1/4:
    dres=defaultdict(int)
    r = slv(1,count=20)
    mx = 0
    mxvs = []
    print("dres:")
    for x in range(sz):
        for y in range(sz):
            p = x+y*1j
            v=dres[p]
            if v>mx:
                mx=v
                mxvs=[p]
            elif mx==v:
                mxvs.append(p)
            print(str(dres[p]).rjust(3),end="")
        print()
    print(len(mxvs))
    for p in mxvs:
        d[p]+=1

