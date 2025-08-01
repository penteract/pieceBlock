#elephant3.py
from z3 import *

def getPath (o,dx,dy,xy):
    # o in C
    assert (dx in [1,-1])
    assert (dy in [1j,-1j])
    yield o
    kx = o+(dx if xy else dy)
    for i in range(4):
        yield kx
        kx+=dx+dy


l = []
from collections import defaultdict
d=defaultdict(int)

sz=3
ker = [i+j*1j for i in range(sz) for j in range(sz)]
#ker = [0j,1j,1+0j,1+1j]
#ker=[0j]

#mx = 3

for x in range(-5,6+sz):
    for y in range(-5,6+sz):
        origin = x+y*1j
        for dx in [1,-1]:
            for dy in [1j,-1j]:
                for xy in [True,False]:
                    path = list(getPath(origin,dx,dy,xy))
                    # sorry that everything becomes strings.
                    l.append([str(x) for x in path])
                    if 0 in path:
                        for p in path:
                            # d[str(p)]+=1
                            for k in ker:
                                d[str(p+k)]+=1
                            

touches0 = set(d)
vs = {}
for pt in touches0:
    vs[pt] = Bool(pt)

def slv(a,b=0,c=0):
    o = Optimize()
    global d2
    d2=defaultdict(int)
    for x in l:
        if all(y in touches0 for y in x):
            o.add(Or(*map(vs.get,x)))
            for y in x: d2[y]+=1
    # o.add(vs[str(0j)]) <- 56
    # o.add(vs[str(1j)]) <- 44
    # o.add(vs[str(1+1j)]) <- 44
    z = 0
    tot = 0
    for pt in touches0:
        k = d[pt]*a+(b if pt in [str(x) for x in ker] else 0) + c
        z = z + If(vs[pt],k,0)
        tot +=k
        d2[pt]=k

    # o.add(Not(vs[str(0j)])) # <- 8
    #for pt in touches0:
    #    z = z + If(vs[pt],1,0)
    mn = o.minimize(z)
    assert o.check() == sat
    val = mn.value()
    print (val) # 40
    iv = int(str(val))
    g = math.gcd(iv,tot)
    print("implicit density", val/tot,f"{iv//g}/{tot//g}", iv/tot)
    m = o.model()
    for x in range(-4,4+sz):
        for y in range(-4,4+sz):
            p = str(x+y*1j)
            if m[Bool(p)]:print("#",end="")
            elif p=="0j": print("*",end="")
            else:print(".",end="")
        print()
for k in ker:
    d[str(k)]-=0
d[str(1+1j)]-=1
slv(1,-17,0)

for x in range(-4,4+sz):
    for y in range(-4,4+sz):
        p = str(x+y*1j)
        print(str(d2[p]).rjust(3),end="")
    print()

"""
for a in range(1,5):
    for b in range(1,5):
        print(a,b)
        slv(a,b)

"""


"""
a,b = Bools("a b")
o.add(Or(a, b))
mx = o.minimize(a*2+b*3)
while o.check() == sat:
    print (mx.value())
    break

o.add_soft(a, 3)                                          
print(k:=o.add_soft(b, 4))

# r = o.minimize(a+b)


x, y = Ints('x y')
opt = Optimize()
opt.set(priority='pareto')
opt.add(x + y == 10, x >= 0, y >= 0)
mx = opt.maximize(x)
my = opt.maximize(y)
while opt.check() == sat:
    print (mx.value(), my.value())
"""
