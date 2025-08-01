from z3 import *

def getPath (o,dx,dy,xy):
    """find the specified path a long janggi e0lephant would take starting at o"""
    # o is a complex number
    assert (dx in [1,-1])
    assert (dy in [1j,-1j])
    yield o
    p = o+(dx if xy else dy)
    for i in range(4):
        yield p
        p+=dx+dy


l = []
from collections import defaultdict

# this is our kernel. We are asking Z3 to calculating the configuration of
# blocking squares which cover the minimal weight of our kernel
d=defaultdict(int) 

sz=1

# just generate all the paths we might care about. This isn't the slow part.
for x in range(-4,3+sz):
    for y in range(-4,3+sz):
        origin = x+y*1j
        aa = min(x+4,3+sz-1-x,3)
        bb = min(y+4,3+sz-1-y,3)
        # If you understand why the below produces an optimal kernel
        # you know more than I do at the time of writing
        d[str(origin)]=((aa+1)*(bb+1)- ((x==0 or x==-1) and (y==0 or y==-1)))*2
        if d[str(origin)]==6*2: d[str(origin)]-=1
        for dx in [1,-1]:
            for dy in [1j,-1j]:
                for xy in [True,False]:
                    path = list(getPath(origin,dx,dy,xy))
                    # sorry that everything becomes strings.
                    l.append([str(x) for x in path])
grid = set(d) # the squares we care about
vs = {}
for pt in grid:
    vs[pt] = Bool(pt) # this variable is true if the cell at p is blocked

o = Optimize()
d2=defaultdict(int)

for x in l:
    if all(y in grid for y in x):
        o.add(Or(*map(vs.get,x))) # at least 1 cell on each path must be blocked

#o.add(And(*[Not(vs[str((1+1j)*i)]) for i in range(4)]))

# we want to know the configuration of blocked cells
# with minimal weight where they overlap our kernel
z = 0
tot = 0
for pt in grid:
    k = d[pt]
    z = z + If(vs[pt],k,0)
    tot +=k
mn = o.minimize(z)
print("checking")
assert o.check() == sat

# see the results
val = mn.value()
iv = int(str(val)) # there's probably a politer way to do this
g = math.gcd(iv,tot)
print("implicit density", val/tot,f"{iv//g}/{tot//g}", iv/tot)

#prettily print a satisfying model
m = o.model()
for x in range(-4,3+sz):
    for y in range(-4,3+sz):
        p = str(x+y*1j)
        if m[Bool(p)]:
            print("#",end="")
        elif p=="0j": print("*",end="")
        else:print(".",end="")
    print()

#print the kernel
for x in range(-4,3+sz):
    for y in range(-4,3+sz):
        p = str(x+y*1j)
        print(str(d[p]).rjust(3),end="")
    print()

"""
tesseract — Yesterday at 4:30 PM
A lower bound on density is given by minimal_satisfying_value(kernel)/sum(kernel), so we just need a kernel for which this equals 1/4

...

praseodymiumspike — Yesterday at 4:56 PM
Why does that formula give a lower bound?
￼
Matthew Bolan — Yesterday at 5:09 PM
Let chi(x) be the indicator function of the blocked squares. We want to minimize the dot product <chi, 1> of chi(x) with the constant function 1. With * as convolution, <chi * f, 1> = <chi, 1><f,1> (edited)
[5:10 PM]
Now I'm ignoring that things are infinite here, so in practice there are some ignorable boundary terms (edited)
￼
praseodymiumspike — Yesterday at 5:13 PM
ok im gonna take a break and come back to this
￼
@Matthew Bolan
Let chi(x) be the indicator function of the blocked squares. We want to minimize the dot product <chi, 1> of chi(x) with the constant function 1. With * as convolution, <chi * f, 1> = <chi, 1><f,1> (edited)
￼
tesseract — Yesterday at 5:16 PM
Thanks for spelling out the argument formally. It's definitely true, but I was struggling to say precisely why (I'm not very familiar with convolutions).
￼
Matthew Bolan — Yesterday at 5:18 PM
Put a less insane way, first imagine the kernel has only 1s. Then the minimizer and tesseract's formula have a very concrete meaning - in this subset of the board with n squares you need m burnt squares, so we cover m/n of this subset. By placing a copy of this subset at every square in some large region with k squares, we cover the board n times (visiting k*n squares up to boundary effects) and in particular count each burnt square n times, so we found (up to boundary effects) km burnt squares, for a fraction of ~ km / kn = m/n (edited)
[5:19 PM]
And it's just fine to change the weights from 0,1 weights to arbitrary weights
[5:19 PM]
Since when you cover the board like that every square gets weighted in every possible way
￼
Matthew Bolan — Yesterday at 5:26 PM
Anyway really in the limit you just care about large bricks of 1s
[5:27 PM]
My kernel is only there to gain a little extra from the boundary
[5:27 PM]
When you convolve it with a brick of 1s you'll get a bit more I suspect, but as the brick gets bigger this effect will go away
￼
@praseodymiumspike
ok im gonna take a break and come back to this
￼
tesseract — Yesterday at 5:45 PM
It was a good question to ask.
Here's my version of the argument:
Put a copy of the kernel centered at each square, and give that square a score equal to the sum of kernel entries over blocked squares.
The average score for a square can't be lower than the minimal allowable score for a kernel. If you switch from giving score to the squares where the kernels are centered to giving score to the blocked squares directly, you don't change the total, and each blocked square gets a score of sum(kernel) because there's a kernel overlaid at every possible relative position. This means minimal_allowable_score * number_of_squares <= number_blocked_squares * sum(kernel). 
Rearranging, we get minimal_allowable_score / sum(kernel) <= number_blocked/total_number = density_blocked. (edited)
[5:46 PM]
To handle the infinites I'm thinking about the limit of larger and larger blocks. The boundary is the square root of the size of the center, so it shouldn't affect the asymptotic density.
￼
...
@tesseract
To handle the infinites I'm thinking about the limit of larger and larger blocks. The boundary is the square root of the size of the center, so it shouldn't affect the asymptotic density.
￼
praseodymiumspike — Yesterday at 7:06 PM
Is there some kind of limit you're implicitly taking, then?
￼
@praseodymiumspike
Is there some kind of limit you're implicitly taking, then?
￼
tesseract — Yesterday at 7:20 PM
The limit of "what is the smallest density of blocked squares needed to prevent elephant movement within an nxn square"  as n tends to infinity.(divide the number of blocked squares by n^2 to give a density). The equalities aren't all exact on finite grids, but they get closer as n gets bigger.
"""
