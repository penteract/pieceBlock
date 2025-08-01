from math import pi, log, sqrt
from numpy import array, hstack, vstack, clip, histogram
from numpy.fft import fft2,fftshift,ifft2
from PIL import Image
import numpy as np
from colorsys import hls_to_rgb

def blur(a):
    kernel = np.array([[1.0,2.0,1.0], [2.0,4.0,2.0], [1.0,2.0,1.0]])
    kernel = kernel / np.sum(kernel)
    arraylist = []
    for y in range(3):
        temparray = np.copy(a)
        temparray = np.roll(temparray, y - 1, axis=0)
        for x in range(3):
            temparray_X = np.copy(temparray)
            temparray_X = np.roll(temparray_X, x - 1, axis=1)*kernel[y,x]
            arraylist.append(temparray_X)

    arraylist = np.array(arraylist)
    arraylist_sum = np.sum(arraylist, axis=0)
    return arraylist_sum

def value_diapason(x, percent=0.95, nbins=100):
    """Use histogram to determine interval, covering 95% of values"""
    counts, bins = histogram(x.ravel(),nbins)
    total = sum(counts)
    accum = 0
    low = bins[-1]
    high = bins[0]
    #enumerate histogram bins starting from the most populated. 
    for i, cnt in sorted(enumerate(counts), 
                          key = (lambda i_c: i_c[1]),
                          reverse=True):
        accum += cnt
        low = min(low, bins[i])
        high = max(high, bins[i+1])
        if accum > percent * total:
            break
    return low, high
    

def toimage(fimg, gamma=1., percent=0.95, extend = 1.1, save=None):
    """Show binary matrix as monochrome image, automatically detecting upper and lower brightness bounds
    """
    low, high = value_diapason(fimg, percent=percent)
    print(low,high)

    mid = (low+high)/2
    low = mid + (low-mid)*extend
    high = mid + (high-mid)*extend
    low=max(low,fimg.min())
    high=min(high,fimg.max())
    print(low,high)
    #global tmp
    tmp = (clip((fimg-low)/(high-low),0,1)**gamma*255).astype(np.uint8)
    tmp = np.array([tmp,tmp,tmp])
    tmp = tmp.swapaxes(0,2)
    image = Image.fromarray(tmp, "RGB")
    if save is not None:
        image.save(save)
        print("Saving image file: {}".format(save))
    return image

def toimagecol(fimg, save=None):
    """save a grid of complex numbers as an image file"""
    colorize(fimg)
    image = Image.fromarray(colorize(fimg).astype(np.uint8), "RGB")
    if save is not None:
        image.save(save)
        print("Saving image file: {}".format(save))
    return image

def colorize(z):
    """colors a grid of complex numbers so that hue indicates argument and brightness indicates modulus"""
    r = np.abs(z)
    high = np.percentile(r,97.5)
    r = clip(r/high,0,1)**1
    arg = np.angle(z)

    #arg,r = np.log(r),arg+pi
    h = (arg )  / (2 * pi)
    l = (r)/2 # 1.0 - 1.0/(1.0 + r**0.1)
    #l=r/2
    s = 0.5

    c = np.vectorize(hls_to_rgb) (h,l*255,s) # --> tuple
    c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
    c = c.swapaxes(0,2)
    return c

def fourimage(a,b,c,d,save=None):
    cors = []
    for im in [a,b,c,d]:
        cors.append(colorize(im))
    r1 = np.concatenate((cors[0],cors[1]))
    r2 = np.concatenate((cors[2],cors[3]))
    r = np.concatenate((r1,r2),axis=1).astype(np.uint8)
    image = Image.fromarray(r, "RGB")

    if save is not None:
        image.save(save)
        print("Saving image file: {}".format(save))
    return image
def embiggen(m,k=2):
    sh2 = tuple(x*k for x in m.shape)
    r = np.zeros(sh2)
    for a in range(m.shape[0]):
        for b in range(m.shape[1]):
            r[a*k:a*k+k,b*k:b*k+k] = m[a,b]
    return r
def takeHalfandRotate(pattern,parity=1):
    """Given a square grid tiling the plane,
    takes pixels with odd coordinates and return a square grid the same size
    containing those points rotatad 45 degrees about the top left corner
    again representing a tiling"""
    sh = pattern.shape
    a,b = np.indices(sh)
    #a-=1
    return pattern[(a-b)%sh[0],(a+b+parity)%sh[1]]

def manyimage(tops,bots,row3=None,save=None):
    t=iter(tops)
    tr=colorize(next(t))
    col = np.array([[[0,255,0]]]*tr.shape[0] )
    for im in t:
        #print("p",tr,"i",im)
        # tr = np.concatenate((tr, colorize(im)),axis=1)
        tr=np.concatenate((tr,col, colorize(im)),axis=1)
    b=iter(bots)
    br=colorize(next(b))
    for im in b:
        br=np.concatenate((br,col, colorize(im)),axis=1)
    if row3 is not None:
        bb=iter(row3)
        bbr=colorize(next(bb))
        for im in bb:
            bbr=np.concatenate((bbr,col, colorize(im)),axis=1)
        
    row = np.array([[[0,255,0]]*tr.shape[1]])
    r = np.concatenate((tr,row,br)).astype(np.uint8)
    if row3 is not None:
        r = np.concatenate((r,row,bbr)).astype(np.uint8)
    image = Image.fromarray(r, "RGB")

    if save is not None:
        image.save(save)
        print("Saving image file: {}".format(save))
    return image

if __name__=="__main__":
    import fourier_pics_elephant as e
    def tomat(grd):
        m = np.zeros((e.sz,e.sz))
        for p in e.crange(e.box):
            m[int(p.real),int(p.imag)] = int(grd[p])
        return m
    #r = e.slv(1)
    import random
    D = tomat(e.dg1)
    tops = [D,tomat(e.dg2),tomat(e.diamond)]
    for i in range(0):
        #print(tops[-1])
        tops.append(tomat(random.choice(e.pats)))
    bots = [ fftshift(fft2(p))  for p in tops]
    bbots = [np.log(np.abs(x)+0.0000001) for x in bots]
    ks = [ fftshift(fft2(tomat(p)))  for p in e.pats]
    #for b in bots[1:]: assert (np.sum(np.abs(bots[0]*b)>.00000001))==1
    #tops.append(np.sum(np.abs(ks[1:]),axis=0))
    #bots.append(ifft2(fftshift(tops[-1])))
    #bots[1] = np.log(np.abs(bots[0])+0.0000001)
    bg=lambda x: embiggen(x,4)
    manyimage(map(bg,tops),map(bg,bots),map(bg,bbots),save="testElephantPats_{}_{}_quad.png".format(e.sz,hex(random.randrange(0x1000))[2:] ))
        
    #assert r==1/4
    #D2 = tomat(e.dres)
    #bg=lambda x: embiggen(x,4)
    #fourimage(bg(D),bg(fftshift(fft2(D))),bg(D2),bg(fftshift(fft2(D2))),save="testElephantbig_{}_{}_quad.png".format(e.sz,hex(random.randrange(0x1000))[2:] ))
    #fourimage(D,fftshift(fft2(D)),D2,fftshift(fft2(D2)),save="testElephant_{}_quad.png".format(e.sz))
    """
    D = np.zeros(8) #solid_dragon(N,0)#toimage(D).show()
    toimage(D, save="hi_res_backgrounds/dragon_{}.png".format(N))
    f = fft2(D)
    fimg = np.abs(f)#(np.abs(f))#np.log(np.abs(f)+1e-100)#np.real(f)
    #f2 = fimg**0.5 #  / blur(fimg)**0.5
    #fimg[1::2,0::2]*=-1
    #fimg[0::2,1::2]*=-1
    D2 = fft2(fimg)
    #fourimage(D,fftshift(f),fftshift(fimg),D2,save="dragon2_{}_quad.png".format(N))
    toimage(np.abs(D2), save="hi_res_backgrounds/dragon_fft_abs_fft_{}.png".format(N))
    #toimage((np.abs(D2)), save="dragon_diagram2_fft_fft{}.png".format(N))
    #toimage(fftshift(np.abs(fimg)), save="dragon_diagram2_{}_fft.png".format(N))
    #f2 = takeHalfandRotate(f)
    sh = D.shape
    x,y = np.indices(sh)
    k = np.exp(1j*pi*(x+y)/sh[0])
    toimagecol(fft2(f*(1+1j)/k),"hi_res_backgrounds/dragon_fft_derainbow_fft_{}.png".format(N))
    #fourimage(f2,takeHalfandRotate(fft2(f2),0),f,fft2(f*(1+1j)/k),save="dragon2_{}_quad2.png".format(N))
    """
