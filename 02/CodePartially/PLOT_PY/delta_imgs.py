#!/usr/bin/env python

import numpy as np
from PIL import Image
from math import sqrt

class cField:
    def __init__(this, filename):
        f = open(filename)
        lines = f.readlines()
        f.close()

        if(len(lines[-1].strip()) == 0):
            lines = lines[:-1]

        n = len(lines)
        this.re = np.zeros(n)
        this.im = np.zeros(n)

        i = 0
        for l in lines:
            re, im = l.strip().split(None, 1)
            this.re[i] = float(re.strip())
            this.im[i] = float(im.strip())

            i += 1

        m = int(sqrt(n))
        this.re = this.re.reshape((m,m))
        this.im = this.im.reshape((m,m))

    def to_img(this, filename):
        m = 0
        b = 0
        if(this.re.max() != this.re.min()):
            m = 255.0/(this.re.max() - this.re.min())
            b = -m*this.re.min()

        re8 = this.re * m + b
        re8 = re8.astype('uint8')

        m = 0
        b = 0
        if(this.im.max() != this.im.min()):
            m = 255.0/(this.im.max() - this.im.min())
            b = -m*this.im.min()

        im8 = this.im * m + b
        im8 = im8.astype('uint8')
        m = this.re.shape[0]

        img = np.zeros((m, m, 3), dtype='uint8')
        img[:,:,0] = re8
        img[:,:,1] = im8

        rbg = Image.fromarray(img, 'RGB')
        print(filename)
        rbg.save(filename)

if __name__ == '__main__':
    from glob import glob
    files = glob('/p/work1/rasmith/Spintronics/134272/delta-*.txt')
    files = files + glob('latest/??-n.txt')
    files.sort()

    for f in files:
        print(f)
        field = cField(f)
        fn = f[:-4] + ".png"
        field.to_img(fn)
