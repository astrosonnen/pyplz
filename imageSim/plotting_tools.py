import numpy as np
from PIL import Image, ImageDraw, ImageFont


def marshall15_pil_format(images, scales=(1., 1., 1.), alpha=1., Q=1.):

    r = images[0] * scales[0]
    g = images[1] * scales[1]
    b = images[2] * scales[2]

    I = (r + g + b)/3.
    f = np.arcsinh(alpha*I*Q)

    R = r*f/(I*Q)
    G = g*f/(I*Q)
    B = b*f/(I*Q)

    M = max(R.max(), G.max(), B.max())

    R = R/M * 255.
    G = G/M * 255.
    B = B/M * 255.

    R[R<0.] = 0.
    G[G<0.] = 0.
    B[B<0.] = 0.

    flatlist = []
    for img in [R, G, B]:
        img = np.uint8(img.round())
        img = np.flipud(img)
        flatlist.append(img.flatten())

    l = []
    for i in range(images[0].size):
        l.append((flatlist[0][i], flatlist[1][i], flatlist[2][i]))

    return l


def lupton04_pil_format(images, scales=(1., 1., 1.), beta=1.):

    r = images[0] * scales[0]
    g = images[1] * scales[1]
    b = images[2] * scales[2]

    I = (r + g + b)/3.
    f = np.arcsinh(I/beta)

    R = r*f/I
    G = g*f/I
    B = b*f/I

    M = max(R.max(), G.max(), B.max())

    R = R/M * 255.
    G = G/M * 255.
    B = B/M * 255.

    R[R<0.] = 0.
    G[G<0.] = 0.
    B[B<0.] = 0.

    flatlist = []
    for img in [R, G, B]:
        img = np.uint8(img.round())
        img = np.flipud(img)
        flatlist.append(img.flatten())

    l = []
    for i in range(images[0].size):
        l.append((flatlist[0][i], flatlist[1][i], flatlist[2][i]))

    return l


def make_crazy_pil_format(data, cuts):

    newlist = []
    for i in range(0, 3):
        flatd = np.flipud(data[i]).flatten()
        flatd[flatd<0.] = 0.
        flatd *= 255./cuts[i]
        flatd[flatd>255.] = 255.
        flatd = np.uint8(flatd.round())
        newlist.append(flatd)

    l = []

    for i in range(0, data[0].size):
        l.append((newlist[0][i], newlist[1][i], newlist[2][i]))

    return l


def make_one_rgb(sci, light_model, source_model, cuts=(99., 99., 99.)):

    auto_cuts = []
    data = []
    lensresid = []
    fullresid = []
    fullmodel = []

    i = 0

    ncol = 5

    for i in range(3):
        data.append(sci[i])
        cut = np.percentile(sci[i], cuts[i])
        auto_cuts.append(cut)

        fullmodel.append(light_model[i] + source_model[i])
        lensresid.append(sci[i] - light_model[i])
        fullresid.append(sci[i] - light_model[i] - source_model[i])

    dlist = make_crazy_pil_format(data, auto_cuts)
    slist = make_crazy_pil_format(source_model, auto_cuts)
    rlist = make_crazy_pil_format(lensresid, auto_cuts)
    fmlist = make_crazy_pil_format(fullmodel, auto_cuts)
    frlist = make_crazy_pil_format(fullresid, auto_cuts)

    s = (data[0].shape[1], data[0].shape[0])
    dim = Image.new('RGB', s, 'black')
    rim = Image.new('RGB', s, 'black')
    sim = Image.new('RGB', s, 'black')
    fmim = Image.new('RGB', s, 'black')
    frim = Image.new('RGB', s, 'black')

    dim.putdata(dlist)
    sim.putdata(slist)
    rim.putdata(rlist)
    fmim.putdata(fmlist)
    frim.putdata(frlist)

    im = Image.new('RGB', (ncol*data[0].shape[0], data[0].shape[1]), 'black')

    im.paste(dim, (0, 0,))
    im.paste(fmim, (1*s[1], 0))
    im.paste(rim, (2*s[1], 0))
    im.paste(sim, (3*s[1], 0))
    im.paste(frim, (4*s[1], 0))

    return im

def make_full_rgb(sci_list, light_list, source_list, outname='rgb.png'):

    rgbsets = []
    ntotbands = len(sci_list)
    if ntotbands == 1:
        rgbsets.append((0, 0, 0))
    elif ntotbands == 2:
        rgbsets.append((1, 1, 0))
        rgbsets.append((0, 0, 0))
        rgbsets.append((1, 1, 1))
    elif ntotbands == 3:
        rgbsets.append((2, 1, 0))
        rgbsets.append((0, 0, 0))
        rgbsets.append((1, 1, 1))
        rgbsets.append((2, 2, 2))
    else:
        nsets = ntotbands - 2
        for i in range(nsets):
            rgbsets.append((i+2, i+1, i))
        for i in range(ntotbands):
            rgbsets.append((i, i, i))

    nsets = len(rgbsets)

    s = (5*sci_list[0].shape[1], sci_list[0].shape[0])
    fullim = Image.new('RGB', (s[0], nsets*s[1]), 'black')

    for i in range(nsets):
        sci_here = []
        light_here = []
        source_here = []
        for ind in rgbsets[i]:
            sci_here.append(sci_list[ind])
            light_here.append(light_list[ind])
            source_here.append(source_list[ind])
        im = make_one_rgb(sci_here, light_here, source_here)

        fullim.paste(im, (0, i*s[1]))

    fullim.save(outname)

