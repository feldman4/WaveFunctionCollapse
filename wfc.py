import numpy as np
from collections import Counter, OrderedDict, deque
from scipy.ndimage import imread
from scipy.ndimage.morphology import binary_erosion


plog2 = [0] + [x*np.log2(x) for x in range(1,1000)]



class SimpleTiledModel(object):
    def __init__(self, N, FMX, FMY, initial, periodic=False):
        self.N = N = int(N)
        self.FMX = FMX = int(FMX)
        self.FMY = FMY = int(FMY)
        self.periodic = periodic

        initial, weights = zip(*sorted(Counter(initial).items()))

        self.M = M = len(initial)
        self.initial = initial
        self.weights = weights
        self.initial_arr = np.array(initial)
        self.offsets = make_offsets(N)
        self.neighbors = make_offset_dict(N, FMX, FMY)
        self.table = compatibility_table(initial, self.offsets, N)
        
        self.wave = np.ones((FMY, FMX, M), dtype=bool)
        self.entropy = self.wave[:,:,0] * -1. # cached version
        self.support = support_table(FMX, FMY, self.table, periodic=periodic)


    def update(self, support=True):
        """Do one round of collapse and propagate. The point to collapse is
        chosen by minimum entropy. Propagation eliminates incompatible 
        coefficients at the neighbors of the collapsed point. 

        If the coefficients at a point change, we keep propagating to all neighbors 
        that are not collapsed, in a breadth-first fashion. 

        Test support:
        from copy import copy, deepcopy
        model0 = deepcopy(model)
        np.random.seed(0); model0.update()
        model1 = deepcopy(model)
        np.random.seed(0); model1.update(support=False)
        """

        wavesum = self.wave.sum(axis=2)
        all_false = zip_where(wavesum==0)
        if all_false:
            raise(ValueError('these are all false:' + str(all_false)))
        if (wavesum>1).sum() == 0:
            print 'all points collapsed'
            return True

        i,j,delta_wave = self.collapse()

        # print 'collapsed',(i,j)
        if support:
            return self.propagate_from_support((i,j), delta_wave)
        return self.propagate_from([(i,j)])

    def propagate_from_support(self, point, delta_wave):
        if not delta_wave.any():
            return

        queue = OrderedDict()
        queue[point] = delta_wave

        while queue:
            (i,j), delta_wave = queue.popitem(last=False)
            self.entropy[i, j] = -1

            new_tests = self.propagate_support((i,j), delta_wave)
            
            for point, delta_wave in new_tests:
                if point in queue:
                    queue[point] = queue[point] | delta_wave
                else:
                    queue[point] = delta_wave
 
    def propagate_support(self, point, d_wave):
        """ Propagate updates based on loss of support from point A to B, rather
        than incompatibility of B with A. Manages changes in support array.
        """
        i,j = point

        # loss in support
        d_support = self.table[:, :, d_wave, :].sum(axis=2) 

        # all gone?
        new_d_waves = self.support[:, :, i, j, :] == d_support 
        new_d_waves &= (d_support > 0)

        self.support[:, :, i, j, :] -= d_support 
        # assert (self.support>=0).all()

        # for each neighbor, pass on coefficients that totally lost support
        # need to ignore offset (0,0)
        neighbors, offsets = self.get_neighbors((i,j))
        queue = []
        for (i2, j2), (a,b) in zip(neighbors, offsets):

            # get the coefficients that were zeroed
            new_d_wave = new_d_waves[a, b, :] & self.wave[i2,j2]

            # breadth-first
            if any(new_d_wave):
                self.wave[i2, j2, new_d_wave] = False
                queue.append(((i2, j2), new_d_wave))

        return queue

    def propagate_from(self, points):
        """Propagation strategy is to check all neighbors of a point if the 
        point is changed. 
        """
        # store points in a deque set
        queue_set = set(points) 
        queue = deque(queue_set)

        modified = 0
        touched = 0
        deleted = 0
        while queue:
            (i,j) = queue.popleft()
            queue_set.remove((i,j))
            neighbors_offsets = zip(*self.get_neighbors((i,j)))
            # fails to identify contradictions
            # neighbors_offsets = [x for x in neighbors_offsets 
            #                        if self.wave[x[0][0],x[0][1],:].sum()>1]
            
            for ((i2, j2), (a,b)) in neighbors_offsets:
                touched += 1
                changed = self.propagate((i,j), (i2,j2), (a,b))
                if changed:
                    if (i2,j2) not in queue_set:
                        queue.append((i2,j2))
                        queue_set.add((i2,j2))
                    modified += 1
                    deleted += changed
                    if self.wave[i2,j2,:].sum() == 0:
                        assert False

        # print 'propagation deleted %d coefficients in %d points; touched %d points' % (deleted, modified, touched)

    def propagate(self, coordinates1, coordinates2, offset):
        """Harmonize tile t1 at coordinates, and neighboring tile t2 at 
        coordinates plus offset. Returns a bool indicating whether the wave
        changed. 

        Each coefficient in the wave at t2 must be compatible with at least one 
        coefficient at t1, or it is set to false. Modifies wave entry for t2 only.
        """

        i,j = coordinates1
        i2,j2 = coordinates2
        a,b = offset
        
        coeffs1 = self.wave[i,j,:].nonzero()[0]
        coeffs2 = self.wave[i2,j2,:].nonzero()[0]

        # ~2/3 of propagation time spent in this line
        agree = self.table[a, b][np.ix_(coeffs1, coeffs2)].any(axis=0)

        self.wave[i2,j2,coeffs2] = agree
        changed = (~agree).sum()

        if changed:
            self.entropy[i2,j2] = -1 # flag to recalculate

        return changed

    def get_neighbors(self, (i,j)):
        """Get the neighbors to a point. Returns neighbors and offsets. 

        In the periodic case the offset is the original offset before applying
        boundary conditions.
        """
        if self.periodic:
            neighbors = [((i + a) % self.FMY, (j + b) % self.FMX) 
                            for a,b in self.offsets]
            return neighbors, self.offsets

        arr = []
        for a,b in self.offsets:
            if (a + i >= 0       and b + j >= 0 and 
                a + i < self.FMY and b + j< self.FMX):
                arr += [[(i+a, j+b), (a,b)]]

        return zip(*arr)

    def collapse(self):
        """Find the point with lowest entropy and pick a single coefficient.
        """
        i,j = self.find_lowest_entropy()
        coefficients = self.wave[i,j,:].copy()
        winner = pick_a_true(coefficients, weights=self.weights)
        self.wave[i,j,:] = False
        self.wave[i,j,winner] = True
        self.entropy[i,j] = -1 # flag to recalculate
        delta_wave = coefficients & ~self.wave[i,j,:]
        return i,j,delta_wave  

    def calc_entropy(self, coordinates):
        i,j = coordinates
        entropy = 0

        # coefficients = self.wave[:,i,j]        
        # colors = self.initial_arr[coefficients].reshape(-1, self.N**2)

        coefficients = self.wave[i,j,:]
        nnz = coefficients.sum()
        colors = self.initial_arr[coefficients].reshape(nnz, self.N**2)
        colors = colors.transpose(1,0)

        

        for pixel_coefs in colors: 
            for c in np.bincount(pixel_coefs):
                entropy -= plog2[c]
            # entropy += shannon_entropy(pixel_coefs)  

        entropy = entropy/nnz + plog2[nnz]
        return entropy

    def get_entropy(self, recalculate=True):
        """Given (M x X x Y) boolean wavefunction and initial states, return
        an (X x Y) array of entropy at each position.
        """
        if recalculate:
            self.entropy[:] = -1

        indices = zip_where(self.entropy == -1)

        for i,j in indices:
            self.entropy[i,j] = self.calc_entropy((i,j))

        return self.entropy

    def find_lowest_entropy(self):
        """Return the index of the first position with the lowest entropy. 
        Positions with zero entropy (i.e., collapsed) are ignored.
        """
        H = self.get_entropy(recalculate=False)
        collapsed = H == 0
        x = H + collapsed*99999
        minima = zip_where(x==x.min())
    
        # return minima[np.random.choice(len(minima))]
        return minima[choice(len(minima))]

    def get_superposition(self, color_map=None):
        N = self.N
        color_null = 0
        superposition = [[color_null for _ in range(self.FMX + N - 1)] 
                                     for _ in range(self.FMY + N - 1)]
        superposition = np.array(superposition)
        weights = np.zeros_like(superposition, dtype=int)
        

        for i in range(self.FMY):
            for j in range(self.FMX):
                for k in range(self.M):
                    if self.wave[i,j,k]:
                        superposition[i:i+N, j:j+N] += self.initial[k]
                        weights[i:i+N, j:j+N] += 1


        superposition = superposition / weights.astype(float)

        return superposition

    def complete(self, support=True):
        i = 0

        while (self.wave.sum(axis=2) > 1).any():
            self.update(support=support)
            i += 1

        if (self.wave.sum(axis=2) != 1).any():
            print 'contradiction after', i, 'updates'
            return None

        return i



# UTILITIES

def memodict(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret 
    return memodict().__getitem__

def rotate(tile, cardinality):

    arr = [tile]

    if cardinality == 2:
        arr = [tile, np.rot90(tile)]

    if cardinality == 4:
        arr = [tile]
        for i in range(1,4):
            arr += [np.rot90(arr[-1])]

    return [array_to_tuple(t) for t in arr]

def array_to_tuple(array):
    return tuple(tuple(x) for x in array)

def initial_from_png(filename, N):
    img = imread(filename)
    if img.ndim == 2:
        img = img[..., None]
    channels = img.shape[-1]
    colors = sorted({tuple(v) for v in img.reshape(-1,channels)})
    colors = {v:i for i,v in enumerate(colors)}
    arr = np.apply_along_axis(lambda v: colors[tuple(v)], 2, img)
    I,J = arr.shape
    tiles = []
    for i in range(I - N + 1):
        for j in range(I - N + 1):
            tile = array_to_tuple(arr[i:i+N, j:j+N])
            tiles += [tile]
    return tiles

def border_points(arr):
    
    border = ~binary_erosion(arr) & arr
    return zip_where(border)

def shannon_entropy(xs):
    """Compute Shannon entropy of a sequence, in bits.
    """
    xs = list(xs)
    N = float(len(xs))
    entropy = 0
    for n in set(xs):
        p = xs.count(n) / N
        entropy += -p * np.log2(p)
    return entropy

def pick_a_true(bools, weights=None):
    """Choose an index corresponding to a true value at random, from a 
    list of true and false values.
    """
    index = np.where(bools)[0]
    if weights is None:
        # return np.random.choice(index)
        return choice(index)
    weights = np.array(weights) * 1.
    weights = weights[index]
    # return np.random.choice(index, p=weights/weights.sum())
    return choice(index, p=weights/weights.sum())

def make_offsets(N):
    """Produce a list of (dI, dJ) offsets for a tile of width N.
    """
    offsets = []
    for i in range(N):
        i -= N/2 
        for j in range(N):
            j -= N/2
            if i or j:
                offsets += [(i,j)]
    return offsets

def make_offset_dict(N, FMX, FMY, periodic=False):
    offsets = make_offsets(N)
    offset_dict = {}
    for i in range(FMY):
        for j in range(FMX):

            arr = []
            for a,b in offsets:
                a += i
                b += j
                if periodic:
                    arr += [(a,b)]
                elif  (a >= 0  and b >= 0 and 
                       a < FMY and b < FMX):
                    arr += [(a,b)]
            offset_dict[(i,j)] = arr
    return offset_dict    

def compatibility_table(tiles, offsets, N):
    """Determine if two tiles are compatible, after applying a (dI, dJ) offset.
    If offsets doesn't include (0,0), table[0,0,:,:] will be all zero.
    """
    M = len(tiles)
    table = np.zeros((N, N, M, M), dtype=bool)
    for k1, t1 in enumerate(tiles):
        for k2, t2 in enumerate(tiles):
            for i,j in offsets:
                agree = (subsquare(t1, i, j) == subsquare(t2, -i, -j)).all()
                table[i, j, k1, k2] = agree
                    
    return table

def support_table(FMX, FMY, table, periodic=True):
    """Describe initial support along edges between points.
    support[offsetI, offsetJ, positionI, positionJ, stateI]
    """

    # sum over origin state to get amount of support in the destination state
    initial_support = table.sum(axis=2).astype(int)
    initial_support = initial_support[:, :, None, None, :]
    support = np.tile(initial_support, [1, 1, FMY, FMX, 1])

    if not periodic:
        support[-1, :, 0, :, :] = 0
        support[:, -1, :, 0, :] = 0
        support[1, :, -1, :, :] = 0
        support[:, 1, :, -1, :] = 0

    return support

def zip_where(arr):
    """
    """
    return zip(*np.where(arr))

def subsquare(tile, i, j):
    """Return a rectangle within tile offset by (i,j) and clipped at
    the edges.
    """
    tile = np.array(tile)
    if i < 0:
        tile = tile[:i,:]
    else:
        tile = tile[i:,:]
    if j < 0:
        tile = tile[:,:j]
    else:
        tile = tile[:,j:]
    return tile

def choice(*args, **kwargs):
    # return np.random.choice(*args, **kwargs)
    return choice_(*args, **kwargs)

def choice_(xs, p=None):
    if isinstance(xs, int):
        xs = range(xs)
    if not (p is None):
        cs = np.cumsum(p)
        n = getRandom() * max(cs)
        for i,p_ in enumerate(cs):
            if p_ > n:
                break 
        return xs[i]
    return xs[int(getRandom() * len(xs))]

ticker = [0]
def getRandom(ticker=ticker):
    ticker[0] += 1
    i = ticker[0]
    return randos[i % len(randos)]

checkerboardA = \
    ( ( 0, 1, 0 )
    , ( 1, 0, 1 )
    , ( 0, 1, 0 )
    )

checkerboardB = \
    ( ( 1, 0, 1 )
    , ( 0, 1, 0 )
    , ( 1, 0, 1 )
    )

checkerboardC = \
    ( ( 1, 0, 0 )
    , ( 0, 1, 1 )
    , ( 1, 0, 0 )
    )

checkerboardD = \
    ( ( 0, 0, 1 )
    , ( 1, 1, 0 )
    , ( 0, 0, 1 )
    )

checkerboardE = \
    ( ( 0, 1, 1 )
    , ( 1, 0, 0 )
    , ( 0, 1, 1 )
    )

checkerboardF = \
    ( ( 1, 1, 0 )
    , ( 0, 0, 1 )
    , ( 1, 1, 0 )
    )

checkerboardABCDEF = \
    [ checkerboardA
    , checkerboardB
    , checkerboardC
    , checkerboardD
    , checkerboardE
    , checkerboardF
    ]


randos= [3.57095269e-01,   6.16476608e-01,   4.42510505e-01,
         5.75952352e-01,   8.48412660e-02,   2.23283599e-01,
         2.17890220e-01,   8.25591994e-01,   4.60241400e-01,
         5.45692371e-01,   7.87609761e-01,   7.39860852e-01,
         7.20097577e-01,   7.97826275e-01,   8.81896387e-01,
         8.21604693e-01,   3.05941324e-01,   7.45375147e-01,
         3.88219250e-01,   4.69667654e-01,   7.10044376e-01,
         3.65126628e-01,   9.27250328e-01,   7.32187105e-01,
         2.43837089e-01,   6.79664391e-01,   5.95789013e-02,
         8.56288178e-01,   4.91198192e-01,   5.32594697e-01,
         1.51655961e-01,   8.04582107e-01,   6.13809025e-01,
         2.49226914e-01,   3.95723003e-01,   5.60075997e-01,
         5.55833098e-01,   5.52143342e-01,   4.62360874e-01,
         6.70656181e-01,   8.83678959e-01,   4.79704318e-01,
         4.16099641e-01,   2.51388180e-01,   9.46254162e-01,
         9.83089900e-01,   5.51854839e-01,   5.49876400e-01,
         4.63910645e-01,   3.85484346e-01,   4.07507369e-01,
         2.85149566e-01,   8.52476098e-04,   4.34611836e-01,
         5.87336538e-01,   5.22781746e-01,   6.81185355e-01,
         8.81993566e-02,   6.29008279e-01,   7.79468082e-01,
         1.96662091e-01,   4.81846389e-01,   3.76311610e-01,
         4.76825510e-02,   3.16277339e-01,   6.49173972e-01,
         2.56410886e-01,   4.53257301e-01,   1.30402898e-01,
         8.79967400e-01,   9.91131913e-01,   9.53496921e-02,
         7.80540215e-01,   3.18601116e-01,   4.35788740e-01,
         6.61891353e-01,   3.22648602e-01,   2.45114552e-02,
         7.01992141e-01,   7.96804980e-01,   1.22742024e-01,
         8.50894394e-01,   9.07070546e-01,   1.97922313e-01,
         8.70998784e-02,   8.67761386e-01,   8.14737827e-01,
         3.03883303e-01,   2.56563996e-02,   7.59538650e-01,
         6.38786620e-01,   5.64430503e-01,   9.09495514e-01,
         6.76696644e-01,   2.83605500e-01,   4.49995756e-01,
         9.92281843e-01,   2.53107448e-01,   8.84983935e-01,
         3.44674163e-01]