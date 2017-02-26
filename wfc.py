import numpy as np
from collections import Counter, deque
from scipy.ndimage import imread
from scipy.ndimage.morphology import binary_erosion



class SimpleTiledModel(object):
    def __init__(self, N, FMX, FMY, initial):
        self.N = int(N)
        self.FMX = int(FMX)
        self.FMY = int(FMY)

        initial, weights = zip(*sorted(Counter(initial).items()))

        self.initial = initial
        self.weights = weights
        self.initial_arr = np.array(initial)
        self.offsets = make_offsets(N)
        self.table = compatibility_table(initial, self.offsets)


        num_colors = len(np.unique(self.initial))
        self.M = len(initial)
        self.wave = np.ones((self.M, FMY, FMX), dtype=bool)

    def update(self):
        """Do one round of collapse and propagate. The point to collapse is
        chosen by minimum entropy. Propagation eliminates incompatible 
        coefficients at the neighbors of the collapsed point. 

        If the coefficients at a point change, we keep propagating to all neighbors 
        that are not collapsed, in a breadth-first fashion. 
        """

        wavesum = self.wave.sum(axis=0)
        all_false = zip_where(wavesum==0)
        if all_false:
            raise(ValueError('these are all false:' + str(all_false)))
        if (wavesum>1).sum() == 0:
            print 'all points collapsed'
            return True

        i,j = self.collapse()

        print 'collapsed',(i,j)

        return self.propagate_from([(i,j)])

    def propagate_from(self, points):

        queue = deque(points)
        modified = 0
        while queue:
            (i,j) = queue.popleft()
            neighbors = self.get_neighbors((i,j))
            neighbors = [(a,b) for a,b in neighbors if sum(self.wave[:,a,b])>1]
            for a,b in neighbors:
                changed = self.propagate((i,j), (a-i, b-j))
                if changed:
                    queue.append((a,b))
                    modified += 1
                    if sum(self.wave[:,a,b]) == 0:
                        assert False

        print 'propagation modified %d points' % modified


    def get_neighbors(self, (i,j)):
        """Get the neighbors to a point. Returns coordinates, not offsets. 

        Could implement periodic or other boundary conditions.
        """
        arr = []
        for a,b in self.offsets:
            a += i
            b += j
            if (a >= 0       and b >= 0 and 
                    a < self.FMY and b < self.FMX):
                arr += [(a,b)]
        return arr

    def collapse(self):
        """Find the point with lowest entropy and pick a single coefficient.
        """
        i,j = self.find_lowest_entropy()
        coefficients = self.wave[:,i,j]
        winner = pick_a_true(coefficients, weights=self.weights)
        self.wave[:,i,j] = False
        self.wave[winner,i,j] = True
        return i,j
         
    def propagate(self, coordinates, offset):
        """Harmonize tile t1 at coordinates, and neighboring tile t2 at 
        coordinates plus offset. Returns a bool indicating whether the wave
        changed.

        Each coefficient in the wave at t2 must be compatible with at least one 
        coefficient at t1, or it is set to false.
        """

        i,j = coordinates
        a,b = offset
        
        # coeffs1 = self.wave[:,i,j].copy()
        # coeffs2 = self.wave[:,a+i,b+j].copy()
        coeffs1 = np.where(self.wave[:,i,j])[0]
        coeffs2 = np.where(self.wave[:,a+i,b+j])[0]

        # t1s = [self.initial[k] for k,x in enumerate(coeffs1) if x]

        
        # print (i,j), coeffs1
        # print (i+a,j+b), coeffs2

        changed = False
        for k2 in coeffs2:
                
            agree = any(self.table[((a,b), k1, k2)] for k1 in coeffs1)

            # print 'set ix', (k,a+i,b+j), 'from', self.wave[k,a+i,b+j], 'to', agree
            
            if not agree:
                self.wave[k2,a+i,b+j] = False
                changed = True

        return changed



    def get_entropy(self):
        """Given (M x X x Y) boolean wavefunction and initial states, return
        an (X x Y) array of entropy at each position.
        """
        entropy = self.wave[0] * 0.
        for i in range(self.FMY):
            for j in range(self.FMX):
                coefficients = self.wave[:,i,j]
                nnz = sum(coefficients)
                colors = self.initial_arr[coefficients].reshape(nnz, self.N**2)
                colors = colors.transpose(1,0)
                # only take the middle
    #             colors = colors[[(self.N**2)/2],:]
                for pixel_coefs in colors:
                    entropy[i,j] += shannon_entropy(pixel_coefs)

        return entropy

    def find_lowest_entropy(self):
        """Return the index of the first position with the lowest entropy. 
        Positions with zero entropy (i.e., collapsed) are ignored.
        """
        H = self.get_entropy()
        collapsed = H == 0
        x = H + collapsed*99999
        minima = zip_where(x==x.min())
        return minima[np.random.choice(len(minima))]

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
                    if self.wave[k,i,j]:
                        superposition[i:i+N, j:j+N] += self.initial[k]
                        weights[i:i+N, j:j+N] += 1


        superposition = superposition / weights.astype(float)

        return superposition


# UTILITIES

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
        return np.random.choice(index)
    weights = np.array(weights) * 1.
    weights = weights[index]
    return np.random.choice(index, p=weights/weights.sum())

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

def compatibility_table(tiles, offsets):
    """Determine if two tiles are compatible, after applying a (dI, dJ) offset.
    """
    d = {}
    for k1, t1 in enumerate(tiles):
        for k2, t2 in enumerate(tiles):
            for i,j in offsets:
                agree = (subsquare(t1, i, j) == subsquare(t2, -i, -j)).all()
                d[((i,j), k1, k2)] = agree
                    
    return d

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

