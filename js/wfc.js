// wfc.js

var _ = require('lodash');
var Deque = require("double-ended-queue");
window._ = _;    
window.wfc = exports
// elm-js api choices
// how much state does js hold?
// could describe a wave by a small record like
// {m, n, compatibility[n, n, m, m]}
// this assumes the same states are present 

// initial: array of states. should all be 
exports.SimpleTiledModel = class SimpleTiledModel {

    constructor(N, initial, FMX, FMY) {
        this.offsets = [-FMX-1, -FMX, -FMX+1, -1, 1, FMX-1, FMX, FMX+1]

        this.N = N // number of points
        this.FMX = FMX
        this.FMY = FMY
        this.neighbors = this.getNeighbors(FMX, FMY, this.offsets)
 
        let initialUnique = initial.reduce(superCountBy, [])
        this.initial = initialUnique.map(x => x[0])
        this.weights = initialUnique.map(x => x[1].length)
        this.M = _.size(initialUnique) // number of states
        this.table = this.makeCompatibilityTable(this.M, this.offsets, FMX)

        this.wave = Array(this.N).fill().map(x => (new Uint8Array(this.M)).fill(1))
        this.entropy = new Float32Array(this.N).fill(-1)

        this.visited = 0
    }


    update() {
        if (this.contradiction()) {
            throw new Error('contradiction before update')
        }
        if (this.done()) {
            return
        }
        let {point, state} = this.randomCollapse()
        console.log(point, state)
        if (point === null) {
            return
        }
        return this.propagate_from([point])
     }

    propagate_from(points) {
        let queue_set = new Set(points)
        let queue = new Deque([...queue_set]);

        while (!queue.isEmpty()) {
            let point = queue.shift()
            queue_set.delete(point)
            
            let neighbors = this.neighbors(point)
                                // this ignores contradictions!
                                .filter(n => sum(this.wave[n.point]) > 1)
            
            for (let n of neighbors) {
                this.visited += 1
                let changed = this.propagate(point, n.point, n.offset)
                if (changed > 0 && !queue_set.has(n.point)) {
                    queue.insertBack(n.point)
                    queue_set.add(n.point)
                }
                if (sum(this.wave[n.point]) == 0) {
                    let message = 'contradiction in propagate_from at ' + [point, n.point]
                    console.log(message)
                    throw new Error('contradiction')
                }
            }

        }
    }

    propagate(point1, point2, offset) {
        let states1 = where(this.wave[point1], x => x == 1)
        let states2 = where(this.wave[point2], x => x == 1)
        let table_entry = this.table[offset]
        var changed = 0
        // core loop, optimize later
        for (let s2 of states2) {
            let row = table_entry[s2]
            var flag = false
            // nicer, but some slower than loop on chrome
            // let flag = states1.some(s1 => row[s1])
            // s2 is not supported by any of states1
            for (let s1 of states1) {
                if (row[s1]) {
                    flag = true
                    break
                }
            }
            // s2 is not supported by any of states1
            if (!flag) {
                this.wave[point2][s2] = 0
                changed += 1
            }
        }
        if (changed > 0) {
            this.entropy[point2] = -1 // flag to recalculate
        }
        return changed
    }

    firstCollapse() {
        // mark bad points: sum of 1 (collapsed) or 0 (contradiction)
        let bignum = 9999 
        let sums = this.wave.map(sum)
                            .map(x => x > 1 ? x : bignum)
        let point = argmin(sums)
        if (sums[point] == bignum) {
            return null
        }
        let state = this.collapse(point)
        return {point: point, state: state}
    }

    randomCollapse() {
        // select among uncollapsed points 
        let candidates = filterMap(this.wave, (xs, i) => sum(Array.from(xs)) > 1 ? i : null)
        if (candidates.length == 0) {return null}
        let point = choice(candidates)
        let state = this.collapse(point)
        return {point: point, state: state}
    }

    collapse(point) {
        // skip entropy for now, just pick first one that works
        // let state = this.wave[point].findIndex(x => x > 0)
        let candidates = filterMap(this.wave[point], (x, i) => x > 0 ? i : null)
        let state = choice(candidates)
        this.wave[point].fill(0)
        this.wave[point][state] = 1
        return state
    }

    contradiction() {
        return this.wave.some(xs => sum(xs) == 0)
    }

    done() {
        return this.wave.every(xs => sum(xs) == 1)
    }

    // making this static requires using SimpleTiledModel.f() instead of this.f()
    getNeighbors(FMX, FMY, offsets) {
        // returns an actual getNeighbors function
        // use integers to index offsets (label equivalence classes)
        // remove filter and map with mod to get periodic boundaries

        let neighbors = (a, b) => (Math.abs(a.i-b.i) < 2) && (Math.abs(a.j-b.j) < 2)

        let fixed_boundaries = n => neighbors(toIJ(n.point, FMX), 
                                              toIJ(n.point - offsets[n.offset], FMX))

        let periodic_boundaries = n => {

            let {i:i0,j:j0} = toIJ(n.origin, FMX)
            let {i:i1,j:j1} = offset_toIJ(offsets[n.offset], FMX)
            // console.log(i0,j0,i1,j1)
            let answer = mod(i0 + i1, FMY) * FMX + mod(j0 + j1, FMX)
            let x = {point: answer, offset: n.offset}
            // console.log(n, x)
            return {point: answer, offset: n.offset}
        }

        return point => offsets.map((x, i) => ({point: mod(point + x, FMX*FMY), offset: i, origin: point}))
                               // .filter(fixed_boundaries)
                               .map(periodic_boundaries)
    }

    makeCompatibilityTable(M, offsets, FMX){
        // TypedArray of compatibilities
        // could require a function ((s1,s2,offset)=>Bool
        // TypedArray.map makes another TypedArray of same type
        // make the inner array state1 so it's fast to check support for each s2
        let states1 = rangeU16(M)
        let states2 = range(M)
        return  offsets.map(offset => 
                states2.map(s2 => 
                states1.map(s1 => this.compatibilityTable(s1, s2, offset, FMX))))
    }

    compatibilityTable(s1, s2, offset, FMX) {
        let {i, j} = offset_toIJ(offset, FMX)
        let tile1 = this.initial[s1]
        let tile2 = this.initial[s2]
        return compatibleTiles(tile1, tile2, i, j)
    }

    compatibilityTableDiagonal(s1, s2, offset, FMX) {
        // function to compute compatibility
        let alternate = (s1 == 0 && s2 == 1) || (s1 == 1 && s2 == 0)
        let same = (s1 == s2) && (s1 == 0 || s1 == 1)
        let a = offset_toIJ(offset, FMX)
        let diagonal = (Math.abs(a.i) + Math.abs(a.j)) == 2
        
        if (alternate) {
            // allow if offset is not diagonal
            return !diagonal
        }
        if (same) {
            // allow if offset is diagonal
            return diagonal
        }
        return false
    }

    superposition() {
        var image = Array.from(this.wave).fill(0)
        var alpha = Array.from(this.wave).fill(0)
        
        let toIJ = offset => offset_toIJ(offset, this.FMX)
        let grab = (tile, ij) => tile[ij.i + 1][ij.j + 1]
        let update = (offset, point, iS) => {
            let tile = this.initial[iS]
            let pixel = grab(tile, toIJ(offset))
            image[point] += pixel
            alpha[point] += 1
            }
        
        this.wave.forEach((states, iW) =>
                    states.forEach((state, iS) => state == 0 ? 0 :  
                    this.neighbors(iW).forEach(x =>
                    update(this.offsets[x.offset], x.point, iS))))
        
        return image.map((px, i) => px / alpha[i])        
    }

    to_2D_wavesum () {
        // binary encoding
        let wavesum = this.wave.map(xs => Array.from(xs).map((x, i) => x * Math.pow(2,i)))
                               .map(sum)
        return this.to_2D(wavesum)
    }

    to_2D(arr_1D) {
        let arr = []
        for (let i = 0; i < this.FMY; i++) {
            arr.push(arr_1D.slice(i*this.FMX, (i+1)*this.FMX))
        }
        return arr
    }
}


function compatibleTiles(tile1, tile2, i, j) {
    return _.isEqual(subsquare(tile1, i, j), subsquare(tile2, -i, -j))
}

function subsquare(tile, i, j) {
    var tile = tile
    if (i < 0) {
        tile = tile.slice(0,i)
    } else {
        tile = tile.slice(i)
    }
    if (j < 0) {
        tile = tile.map(xs => xs.slice(0,j))
    } else {
        tile = tile.map(xs => xs.slice(j))
    }
    return tile
}

function filterMap(xs, f) {
    return Array.from(xs).map(f).filter(x => !(x == null))
}

let range_ = arr => arr.fill().map((x,i) => i)

let range = n => range_(Array(n))

// TypedArray.map coerces results into the same type, so
// rangeU8(257)[256] == 0
let rangeU8 = n => range_(new Uint8Array(n))

let rangeU16 = n => range_(new Uint16Array(n))

let watch = x => {console.log(x); return x}

let toIJ = (point, X) => ({i: mod(Math.floor(point / X), X), 
                           j: mod(point, X)})

let offset_toIJ = (offset, X) => {
    var {i, j} = toIJ(offset + X + 1, X)
    return {i: i - 1, j: j - 1}
}

let choice = xs =>
    xs[Math.floor(getRandom()*xs.length)]

let getRandom = getRandomStore()
// let getRandom = Math.random

function getRandomStore() {
    var i = 0
    return () => {i += 1; return randos[i % randos.length]}

}

let mod = (x,y) => x - y * Math.floor(x / y)

let sum = xs => xs.reduce((a,b) => a + b)

let where = (xs, test) => {
    // for loop in case array.map is weird
    indices = []
    for (i = 0; i < xs.length; i++) {
        if (test(xs[i])) {indices.push(i)}
    }
    return indices
}

function superCountBy(acc, x, i) {
    // usage: stuff.reduce(superCountBy, [])
    //
    // employs _.isEqual rather than String or whatever
    // to do the counting. results returned as an array of
    // [value, [keys]]. all this could be averted if 
    // javascript had tuples.
    for (let j = 0; j < acc.length; j++) {
        if (_.isEqual(x, acc[j][0])) {
            acc[j][1].push(i)
            return acc
        }
    }
    // no match, add to acc
    acc.push([x, [i]])
    return acc
}

function argmin(xs) {
    let n = xs.length
    var dmin = Infinity
    var imin = -1;

    for (let i = 0; i < n; i++) {
        if (xs[i] < dmin) {
            dmin = xs[i]
            imin = i
        }
    }

    return imin
}

exports.mod = mod
exports.range = range
exports.toIJ = toIJ
exports.offset_toIJ = offset_toIJ
exports.sum = sum
exports.subsquare = subsquare
exports.compatibleTiles = compatibleTiles
exports.superCountBy = superCountBy
exports.watch = watch
exports.getRandom = getRandom


// function makeOffsets(N) {
//  let offsets = []
//  for (var i = 0; i < N; i++) {
//      i -= N / 2
//      for (var j = 0; j < N; j++) {
//          j -= N / 2
//          if (i != 0 || j != 0) {
//              offsets.push([i, j])
//          }
//      }

//  }
//  return offsets
// }


let checkerboardA = 
    [ [ 0, 1, 0 ]
    , [ 1, 0, 1 ]
    , [ 0, 1, 0 ]
    ]

let checkerboardB = 
    [ [ 1, 0, 1 ]
    , [ 0, 1, 0 ]
    , [ 1, 0, 1 ]
    ]

let checkerboardC = 
    [ [ 1, 0, 0 ]
    , [ 0, 1, 1 ]
    , [ 1, 0, 0 ]
    ]

let checkerboardD = 
    [ [ 0, 0, 1 ]
    , [ 1, 1, 0 ]
    , [ 0, 0, 1 ]
    ]

let checkerboardE = 
    [ [ 0, 1, 1 ]
    , [ 1, 0, 0 ]
    , [ 0, 1, 1 ]
    ]

let checkerboardF = 
    [ [ 1, 1, 0 ]
    , [ 0, 0, 1 ]
    , [ 1, 1, 0 ]
    ]

let checkerboardABCDEF = 
    [ checkerboardA
    , checkerboardB
    , checkerboardC
    , checkerboardD
    , checkerboardE
    , checkerboardF
    ]

exports.checkerboardABCDEF = checkerboardABCDEF

let randos = 
      [  3.57095269e-01,   6.16476608e-01,   4.42510505e-01,
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

let knot = 
    [[[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 0, 0], [1, 0, 1]],
     [[1, 1, 1], [0, 0, 0], [0, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 0]],
     [[1, 1, 1], [0, 0, 1], [1, 0, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
     [[1, 0, 0], [1, 0, 1], [1, 0, 1]],
     [[0, 0, 0], [0, 1, 1], [0, 1, 0]],
     [[0, 0, 0], [1, 1, 1], [1, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 1]],
     [[0, 0, 0], [1, 1, 0], [0, 1, 0]],
     [[0, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 1], [0, 1, 0], [0, 1, 0]],
     [[1, 1, 1], [1, 0, 0], [1, 0, 1]],
     [[1, 1, 1], [0, 0, 0], [0, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 0]],
     [[1, 1, 1], [0, 0, 1], [1, 0, 1]],
     [[1, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 0], [1, 0, 1], [1, 0, 1]],
     [[0, 0, 0], [0, 1, 1], [0, 1, 1]],
     [[0, 0, 0], [1, 1, 1], [1, 1, 1]],
     [[0, 0, 0], [1, 1, 0], [1, 1, 0]],
     [[0, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
     [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
     [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
     [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
     [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
     [[1, 0, 1], [1, 0, 1], [0, 0, 0]],
     [[0, 1, 0], [0, 1, 0], [0, 0, 0]],
     [[1, 0, 1], [1, 0, 1], [0, 0, 0]],
     [[0, 1, 1], [0, 1, 1], [0, 0, 0]],
     [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
     [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
     [[1, 0, 1], [1, 0, 1], [0, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[1, 1, 1], [1, 0, 0], [1, 0, 1]],
     [[1, 1, 1], [0, 0, 0], [0, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 0], [0, 0, 0], [1, 1, 1]],
     [[1, 0, 1], [0, 0, 0], [1, 1, 1]],
     [[0, 1, 0], [0, 0, 0], [1, 1, 1]],
     [[1, 0, 1], [0, 0, 0], [1, 1, 1]],
     [[0, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 0], [0, 0, 0], [1, 1, 1]],
     [[1, 0, 1], [0, 0, 1], [1, 1, 1]],
     [[0, 1, 0], [0, 1, 0], [1, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[1, 0, 0], [1, 0, 1], [1, 0, 1]],
     [[0, 0, 0], [0, 1, 1], [0, 1, 0]],
     [[0, 0, 0], [1, 1, 1], [1, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 1], [1, 1, 1], [0, 0, 0]],
     [[0, 1, 0], [1, 1, 0], [0, 0, 0]],
     [[1, 0, 1], [1, 0, 1], [0, 0, 1]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 1], [0, 1, 0], [0, 1, 0]],
     [[1, 1, 1], [1, 0, 0], [1, 0, 1]],
     [[1, 1, 1], [0, 0, 0], [0, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 0]],
     [[1, 1, 1], [0, 0, 0], [1, 0, 1]],
     [[1, 1, 1], [0, 0, 0], [0, 1, 0]],
     [[1, 1, 1], [0, 0, 0], [1, 0, 1]],
     [[1, 1, 1], [0, 0, 0], [0, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 0], [0, 0, 0], [1, 1, 1]],
     [[1, 0, 1], [0, 0, 1], [1, 1, 1]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 0], [1, 0, 1], [1, 0, 1]],
     [[0, 0, 0], [0, 1, 1], [0, 1, 1]],
     [[0, 0, 0], [1, 1, 1], [1, 1, 1]],
     [[0, 0, 0], [1, 1, 0], [1, 1, 0]],
     [[0, 0, 0], [1, 0, 1], [1, 0, 1]],
     [[0, 0, 0], [0, 1, 0], [0, 1, 0]],
     [[0, 0, 0], [1, 0, 1], [1, 0, 1]],
     [[0, 0, 0], [0, 1, 1], [0, 1, 1]],
     [[0, 0, 0], [1, 1, 1], [1, 1, 1]],
     [[0, 0, 0], [1, 1, 1], [1, 1, 1]],
     [[0, 0, 0], [1, 1, 1], [1, 1, 1]],
     [[0, 0, 0], [1, 1, 1], [1, 1, 1]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 0]],
     [[0, 1, 1], [0, 1, 1], [0, 0, 0]],
     [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
     [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
     [[1, 0, 1], [1, 0, 1], [0, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 0], [0, 1, 0], [0, 1, 1]],
     [[1, 0, 1], [1, 0, 0], [1, 1, 1]],
     [[0, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 0], [0, 0, 0], [1, 1, 1]],
     [[1, 0, 1], [0, 0, 1], [1, 1, 1]],
     [[0, 1, 0], [0, 1, 0], [1, 1, 0]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 1]],
     [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 0, 1], [1, 0, 1], [1, 0, 0]],
     [[0, 1, 0], [0, 1, 1], [0, 0, 0]],
     [[1, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
     [[0, 0, 1], [1, 1, 1], [0, 0, 0]],
     [[0, 1, 0], [1, 1, 0], [0, 0, 0]],
     [[1, 0, 1], [1, 0, 1], [0, 0, 1]],
     [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 0, 1], [1, 0, 0], [1, 1, 1]],
     [[0, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 1], [0, 0, 0], [1, 1, 1]],
     [[1, 1, 0], [0, 0, 0], [1, 1, 1]],
     [[1, 0, 1], [0, 0, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]

exports.knot = knot

let dungeon = 
[[[0, 0, 0], [0, 0, 0], [0, 0, 1]],
 [[0, 0, 0], [0, 0, 0], [0, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 0], [0, 1, 1], [0, 1, 0]],
 [[0, 0, 0], [1, 1, 1], [1, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 1]],
 [[0, 0, 0], [1, 1, 1], [0, 1, 0]],
 [[0, 0, 0], [1, 1, 1], [1, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 1]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 1, 1], [0, 1, 0], [0, 1, 0]],
 [[1, 1, 1], [1, 0, 0], [1, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 1], [0, 0, 1]],
 [[1, 1, 1], [0, 1, 0], [0, 1, 0]],
 [[1, 1, 1], [1, 0, 0], [1, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
 [[1, 0, 0], [1, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
 [[1, 0, 0], [1, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 1], [0, 0, 1], [1, 1, 1]],
 [[0, 1, 0], [0, 1, 0], [1, 1, 0]],
 [[1, 0, 0], [1, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 1, 0], [0, 1, 0], [0, 1, 0]],
 [[1, 0, 0], [1, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 1], [1, 1, 1], [0, 0, 0]],
 [[0, 1, 0], [1, 1, 0], [0, 0, 0]],
 [[1, 0, 0], [1, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 1]],
 [[0, 0, 1], [0, 0, 1], [0, 1, 1]],
 [[0, 1, 0], [0, 1, 0], [1, 1, 1]],
 [[1, 0, 0], [1, 0, 0], [1, 1, 0]],
 [[0, 0, 0], [0, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 0], [0, 0, 0], [0, 0, 0]],
 [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 1], [0, 1, 1], [0, 1, 1]],
 [[0, 1, 0], [1, 1, 1], [1, 1, 1]],
 [[1, 0, 0], [1, 1, 0], [1, 1, 0]],
 [[0, 0, 0], [1, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
 [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
 [[1, 0, 0], [1, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 1]],
 [[0, 0, 0], [0, 0, 0], [0, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 1], [0, 0, 1], [1, 1, 1]],
 [[0, 1, 1], [0, 1, 1], [1, 1, 1]],
 [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 0], [1, 1, 0], [1, 1, 0]],
 [[1, 0, 0], [1, 0, 0], [1, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 1], [0, 0, 1]],
 [[0, 0, 0], [0, 1, 1], [0, 1, 1]],
 [[0, 0, 0], [1, 1, 1], [1, 1, 1]],
 [[0, 0, 1], [1, 1, 1], [1, 1, 1]],
 [[0, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 0], [1, 1, 0], [1, 1, 1]],
 [[1, 0, 0], [1, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [1, 1, 1]],
 [[0, 0, 1], [0, 0, 1], [1, 1, 1]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
 [[0, 1, 1], [0, 1, 1], [0, 1, 1]],
 [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
 [[1, 1, 0], [1, 1, 1], [1, 1, 0]],
 [[1, 0, 0], [1, 1, 1], [1, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 1], [1, 1, 1], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 1], [0, 0, 0]],
 [[0, 1, 1], [0, 1, 1], [0, 0, 0]],
 [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
 [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
 [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
 [[1, 1, 1], [1, 1, 1], [0, 0, 0]],
 [[1, 1, 1], [1, 1, 0], [0, 0, 0]],
 [[1, 1, 1], [1, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 1], [0, 0, 0], [0, 0, 0]],
 [[0, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 1], [0, 0, 0], [0, 0, 0]],
 [[1, 1, 0], [0, 0, 0], [0, 0, 0]],
 [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
 [[0, 0, 0], [0, 0, 0], [0, 0, 0]]]

 exports.dungeon = dungeon