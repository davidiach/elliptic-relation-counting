"""Decision from weighted Gram moments on sign-paired slots of odd-order curves."""
from flint import nmod_poly
from itertools import combinations


def power_sums(xs, weights, length, p):
    sums = [0]*length
    for x,weight in zip(xs,weights,strict=True):
        term = weight
        for j in range(length):
            sums[j] += term
            term = term*x % p
    return [entry % p for entry in sums]


def moment_channels(slots, tags, p):
    """Each channel weights one slot by x^j and the other two by their tags."""
    tagged = [power_sums(xs,weights,5,p) for xs,weights in zip(slots,tags,strict=True)]
    channels = []
    for i,xs in enumerate(slots):
        quadratics = []
        for h in ((i+1)%3,(i+2)%3):
            s = tagged[h]
            quadratics.append([nmod_poly([s[r+2],-2*s[r+1],s[r]],p) for r in range(3)])
        a,b = quadratics
        weight = a[2]*b[0]-2*a[1]*b[1]+a[0]*b[2]
        sums = power_sums(xs,[1]*len(xs),2*len(xs)+4,p)
        channels.append([sum(int(weight[r])*sums[j+r] for r in range(5)) % p
                         for j in range(2*len(xs))])
    return channels


def visible_locator(moments,p):
    """Berlekamp--Massey, including a possible support abscissa at zero."""
    connection, previous = [1], [1]
    degree, delay, pivot = 0, 1, 1
    for n, val in enumerate(moments):
        diff = (val+sum(connection[j]*moments[n-j] for j in range(1,degree+1))) % p
        if not diff:
            delay += 1
            continue
        saved = connection.copy()
        if len(connection) < len(previous)+delay:
            connection.extend([0]*(len(previous)+delay-len(connection)))
        scale = diff*pow(pivot,-1,p) % p
        for j,coefficient in enumerate(previous):
            connection[j+delay] = (connection[j+delay]-scale*coefficient) % p
        if 2*degree <= n:
            degree,previous,pivot,delay = n+1-degree,saved,diff,1
        else:
            delay += 1
    return nmod_poly(connection[:degree+1][::-1],p)


def complete_locator(visible,size,p,repair_sums=None):
    degree = visible.degree()
    missing = size-degree
    if missing == 0:
        return visible
    if repair_sums is None or missing not in (1,2):
        return None
    # Newton sums for the visible roots. Leading coefficient is one.
    s1 = -int(visible[degree-1]) if degree else 0
    s2 = s1*s1 - (2*int(visible[degree-2]) if degree>=2 else 0)
    u = (repair_sums[0]-s1) % p
    v = (repair_sums[1]-s2) % p
    if missing==1:
        return visible*nmod_poly([-u,1],p)
    return visible*nmod_poly([(u*u-v)*pow(2,-1,p),-u,1],p)


def translate(locator_poly, t, a, c, p):
    """Product of S_3(w,t,X) over locator roots w, without extracting the roots."""
    X = nmod_poly([0,1],p)
    L = (X-t)**2
    M = -2*((X+t)*(X*t+a)+2*c)
    C = (X*t-a)**2-4*c*(X+t)
    LC = L*C

    # Work with a+bJ where J^2 = -M J-LC. Homogenization avoids inverting L.
    def multiply(left,right):
        u,v = left
        s,r = right
        vr = v*r
        return u*s-vr*LC, u*r+v*s-vr*M

    zero = nmod_poly([],p)
    one = nmod_poly([1],p)
    coeffs = locator_poly.coeffs()
    level = [(nmod_poly([int(v)],p),zero) for v in coeffs]
    width = 1 << (len(level)-1).bit_length()
    level.extend([(zero,zero)]*(width-len(level)))
    j_power, l_power = (zero,one), L
    while len(level)>1:
        merged = []
        for lo,hi in zip(level[::2],level[1::2],strict=True):
            u,v = multiply(hi,j_power)
            merged.append((lo[0]*l_power+u,lo[1]*l_power+v))
        level = merged
        if len(level)>1:
            j_power = multiply(j_power,j_power)
            l_power *= l_power
    u,v = level[0]
    norm = u*u-M*u*v+LC*v*v
    quotient,remainder = divmod(norm,L**(2*(width-1)-locator_poly.degree()))
    if not remainder.is_zero():
        raise ArithmeticError("Homogeneous resultant has a nonzero division remainder")
    return quotient


def decide(sizes, channels, t, a, c, p, repair_sums=None, seen_locators=None):
    """Read only sizes, weighted moments, curve and target. None requests a retry."""
    if 0 in sizes:
        return False
    if repair_sums is None:
        repair_sums = [None]*len(sizes)
    supports = []
    for i,(seq,k,sums) in enumerate(zip(channels,sizes,repair_sums,strict=True)):
        poly=visible_locator(seq,p)
        if seen_locators is not None:
            previous=seen_locators[i]
            poly=(poly//poly.gcd(previous))*previous
            seen_locators[i]=poly
        supports.append(complete_locator(poly,k,p,sums))
    if any(poly is None for poly in supports):
        return None
    return join(supports,t,a,c,p)


def duplicate_probes(u,t,a,c,p):
    """Possible finite third abscissas when two slots each force the same u."""
    f=(u**3+a*u+c)%p
    twice=((3*u*u+a)**2*pow(4*f,-1,p)-2*u)%p
    shifted=translate(nmod_poly([-twice,1],p),t,a,c,p)
    return {t,*(int(z) for z,multiplicity in shifted.roots())}


def decide_slots(slots,t,a,c,p,rng,contains):
    """Complete three-slot decision. contains(i,x) tests the original slot hash."""
    sizes=list(map(len,slots))
    if 0 in sizes:
        return False
    sums=[(sum(s)%p,sum(x*x for x in s)%p) for s in slots]
    for i,j in combinations(range(3),2):
        if sizes[i]==sizes[j]==1 and sums[i][0]==sums[j][0]:
            return any(contains(3-i-j,z) for z in duplicate_probes(sums[i][0],t,a,c,p))
    seen=[nmod_poly([1],p) for s in slots]
    tags=[[1]*k for k in sizes]
    while True:
        result=decide(sizes,moment_channels(slots,tags,p),t,a,c,p,sums,seen)
        if result is not None:
            return result
        tags=[[rng.randrange(1,p) for x in s] for s in slots]


def join(supports,t,a,c,p):
    """Match all finite pair-sum abscissas and the common point at infinity."""
    smallest,middle,largest = sorted(supports,key=lambda poly:poly.degree())
    if smallest(t)==0 and middle.gcd(largest).degree()>0:
        return True
    probes = translate(smallest,t,a,c,p)
    for z,multiplicity in probes.roots():
        translated = translate(middle,int(z),a,c,p)
        if largest.gcd(translated).degree()>0:
            return True
    return False
