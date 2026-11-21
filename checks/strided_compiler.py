"""Specialize the existing quadratic graph to public fixed bits on all but one slot."""
from copy import deepcopy

from c38_network import Compiled, Pinned, arith, sc


def free_positions(b,k,owners,ell):
    """Few owners, and a single contiguous free run per group of owner limbs."""
    if k<1 or owners<1:
        raise ValueError("A limb and an owner group must both have positive width.")
    if not 0 <= ell <= b:
        raise ValueError("The width exponent must lie within the abscissa word.")
    ordered=sorted(range(b),key=lambda bit:((bit//k)%owners,bit))
    return tuple(sorted(ordered[:ell]))


class FixedBits:
    def __init__(self,base,patterns,free_mask):
        self.base,self.patterns,self.free_mask=base,patterns,free_mask
        self.uid,self.ports=base.uid,base.ports
        self.alphabet,self.pc,self.axes=base.alphabet,base.pc,base.axes

    def entry(self,symbols,ring):
        value=self.base.entry(symbols,ring)
        if value==ring.zero:
            return value
        assignment=self.base.local_assignment(symbols)
        par=self.base.data["par"]
        row=self.base.index
        for name,pattern in self.patterns.items():
            words=[arith.read_word(w,par) for w in assignment["hs"][name][row]]
            for j,limb in enumerate(arith.H_inverse_row(words,row,par)):
                shift=par["k"]*(j*par["L2"]+row)
                fixed=(~self.free_mask >> shift) & (par["B"]-1)
                if ((limb ^ (pattern >> shift)) & fixed)!=0:
                    return ring.zero
        return value


def specialize(data,dimension,positions,patterns):
    """Build from public seed patterns. Slot x1 keeps its ordinary affine hash."""
    data=deepcopy(data)
    par=data["par"]
    if set(patterns)!={f"x{i}" for i in range(2,data["m"]+1)}:
        raise ValueError("Expected one public pattern per fixed-bit slot.")
    positions=tuple(positions)
    if len(set(positions))!=len(positions) or any(i<0 or i>=par["b"] for i in positions):
        raise ValueError("Invalid free bit position.")
    free_mask=sum(1<<i for i in positions)
    if any(c<0 or c>=1<<par["b"] or c&free_mask for c in patterns.values()):
        raise ValueError("Patterns must specify just the fixed abscissa bits.")
    for i in range(2,data["m"]+1):
        # The old hash-only expressions now have identically zero operands.
        data["constants"][f"alpha{i}"]=data["constants"][f"beta{i}"]=0
    base=Compiled(data,dimension)
    active={(i//par["k"])%par["L2"] for i in positions}
    H={name:arith.H(arith.limbs(c,par),par) for name,c in patterns.items()}
    zero_carry=arith.signed_word(0,par["carry_bound"],par["carry_bits"])
    pins={}
    first_unused=arith.arithmetic_count(data["m"])+1
    for edge,path in base.paths.items():
        if path[0]=="outs" and path[1]>=first_unused:
            pins[edge]=0
        elif path[0]=="carry_ports" and path[1]>=first_unused:
            pins[edge]=zero_carry[path[-1]]
        elif path[0]=="hs":
            _,name,row,freq,digit=path
            if name in patterns and row not in active:
                pins[edge]=(H[name][row][freq]//base.codec.Q**digit)%base.codec.Q
    rows,cols=[],[]
    for group,out in ((base.net.rows,rows),(base.net.cols,cols)):
        for tensor in group:
            source=(FixedBits(tensor,patterns,free_mask)
                    if getattr(tensor,"kind",None)=="first" else tensor)
            restricted=Pinned(source,pins)
            if restricted.ports:
                out.append(restricted)
            elif restricted.entry((),sc.Ring(1))!=sc.Ring(1).one:
                raise ValueError("An isolated public predicate rejects.")
    remaining=tuple(e for e in base.net.edges if e not in pins)
    net=sc.Network("C38-fixed-coordinate-bits",rows,cols,remaining,base.codec.alphabet)
    return base,net,pins
