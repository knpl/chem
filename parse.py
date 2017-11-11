import re
from functools import reduce
from itertools import chain
import numpy as np

class Token:
    def __init__(self, tokentype, pattern, start, end):
        self._type = tokentype
        self._pattern = pattern
        self._start = start
        self._end = end

    def type(self):
        return self._type

    def pattern(self):
        return self._pattern

    def start(self):
        return self._start

    def end(self):
        return self._end

    def span(self):
        return (self._start, self._end)

class Tokenizer:
    tokens = [
        ('ELEM', r'[A-Z][a-z]{0,2}'),
        ('NUM', r'[1-9][0-9]*'),
        ('LPAREN', r'\('),
        ('RPAREN', r'\)'),
        ('PLUS', r'\+'),
        ('ARROW', r'->'),
        ('SPACE', r'[ \t]+'),
        ('INVALID', r'.')]

    matcher = re.compile('|'.join('(?P<{}>{})'. \
                                  format(*pair) for pair in tokens))

    def __init__(self, pattern):
        self._mobs = Tokenizer.matcher.finditer(pattern)
        l = len(pattern)
        self._end_token = Token('END', '', l, l)
        self._peek = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._peek:
            token = self._peek
            self._peek = None # Clear peek
            return token

        token = None
        try:
            mob = next(self._mobs)
            group = mob.lastgroup
            token = Token(group, mob.group(group), *mob.span())
        except StopIteration:
            token = self._end_token

        return token

    def peek(self):
        if self._peek:
            return self._peek # Keep returning peek until it is consumed by __next__

        try:
            mob = next(self._mobs)
            group = mob.lastgroup
            self._peek = Token(group, mob.group(group), *mob.span())
        except StopIteration:
            self._peek = self._end_token

        return self._peek

class Parser:
    def __init__(self, pattern):
        self._pattern = pattern
        self._tokens = Tokenizer(pattern)
        # initialize left-hand-side and right-hand-side molecules
        self._lhs = []
        self._rhs = []
        self._use_rhs = False
        self._mol = None
        self._molspan = None

    def get_lhs(self):
        return self._lhs

    def get_rhs(self):
        return self._rhs

    def parse(self):
        self.formula()

    def formula(self):
        # Parse left hand side
        self.molecules()
        self.token('ARROW')
        # Parse right hand side
        self._use_rhs = True
        self.molecules()
        # Parse end
        self.token('END')

    def molecules(self):
        l = self._rhs if self._use_rhs else self._lhs
        while True:

            if self.test('SPACE'):
                self.token('SPACE')
            
            self.molecule()
            (start, end) = self._molspan
            l.append((self._pattern[start:end], self._mol))
        
            if self.test('SPACE'):
                self.token('SPACE')
            
            if self.test('PLUS'):
                self.token('PLUS')
            else:
                return True

    def molecule(self):
        mol = {}
        start = self._tokens.peek().start()
        end = start
        while True:
            if self.test('LPAREN'):            
                self.token('LPAREN')
                self.molecule()
                submol = self._mol
                end = self.token('RPAREN').end()

                count = 1
                if self.test('NUM'):
                    tok = self.token('NUM')
                    count = int(tok.pattern())
                    end = tok.end()

                muladd(mol, submol, count)

            elif self.test('ELEM'):
                tok = self.token('ELEM')
                elem = tok.pattern()
                end = tok.end()

                count = 1
                if self.test('NUM'):
                    tok = self.token('NUM')
                    count = int(tok.pattern())
                    end = tok.end()

                mol[elem] = count + (mol[elem] if elem in mol else 0)
            else:
                break

        
        if start == end:
            tok = self._tokens.peek()
            raise RuntimeError(
                'Parse Error: "{}" (type {}) at {} (expected ELEM or LPAREN token).' \
                .format(tok.pattern(), tok.type(), tok.start()))
        
        self._mol = mol
        self._molspan = (start, end)
    
    def token(self, tokentype):
        token = next(self._tokens)
        if token.type() != tokentype:
            raise RuntimeError('Parse Error: "{}" (type: {}) at {} (expected: {}).' \
                               .format(token.pattern(), token.type(), token.start(), tokentype))
        return token

    def test(self, tokentype):
        return self._tokens.peek().type() == tokentype

    def muladd(mola, molb, n):
        for elem, count in molb.items():
            mola[elem] = n * count + (mola[elem] if elem in mola else 0)

def lcm(a, b):
    if a == 0 or b == 0:
        return 0
    (a, b) = (abs(a), abs(b))
    (p, q) = (a, b)
    while q != 0:
        (p, q) = (q, p % q)
    return (a // p) * b

def gcd(a, b):
    (a, b) = (abs(a), abs(b))
    if a == 0:
        return b
    while b != 0:
        (a, b) = (b, a % b)
    return a

# Compute row-reduced system of equations and compute its nullspace.
# Returns base of nullspace.
def nullspace(system):
    (m, n) = system.shape
    pivot_cols = rowreduce(system)
    pivot_vals = [system[i, j] for i, j in enumerate(pivot_cols)]

    nullbase = np.zeros((n - len(pivot_cols), n), dtype=np.int32)

    plcm = 1
    pcount = 0
    for j in range(n):
        if pcount < len(pivot_cols) and j == pivot_cols[pcount]:
            # j-th column is pivot
            plcm = lcm(plcm, pivot_vals[pcount])
            pcount += 1            
        else: 
            # j-th column is no pivot
            solution = nullbase[j - pcount]

            pgcd = plcm
            for k, v in enumerate(system[:pcount, j]):
                (pcol, pval) = (pivot_cols[k], pivot_vals[k])
                entry = -(plcm//pval)*v
                solution[pcol] = entry
                pgcd = gcd(pgcd, entry)

            for pcol in pivot_cols:
                solution[pcol] //= pgcd
            solution[j] = plcm // pgcd

    return nullbase

# Compute row-reduced system of equations. 
# Returns list containing pivot column indices.
def rowreduce(system):
    # m rows, n columns
    (m, n) = system.shape
    
    pcount = 0
    pcols = []

    j = 0
    while pcount < m and j < n:
        # Scan j-th column for pivot.
        pidx = -1
        for i in range(pcount, m):
            if system[i, j] != 0:
                pidx = i
                break

        if pidx < 0:
            # j-th column is no pivot column.
            j += 1
            continue

        # j-th column is pivot column. Clear pivot column.
        pcols.append(j)

        prow = system[pidx]
        pval = prow[j]
        for i in chain(range(pcount), range(pidx+1, m)):
            row = system[i]
            val = row[j]
            if val != 0:
                g = gcd(val, pval)
                row *= pval//g
                row[j:] -= (val//g)*prow[j:]

        # Swap rows if necessary
        if pidx != pcount:
            system[[pcount, pidx]] = system[[pidx, pcount]]

        pcount += 1
        j += 1

    # Make pivot values positive and rows relprime.
    for i, j in enumerate(pcols):
        row = system[i][j:]
        g = reduce(gcd, row)
        row //= g if row[0] >= 0 else -g

    return pcols

def solve(lhs, rhs):
    elems = set()
    for mol in lhs:
        for elem in mol.keys():
            elems.add(elem)
    for mol in rhs:
        for elem in mol.keys():
            elems.add(elem)
    index = list(sorted(elems))
    revindex = {elem: idx for idx, elem in enumerate(index)}

    system = np.zeros((len(index), len(lhs) + len(rhs)), dtype=np.int32)
    for col, mol in enumerate(lhs):
        for elem, count in mol.items():
            system[revindex[elem], col] = count
    for col, mol in enumerate(rhs, len(lhs)):
        for elem, count in mol.items():
            system[revindex[elem], col] = -count

    return nullspace(system)

if __name__ == '__main__':
    formulas = ['C3H8 + O2 -> H2O + CO2']
    for formula in formulas:
        print('{}'.format(formula))
        p = Parser(formula)

        p.parse()
        lhs = p.get_lhs()
        rhs = p.get_rhs()
        (lhs_formulas, lhs_molecules) = zip(*lhs)
        (rhs_formulas, rhs_molecules) = zip(*rhs)

        nullbase = solve(lhs_molecules, rhs_molecules)
        for row in nullbase:
            lhs_str = ' + '.join('{}{}'.format(n, s) for n, s in zip(row[:len(lhs)], lhs_formulas))
            rhs_str = ' + '.join('{}{}'.format(n, s) for n, s in zip(row[len(lhs):], rhs_formulas))
            print(lhs_str + ' -> ' + rhs_str)
