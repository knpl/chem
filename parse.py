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
            token = Token('END', '', -1, -1)

        return token

    def peek(self):
        if self._peek:
            return self._peek # Keep returning peek until it is consumed by __next__

        try:
            mob = next(self._mobs)
            group = mob.lastgroup
            self._peek = Token(group, mob.group(group), *mob.span())
        except StopIteration:
            self._peek = Token('END', '', -1, -1)

        return self._peek

class Parser:
    def __init__(self, pattern):
        self._pattern = pattern
        self._tokens = Tokenizer(pattern)
        # initialize left-hand-side and right-hand-side molecules
        self._lhs = []
        self._rhs = []
        self._use_rhs = False

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
            
            l.append(self.molecule())
        
            if self.test('SPACE'):
                self.token('SPACE')
            
            if self.test('PLUS'):
                self.token('PLUS')
            else:
                return True

    def molecule(self):
        mol = {}
        while True:
            if self.test('LPAREN'):            
                self.token('LPAREN')
                submol = self.molecule()
                self.token('RPAREN')

                count = 1
                if self.test('NUM'):
                    count = int(self.token('NUM').pattern())

                muladd(mol, submol, count)

            elif self.test('ELEM'):
                elem = self.token('ELEM').pattern()
                count = 1
                if self.test('NUM'):
                    count = int(self.token('NUM').pattern())

                mol[elem] = count + (mol[elem] if elem in mol else 0)
            else:
                break

        
        if len(mol) == 0:
            tok = self._tokens.peek()
            raise RuntimeError(
                'Parse Error: "{}" (type {}) at {} (expected ELEM or LPAREN token).' \
                .format(tok.pattern(), tok.type(), tok.start()))

        return mol
    
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

def solve(lhs, rhs):
    elems = set()
    for mol in lhs:
        for elem in mol.keys():
            elems.add(elem)
    for mol in rhs:
        for elem in mol.keys():
            elems.add(elem)
    index = list(sorted(elems))
    revindex = {elem:idx for idx,elem in enumerate(index)}
    n = len(index)

    print('Element index:')
    print(index)
    print('Reverse index:')
    print(revindex)

    system = np.zeros((len(index), len(lhs) + len(rhs)), dtype=np.int32)
    for col, mol in enumerate(lhs):
        for elem, count in mol.items():
            system[revindex[elem], col] = count
    for col, mol in enumerate(rhs, len(lhs)):
        for elem, count in mol.items():
            system[revindex[elem], col] = -count
    print('System of equations:')
    print(system)
    print('Row-reduced:')
    rowreduce(system)
    print(system)

def gcd(a,b):
    (a, b) = (abs(a), abs(b))
    if a == 0:
        return b
    while b != 0:
        (a, b) = (b, a % b)
    return a

def rowreduce(system):
    # n rows, m columns
    (n, m) = system.shape
    col = 0
    pivots = 0
    pcols = []
    while pivots < n and col < m:
        # Find pivot in current column
        pivot = -1
        for k in range(pivots, n):
            if system[k, col] != 0:
                pivot = k
                break
        if pivot == -1:
            col += 1
            continue

        # Found pivot. Clear remaining column
        pcols.append(col)
        prow = system[pivot]
        pval = prow[col]
        for k in chain(range(pivots), range(pivot+1, n)):
            row = system[k]
            val = row[col]
            if val != 0:
                g = gcd(val, pval)
                row *= pval//g
                row[col:] -= (val//g)*prow[col:]

        # Swap rows if necessary
        if pivot != pivots:
            system[[pivots, pivot]] = system[[pivot, pivots]]

        pivots += 1
        col += 1

    for i, j in enumerate(pcols):
        row = system[i][j:]
        row //= reduce(gcd, row)
        if row[0] < 0:
            system[i][j:] = -row
        
if __name__ == '__main__':
    p = Parser('C3H8 + O2 -> H2O + CO2')
    p.parse()
    lhs = p.get_lhs()
    rhs = p.get_rhs()
    print('Left hand side:')
    for d in lhs:
        print(d)
    print('\nRight hand side:')
    for d in rhs:
        print(d)

    solve(lhs, rhs)
            
