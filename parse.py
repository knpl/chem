import re
import argparse

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
        while True:

            if self.test('SPACE'):
                self.token('SPACE')
            
            self.molecule()
        
            if self.test('SPACE'):
                self.token('SPACE')
            
            if self.test('PLUS'):
                self.token('PLUS')

            else:
                return True

    def molecule(self, mol={}):
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

        
if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('pattern', metavar='PATTERN', type=str, help='The pattern.')
    args = argparser.parse_args()

    if args.pattern:
        p = Parser(args.pattern)
        p.parse()
        print('Success')
#        for token in Tokenizer(args.pattern):
#            if token.type() == 'INVALID':
#                raise RuntimeError('Invalid token: "{}"\n{}\n{}^' \
#                                   .format(token.pattern(), args.pattern, ' '*token.start()))
#            print('{}: "{}"'.format(token.type(), token.pattern()))
            
