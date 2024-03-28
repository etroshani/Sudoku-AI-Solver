## Solve Every Sudoku Puzzle
##Travail complete par:
##Eran TROSHANI 20187276
##et
##Qiwu WEN 20230961

## See http://norvig.com/sudoku.html

## Throughout this program we have:
##   r is a row,    e.g. 'A'
##   c is a column, e.g. '3'
##   s is a square, e.g. 'A3'
##   d is a digit,  e.g. '9'
##   u is a unit,   e.g. ['A1','B1','C1','D1','E1','F1','G1','H1','I1']
##   grid is a grid,e.g. 81 non-blank chars, e.g. starting with '.18...7...
##   values is a dict of possible values, e.g. {'A1':'12349', 'A2':'8', ...}

def cross(A, B):
    "Cross product of elements in A and elements in B."
    return [a + b for a in A for b in B]


digits = '123456789'
rows = 'ABCDEFGHI'
cols = digits
squares = cross(rows, cols)
unitlist = ([cross(rows, c) for c in cols] +
            [cross(r, cols) for r in rows] +
            [cross(rs, cs) for rs in ('ABC', 'DEF', 'GHI') for cs in ('123', '456', '789')])
units = dict((s, [u for u in unitlist if s in u])
             for s in squares)
peers = dict((s, set(sum(units[s], [])) - set([s]))
             for s in squares)


################ Unit Tests ################

def test():
    "A set of tests that must pass."
    assert len(squares) == 81
    assert len(unitlist) == 27
    assert all(len(units[s]) == 3 for s in squares)
    assert all(len(peers[s]) == 20 for s in squares)
    assert units['C2'] == [['A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2'],
                           ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'],
                           ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']]
    assert peers['C2'] == set(['A2', 'B2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2',
                               'C1', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9',
                               'A1', 'A3', 'B1', 'B3'])
    print ('All tests pass.')


################ Parse a Grid ################

def parse_grid(grid):
    """Convert grid to a dict of possible values, {square: digits}, or
    return False if a contradiction is detected."""
    ## To start, every square can be any digit; then assign values from the grid.
    values = dict((s, digits) for s in squares)
    for s, d in grid_values(grid).items():
        if d in digits and not assign(values, s, d):
            return False  ## (Fail if we can't assign d to square s.)
    return values


def grid_values(grid):
    "Convert grid into a dict of {square: char} with '0' or '.' for empties."
    chars = [c for c in grid if c in digits or c in '0.']
    assert len(chars) == 81
    return dict(zip(squares, chars))


################ Constraint Propagation ################

def assign(values, s, d):
    """Eliminate all the other values (except d) from values[s] and propagate.
    Return values, except return False if a contradiction is detected."""
    other_values = values[s].replace(d, '')
    if all(eliminate(values, s, d2) for d2 in other_values):
        return values
    else:
        return False


def eliminate(values, s, d):
    """Eliminate d from values[s]; propagate when values or places <= 2.
    Return values, except return False if a contradiction is detected."""
    if d not in values[s]:
        return values  ## Already eliminated
    values[s] = values[s].replace(d, '')
    ## (1) If a square s is reduced to one value d2, then eliminate d2 from the peers.
    if len(values[s]) == 0:
        return False  ## Contradiction: removed last value
    elif len(values[s]) == 1:
        d2 = values[s]
        if not all(eliminate(values, s2, d2) for s2 in peers[s]):
            return False
    ## (2) If a unit u is reduced to only one place for a value d, then put it there.
    for u in units[s]:
        dplaces = [s for s in u if d in values[s]]
        if len(dplaces) == 0:
            return False  ## Contradiction: no place for this value
        elif len(dplaces) == 1:
            # d can only be in one place in unit; assign it there
            if not assign(values, dplaces[0], d):
                return False
    return values


################ Display as 2-D grid ################

def display(values):
    "Display these values as a 2-D grid."
    width = 1 + max(len(values[s]) for s in squares)
    line = '+'.join(['-' * (width * 3)] * 3)
    for r in rows:
        print ('').join(values[r + c].center(width) + ('|' if c in '36' else '')
                        for c in cols)
        if r in 'CF': print (line)
    print


################ Search ################

def solve(grid): return search(parse_grid(grid))

def worst_solve(grid): return worstSearch(parse_grid(grid))

def random_solve(grid): return randomSearch(parse_grid(grid))

def nakedPair_solve(grid):
    values = naked_pairs(parse_grid(grid))
    return search(values)

def nakedTriple_solve(grid):
    values = naked_triples(parse_grid(grid))
    return search(values)


def search(values):
    "Using depth-first search and propagation, try all possible values."
    if values is False:
        return False  ## Failed earlier
    if all(len(values[s]) == 1 for s in squares):
        return values  ## Solved!
    ## Chose the unfilled square s with the fewest possibilities
    n, s = min((len(values[s]), s) for s in squares if len(values[s]) > 1)
    return some(search(assign(values.copy(), s, d))
                for d in values[s])


def worstSearch(values):  # Question 2
    "Using depth-first search and propagation without using the most of possibilities, try all possible values."
    if values is False:
        return False  ## Failed earlier
    if all(len(values[s]) == 1 for s in squares):
        return values  ## Solved!

    n, s = max((len(values[s]), s) for s in squares if len(values[s]) > 1) #most of possibilities candidats first

    return some(worstSearch(assign(values.copy(), s, d))
                for d in values[s])


def randomSearch(values):  # Question 2 second ordering
    "Using depth-first search and propagation without using the most of possibilities, try all possible values."
    if values is False:
        return False  ## Failed earlier
    if all(len(values[s]) == 1 for s in squares):
        return values  ## Solved!
    s = random.choice(squares) #choice randomly the next candidat
    n = values[s]
    return some(randomSearch(assign(values.copy(), s, d))
                        for d in values[s])




################ Utilities ################

def find_random_unit(values,unit): #insert random values to unit
    exist = set()
    notExist = list(digits)
    changed = list() #list of no fixed squares
    for s in unit: # add exist elems to set
        if values[s] != '0' and values[s] != '.':
            exist.update(values[s])
            notExist.remove(values[s])

    for s in unit: #random fill
        if(values[s] not in exist):
            d = random.choice(notExist)
            values[s] = d
            notExist.remove(d)
            changed.append(s)
    return changed

def fill_all_units(values): #insert units without considering rows and colonnes
                            #returns a list of changedSquares
    unit1 = units['A1'][2]
    unit2 = units['A4'][2]
    unit3 = units['A7'][2]
    unit4 = units['D1'][2]
    unit5 = units['D4'][2]
    unit6 = units['D7'][2]
    unit7 = units['G1'][2]
    unit8 = units['G4'][2]
    unit9 = units['G7'][2]

    unitList = [unit1,unit2,unit3,unit4,unit5,unit6,unit7,unit8,unit9] #all units 3*3
    changedSquares = [] #list of no fixed squares

    for unit in unitList:
        valuesChanged = find_random_unit(values,unit)
        changedSquares.append(valuesChanged)

    return changedSquares


def conflit2(values): #returns global number of conflicts
    rs = [cross(rows, c) for c in cols]
    rowsError = 0
    cls = [cross(r, cols) for r in rows]
    colsError = 0
    for u in rs:
        rowsError += check(values,u)

    for u in cls:
        colsError += check(values,u)
    return colsError + rowsError

def check(values,unit): # returns number of conflicts in an unit
    setDigit = set()
    num = 0
    for s in unit:
        d = values[s]
        if d not in setDigit:
            setDigit.add(d)
        else: num += 1
    return num

def find_random_neighbor(values,changed): #returns a random neighbor values

    copy = values.copy() #returnValue
    copyChange = changed
    for unit in changed: #removing 1 elem unit
        if len(unit) == 1:
            copyChange.remove(unit)
    if len(copyChange) == 0:
        return values #no change

    randomUnit = random.choice(copyChange)
    if len(randomUnit) > 1:
        a, b = random.sample(randomUnit,2)
        copy[a], copy[b] = copy[b] , copy[a]

    return copy

def hillclimbing(grid):

    randGrid, existingValues = randomizedGridFill(grid)
    baseContradiction = conflit2(randGrid) #number of global conflicts

    for unit in unitlist[18:]:
        for i in range(len(unit)):
            if unit[i] not in existingValues:
                for j in range(len(unit)):
                    if unit[j] not in existingValues:
                        initValue = randGrid[unit[i]]
                        randGrid[unit[i]] = randGrid[unit[j]]
                        randGrid[unit[j]] = initValue

                        newContradiction = conflit2(randGrid)

                        if newContradiction < baseContradiction:
                            baseContradiction = newContradiction
                        else:
                            randGrid[unit[j]] = randGrid[unit[i]]
                            randGrid[unit[i]] = initValue
    return randGrid


def randomizedGridFill(grid): #returns a pair of (values, fixed squares) same function as fill_all_units, but different structure
    existingValues = {}
    for unit in unitlist[18:]:
        availableValues = [i for i in digits]
        for k in unit:
            if grid[k] in availableValues:
                existingValues[k] = grid[k]
                availableValues.remove(grid[k])
            if grid[k] == "0":
                val = random.choice(availableValues)
                grid[k] = str(val)
                availableValues.remove(val)

    return grid, existingValues #pair of values and existingValues(same as fill_all_units)


def calcul_temp(t_start):
    alpha = 0.95
    return t_start * alpha

def simulated_annealing(values):
    temp = 6.0 #initial temperature
    changed = fill_all_units(values) #random fill
    epsilon = 5e-100 #very small temp, for end condition of while boucle

    while temp > epsilon:
        temp = calcul_temp(temp)

        if temp == 0:
            return values
        next = find_random_neighbor(values, changed)# values of a random neighbor
        diffrence = conflit2(values) - conflit2(next)

        if diffrence > 0:
            values = next
        else:
            random_float = random.random()
            probability = math.exp(diffrence/temp)
            if random_float <= probability:
                values = next

            if conflit2(values) == 0:
                return values
    return False





def naked_pairs(values): #heuristique, nakedPair
    for unit in unitlist:
        pairs = [s for s in unit if len(values[s]) == 2] #number of possibilities = 2 in a square
        pair_values = [(values[s], s) for s in pairs] #
        seen = {} #check dict
        for pair in pair_values: #check
            if pair[0] in seen:
                for s in unit:
                    if s != pair[1] and s != seen[pair[0]]: #found
                        for digit in pair[0]:
                            eliminate(values,s,digit) #eliminate values[s] in peer squares
            else:
                seen[pair[0]] = pair[1]
    return values

def naked_triples(values): #heuristique nakedTriple
    for unit in unitlist:
        checkResult = check3candidates(values, unit) #check if there is any naked triple in the same unit

        if checkResult == False : #0 solution
            return values
        else: #there are some nakedTriple to eliminate
            for tripleSquares in checkResult:
                othersquares = set(unit).copy()
                digitResult = set()

                for s in tripleSquares:
                    othersquares.discard(s)
                    digitResult.add(values[s])
                all(eliminate(values,s,d) for s in othersquares for d in digitResult) #eliminate

    return values


def check3candidates(values, unit): #check 3 candidates of naked triple
    possibleResult = {} #dict(digit: square)
    squaresResult = []

    for s in unit:
        digit = values[s]
        if len(digit) <= 3 and len(digit) > 1: #number of possible result
            if digit not in possibleResult:
                possibleResult[digit] = [s]
            else:
                possibleResult[digit].append(s)

    #Check if there is any triple in possibleResult
    for d, sList in possibleResult.items():
        if len(d) == 3 and len(sList) == 3:

            for s in sList:
                squaresResult.append(s)

    if len(squaresResult) == 0:
        return False

    return squaresResult


def some(seq):
    "Return some element of seq that is true."
    for e in seq:
        if e: return e
    return False


def from_file(filename, sep='\n'):
    "Parse a file into a list of strings, separated by sep."
    return file(filename).read().strip().split(sep)


def shuffled(seq):
    "Return a randomly shuffled copy of the input sequence."
    seq = list(seq)
    random.shuffle(seq)
    return seq


################ System test ################

import time, random
import math

def solve_all(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids. Report results.
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""

    def time_solve(grid):
        start = time.clock()
        values = solve(grid)
        t = time.clock() - start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print ('(%.2f seconds)\n' % t)
        return (t, solved(values))

    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print ("Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times) / N, N / sum(times), max(times)))


def solve_all_using_worst(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids in the worst order. Report results.
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""

    def time_solve(grid):
        start = time.clock()
        values = worst_solve(grid)
        t = time.clock() - start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print ('(%.2f seconds)\n' % t)
        return (t, solved(values))

    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print ("Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times) / N, N / sum(times), max(times)))


def solve_all_using_random(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids in the worst order. Report results.
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""

    def time_solve(grid):
        start = time.clock()
        values = random_solve(grid)
        t = time.clock() - start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print ('(%.2f seconds)\n' % t)
        return (t, solved(values))

    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print ("Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times) / N, N / sum(times), max(times)))

def solve_all_using_nakedPair(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids. Report results. Using nakedPair heuristique
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""

    def time_solve(grid):
        start = time.clock()
        values = nakedPair_solve(grid)
        t = time.clock() - start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print ('(%.2f seconds)\n' % t)
        return (t, solved(values))

    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print ("Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times) / N, N / sum(times), max(times)))


def solve_all_using_nakedTriple(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids. Report results. using nakedTriple heuristique
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""

    def time_solve(grid):
        start = time.clock()
        values = nakedTriple_solve(grid)
        t = time.clock() - start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print ('(%.2f seconds)\n' % t)
        return (t, solved(values))

    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print ("Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times) / N, N / sum(times), max(times)))

def solve_all_using_annealing(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids. Report results. using simulated annealing
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""


    def time_solve(grid):
        start = time.clock()

        values = grid_values(grid1)
        values = simulated_annealing(values)

        t = time.clock() - start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print ('(%.2f seconds)\n' % t)

        return (t, solved(values))

    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print ("Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times) / N, N / sum(times), max(times)))


def solve_all_using_hill(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids. Report results. Using hill-climbing
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""


    def time_solve(grid):
        start = time.clock()
        values = grid_values(grid)
        values = hillclimbing(values)

        t = time.clock() - start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print ('(%.2f seconds)\n' % t)

        return (t, solved(values))

    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print ("Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times) / N, N / sum(times), max(times)))

def solved(values):
    "A puzzle is solved if each unit is a permutation of the digits 1 to 9."

    def unitsolved(unit): return set(values[s] for s in unit) == set(digits)

    return values is not False and all(unitsolved(unit) for unit in unitlist)


def random_puzzle(N=17):
    """Make a random puzzle with N or more assignments. Restart on contradictions.
    Note the resulting puzzle is not guaranteed to be solvable, but empirically
    about 99.8% of them are solvable. Some have multiple solutions."""
    values = dict((s, digits) for s in squares)
    for s in shuffled(squares):
        if not assign(values, s, random.choice(values[s])):
            break
        ds = [values[s] for s in squares if len(values[s]) == 1]
        if len(ds) >= N and len(set(ds)) >= 8:
            return ''.join(values[s] if len(values[s]) == 1 else '.' for s in squares)
    return random_puzzle(N)  ## Give up and make a new puzzle


grid1 = '003020600900305001001806400008102900700000008006708200002609500800203009005010300'
grid2 = '4.....8.5.3..........7......2.....6.....8.4......1.......6.3.7.5..2.....1.4......'
hard1 = '.....6....59.....82....8....45........3........6..3.54...325..6..................'


if __name__ == '__main__':
    test()

    # tests for easy puzzles
    print "easy puzzles:"
    solve_all(from_file("100sudoku.txt"), "100sudoku", None)  # Q1 test
    solve_all(from_file("1000sudoku.txt"), "1000sudoku",None) #Q2 test
    solve_all([random_puzzle() for _ in range(1000)], "random", 100.0) #Q1 random puzzle
    solve_all_using_random(from_file("100sudoku.txt"), "100random", None)  # Q2 test1
    solve_all_using_worst(from_file("100sudoku.txt"), "100worst", None)  # Q2 test2
    print

    # tests for hard puzzles.
    print "hard puzzles:"
    solve_all(from_file("top95.txt"), "95top", None)  # Q1
    solve_all_using_random(from_file("top95.txt"), "95random", None)  # Q2 test1
    solve_all_using_worst(from_file("top95.txt"), "95worst", None)  # Q2 test2

    print
    #tests for nakedPair and nakedTriple
    print "heuristique researches:"
    solve_all(from_file("top95.txt"), "95top", None)  # Q3
    solve_all_using_nakedPair(from_file("top95.txt"), "95topNakedPair", None)  # Q3
    solve_all_using_nakedTriple(from_file("top95.txt"), "95topNakedTriple", None) #Q3

    print
    #tests for local research
    print "local researches:"
    solve_all_using_hill(from_file("100sudoku.txt"), "100hill", None)
    solve_all_using_annealing(from_file("100sudoku.txt"), "100annealing", None )
    solve_all_using_annealing(from_file("top95.txt"), "95topAnnealing", None)
## References used:
## http://www.scanraid.com/BasicStrategies.htm
## http://www.sudokudragon.com/sudokustrategy.htm
## http://www.krazydad.com/blog/2005/09/29/an-index-of-sudoku-strategies/
## http://www2.warwick.ac.uk/fac/sci/moac/currentstudents/peter_cock/python/sudoku/
