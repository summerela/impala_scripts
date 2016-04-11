"Identity job: Map, Reduce, and Combine as Classes.  Combiner is the same as the Reducer."
import hadoopy

class Mapper(object):
    def __init__(self):
        pass

    def map(self, key, value):
        yield key, value

class Reducer(object):
    def __init__(self):
        pass

    def reduce(self, key, values):
        for value in values:
            yield key, value


if __name__ == '__main__':
    hadoopy.run(Mapper, Reducer, Reducer)
