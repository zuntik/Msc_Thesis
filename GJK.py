import numpy as np
import matplotlib.pyplot as plt

class Shape:

    def getFarthestPointsInDirection(d):
        pass

    def __add__():
        pass

    def __sub__():
        pass

class Polygon(Shape):

    def __init__(self, points):
        self.points = points
        self.dim = points.shape[1]
        # points is N * dim 

    def getFarthestPointsInDirection(self, d):
        dot_prods = self.points @ np.array(d).reshape(-1,1)
        dot_prods= dot_prods.flatten()
        farthest = np.argwhere(dot_prods == np.amax(dot_prods)).flatten()
        return self.points[farthest,:]

    def plot(self):
        plt.scatter(self.points[:,0],self.points[:,1])


class Simplex:
    def __init__(self):
        self.list_of_points = list()

    def add(self,p):
        if not self.list_of_points:
            self.list_of_points.append(p)
            self.dim = len(p)
        else:
            self.list_of_points.append(p)

    def getA(self):
        return self.list_of_points[-1]

    def getB(self):
        return self.list_of_points[-2]

    def getC(self):
        return self.list_of_points[-3]

def tripleProduct(a,b,c):
    #print(a)
    #print(b)
    #print(c)
    a = a.flatten()
    b = b.flatten()
    c = c.flatten()
    #return b.dot(c.dot(a)) - a.dot(c.dot(b))
    return -a.dot(c.dot(b)) + b.dot(c.dot(a))


def do_simplex(simplex,d):
    print("simplex")
    a = simplex.getA()
    b = simplex.getB()
    ao = - a
    if len(simplex.list_of_points) == 3:
        c = simplex.getC()
        ab = b-a
        ac = c-a
        abPrep=tripleProduct(ac,ab,ab)
        acPrep=tripleProduct(ab,ac,ac)

        if abPrep.dot(ao) > 0:
            simplex.removeC()
            d = abPrep
        elif acPrep.dot(ao) > 0:
            simplex.removeB()
            d = acPrep
        else:
            return None

    else:
        ab = b-a
        abPrep = tripleProduct(ab,ao,ab)
        d = abPrep

    return d


def support(sa,sb,d):
    pointsa = sa.getFarthestPointsInDirection( d)
    pointsb = sb.getFarthestPointsInDirection(-d)
    return pointsa[0,:] - pointsb[0,:]

def supportv2(sa,sb,d):
    pointsa = sa.getFarthestPointsInDirection( d)
    pointsb = sb.getFarthestPointsInDirection(-d)
    all_supports = np.array([ pointa + (-pointb)  for pointa in pointsa for pointb in pointsb ])
    # find supports with smallest norm
    norms = np.linalg.norm(all_supports, axis=1)
    return all_supports[np.argmin(norms),:]

def GJK(s1,s2):

    # initial direction
    d = np.array([1 , 0 ]).reshape(-1,1)
    d = np.array([-1 , 1 ]).reshape(-1,1)

    # initial resulting support
    s = support(s1,s2,d)
    print("first support: "+ str(s.flatten()))

    # creat the simplex
    simplex = Simplex()
    simplex.add(s)

    # d is the negated s
    d = -s 

    while True:
        a = support(s1,s2,d)
        if a.dot(d) <= 0:
            return False
        simplex.add(a)

        d = do_simplex(simplex,d)

        if d is None:
            return True



if __name__== "__main__":
    s1 = Polygon(np.array([ [0,0], [0,2], [1,0], [1,1],[1,2] ]))
    s2 = Polygon(np.array([ [2,1], [2,2], [3,1],[3,2] ]))

    s1 = Polygon(np.array([ [1,0], [0,1], [2,1], [1,2] ]))
    s2 = Polygon(np.array([ [5,0], [4,1], [6,1], [5,2] ]))

    s1 = Polygon(np.array([ [1,0], [0,1], [2,1], [1,2] ]))
    s2 = Polygon(np.array([ [2,0], [1,1], [3,1], [2,2] ]))

    #s1 = Polygon(np.array([ [4,11], [9,9], [4,5] ]))
    #s2 = Polygon(np.array([ [5,7], [12,7], [10,2], [7,3] ]))

    s1.plot()
    s2.plot()

    plt.show()

    if GJK(s1,s2):
        print("They intersect")
    else:
        print("They do not intersect")


