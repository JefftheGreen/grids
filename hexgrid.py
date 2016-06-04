# /usr/bin/env python
# -*- coding: utf-8 -*-


import decimal

import numpy as np
import scipy as sp
import scipy.spatial

import new_math as math
from utilities.variables import Identifier
from . import grid


class HexGrid(grid.Grid):
    # Values are the length of coordinates used for internal representation.
    COORD_LEN = {'face': 2, 'edge': 3, 'vertex': 3}

    N, NE, SE, S, SW, NW, NNE, E, SSE, SSW, W, NNW = \
        [Identifier(name) for name in
         ('N', 'NE', 'SE', 'S', 'SW', 'NW', 'NNE', 'E', 'SSE', 'SSW', 'W',
          'NNW')]

    # Values should be lists of identifiers acceptable for that element type.
    # Two types should not share an identifier.
    # The edge and vertex identifiers should be arranged in such a fashion
    # that 'edge'[n] is between 'vertex'[n-1] and 'vertex'[n] and, conversely,
    # that 'vertex'[n] is between 'edge'[n] and 'edge'[n+1]. (Modulo the
    # respective lengths.)
    COORD_DIRS = {'face': [],
                  'edge': [N, NE, SE, S, SW, NW],
                  'vertex': [NNE, E, SSE, SSW, W, NNW]}

    RELATIVE_DIRS = {'face': [N, NE, SE, S, SW, NW],
                     'vertex': {E: [NNW, E, SSW],
                                W: [NNE, W, SSE]}}

    ALL_DIRS = [N, NE, SE, S, SW, NW, NNE, E, SSE, SSW, W, NNW]

    PROPER_DIRS = [N, NE, SE, E, W]

    NEIGHBORS = {'face': {N: (0, 1),
                          S: (0, -1),
                          NE: (1, 0),
                          SW: (-1, 0),
                          SE: (1, -1),
                          NW: (-1, 1)},
                 'edge': {N: {NNE: (0, 1, SE),
                              SSE: (0, 0, NE),
                              SSW: (-1, 1, SE),
                              NNW: (-1, 1, NE)},
                          NE: {NNE: (0, 1, SE),
                               W: (0, 0, N),
                               E: (1, -1, N),
                               SSW: (0, 0, SE)},
                          SE: {NNW: (0, 0, NE),
                               E: (1, -1, N),
                               W: (0, -1, N),
                               SSE: (0, -1, NE)}},
                 'vertex': {E: {NNW: (1, 0, W),
                                E: (2, -1, W),
                                SSW: (1, -1, W)},
                            W: {NNE: (-1, 1, E),
                                W: (-2, 1, E),
                                SSE: (-1, 0, E)}}}

    ADJACENCY = {'face': {'face': {(0, 1): ((0, 0, N), N),
                                   (1, 0): ((0, 0, NE), NE),
                                   (1, -1): ((0, 0, SE), SE),
                                   (0, -1): ((0, -1, N), S),
                                   (-1, 0): ((-1, 0, NE), SW),
                                   (-1, 1): ((-1, 1, SE), SW)},
                          'edge': {(0, 0, NE): ((1, 0), NE),
                                   (-1, 0, NE): ((-1, 0), SW),
                                   (-1, 1, SE): ((-1, 1), NW),
                                   (0, -1, N): ((0, -1), S),
                                   (0, 0, N): ((0, 1), N),
                                   (0, 0, SE): ((1, -1), SE)},
                          'vertex': {(1, 0, W): ((0, 1, SE), NNE),
                                     (0, 0, E): ((1, -1, N), E),
                                     (1, -1, W): ((0, -1, NE), SSE),
                                     (-1, 0, E): ((-1, 0, SE), SSW),
                                     (0, 0, W): ((-1, 0, N), W),
                                     (-1, 1, E): ((-1, 1, NE), NNW)}},
                 'edge': {'face': {N: {(0, 0): ((0, -1, N), S),
                                       (0, 1): ((0, 1, N), N)},
                                   NE: {(0, 0): ((-1, 0, NE), SW),
                                        (1, 0): ((1, 0, NE), NE)},
                                   SE: {(0, 0): ((-1, 1, SE), NW),
                                        (1, -1): ((1, -1, SE), SE)}},
                          'edge': {N: {(0, 0, NE): ((1, 0, W), SSW),
                                       (0, 1, SE): ((1, 0, W), NNW),
                                       (-1, 1, NE): ((-1, 1, E), NNW),
                                       (-1, 1, SE): ((-1, 1, E), SSW)},
                                   NE: {(0, 0, N): ((1, 0, W), W),
                                        (0, 1, SE): ((1, 0, W), NNW),
                                        (1, -1, N): ((0, 0, E), E),
                                        (0, 0, SE): ((0, 0, E), SSW)},
                                   SE: {(0, 0, NE): ((0, 0, E), NNW),
                                        (1, -1, N): ((0, 0, E), E),
                                        (0, -1, N): ((1, -1, W), W),
                                        (0, -1, NE): ((1, -1, W), SSE)}},
                          'vertex': {N: {(1, 0, W): ((1, 0), E),
                                         (-1, 1, E): ((-1, 1), W)},
                                     NE: {(1, 0, W): ((0, 1), NNW),
                                          (0, 0, E): ((1, -1), SSE)},
                                     SE: {(0, 0, E): ((1, 0), NNE),
                                          (1, -1, W): ((0, -1), SSW)}}},
                 'vertex': {'face': {E: {(0, 0): ((0, 0, W), W),
                                         (1, 0): ((2, 0, W), NNE),
                                         (1, -1): ((2, -2, W), SSE)},
                                     W: {(0, 0): ((0, 0, E), E),
                                         (-1, 1): ((-2, 2, E), NNW),
                                         (-1, 0): ((-2, 0, E), SSW)}},
                            'edge': {E: {(0, 0, NE): ((1, 0, W), NNW),
                                         (0, 0, SE): ((1, -1, W), SSW),
                                         (1, -1, N): ((2, -1, W), E)},
                                     W: {(-1, 0, N): ((-2, 1, E), W),
                                         (-1, 0, NE): ((-1, 0, E), SSE),
                                         (-1, 1, SE): ((-1, 1, E), NNE)}},
                            'vertex': {W: {(-1, 0, E): ((-1, 0, NE), SSE),
                                           (-2, 1, E): ((-1, 0, N), W),
                                           (-1, 1, E): ((-1, 1, SE), NNE)},
                                       E: {(1, 0, W): ((0, 0, NE), NNW),
                                           (1, -1, W): ((0, 0, SE), SSW),
                                           (2, -1, W): ((1, -1, N), E)}}
                            }}

    OFF_FROM_CENTER = {N: (0, math.sqrt(3) / 2),
                       NE: (0.75, math.sqrt(3) / 4),
                       SE: (0.75, -math.sqrt(3) / 4),
                       E: (1, 0),
                       W: (-1, 0)}

    DIR_ANGLES = {N: -math.pi/2,
                  NE: -math.pi/6,
                  SE: math.pi/6,
                  S: math.pi/2,
                  SW: 5*math.pi/6,
                  NW: -5*math.pi/6,
                  NNE: -math.pi/3,
                  E: 0,
                  SSE: math.pi/3,
                  SSW: 2*math.pi/3,
                  W: math.pi,
                  NNW: -2*math.pi/3}

    @staticmethod
    def to_pixel_face_center(face):
        q, r = face
        x = 3 / 2 * q
        y = math.sqrt(3) * (r + q / 2)
        return x, y

    @staticmethod
    def to_face_pixel(pixel, rounding=None):
        x, y = pixel
        q = x * 2 / 3
        r = -x / 3 + math.sqrt(3) / 3 * y
        if rounding == 'half ceiling':
            q, r = ((int(decimal.Decimal(i).quantize(1,
                                                     decimal.ROUND_HALF_UP))
                     if i > 0 else
                     int(decimal.Decimal(i).quantize(1,
                                                     decimal.ROUND_HALF_DOWN)))
                    for i in (q, r))
        elif rounding == 'half floor':
            q, r = ((int(decimal.Decimal(i).quantize(1,
                                                     decimal.ROUND_HALF_UP))
                     if i < 0 else
                     int(decimal.Decimal(i).quantize(1,
                                                     decimal.ROUND_HALF_DOWN)))
                    for i in (q, r))
        elif rounding is not None:
            q, r = (int(decimal.Decimal(i).quantize(1, rounding))
                    for i in (q, r))
        return q, r

    @classmethod
    def to_pixel_face_vertices(cls, face):
        vertices = cls.neighbors(face, 'vertex')
        vertices = [cls.to_pixel_vertex(v) for v in vertices]
        hull = sp.spatial.ConvexHull(np.array(vertices))
        poly = hull.points[hull.vertices]
        poly = [tuple(p) for p in poly]
        return tuple(poly)

    @classmethod
    def to_pixel_edge_mid(cls, edge):
        center = cls.to_pixel_face_center(edge[:-1])
        return math.psum(center, cls.OFF_FROM_CENTER[edge[-1]])

    @classmethod
    def to_pixel_edge_ends(cls, edge):
        p1, p2 = tuple(cls.neighbors(edge, 'vertex').keys())
        p1, p2 = (cls.to_pixel_vertex(p) for p in (p1, p2))
        return p1, p2

    @classmethod
    def to_pixel_vertex(cls, vertex):
        center = cls.to_pixel_face_center(vertex[:-1])
        return math.psum(center, cls.OFF_FROM_CENTER[vertex[-1]])

    @classmethod
    def edge_angle(cls, edge):
        d = cls.proper(edge)[-1]
        if d == cls.N:
            return 0
        elif d == cls.NE:
            return math.pi/3
        elif d == cls.SE:
            return -math.pi/3

    @classmethod
    def get_adjacency_map(cls, element, neighbor_type):
        element_type = cls.element_type(element)
        if element_type == 'face':
            return cls.ADJACENCY['face'][neighbor_type]
        else:
            return cls.ADJACENCY[element_type][neighbor_type][element[-1]]

    @classmethod
    def aliases(cls, element):
        element_type = cls.element_type(element)
        element = cls.proper(element)
        if element_type is None:
            return
        elif element_type == 'face':
            return element,
        elif element_type == 'edge':
            if element[-1] == cls.N:
                return element, cls.neighbor(element[:-1], cls.N) + (cls.S,)
            elif element[-1] == cls.NE:
                return element, cls.neighbor(element[:-1], cls.NE) + (cls.SW,)
            else:
                return element, cls.neighbor(element[:-1], cls.SE) + (cls.NE,)
        elif element_type == 'vertex':
            if element[-1] == cls.E:
                return (element,
                        cls.neighbor(element[:-1], cls.NE) + (cls.SSW,),
                        cls.neighbor(element[:-1], cls.SE) + (cls.NNW,))
            else:
                return (element,
                        cls.neighbor(element[:-1], cls.NW) + (cls.SSE,),
                        cls.neighbor(element[:-1], cls.SW) + (cls.NNE,))

    @classmethod
    def proper(cls, element):
        element_type = cls.element_type(element)
        if element_type is None:
            return
        elif element_type == 'face' or element[-1] in cls.PROPER_DIRS:
            return element
        elif element_type == 'edge':
            if element[-1] == cls.S:
                return cls.neighbor(element[:-1], cls.S) + (cls.N,)
            elif element[-1] == cls.NW:
                return cls.neighbor(element[:-1], cls.NW) + (cls.SE,)
            else:
                return cls.neighbor(element[:-1], cls.SW) + (cls.NE,)
        elif element_type == 'vertex':
            if element[-1] == cls.NNE:
                return cls.neighbor(element[:-1], cls.NE) + (cls.W,)
            elif element[-1] == cls.SSE:
                return cls.neighbor(element[:-1], cls.SE) + (cls.W,)
            elif element[-1] == cls.SSW:
                return cls.neighbor(element[:-1], cls.SW) + (cls.E,)
            else:
                return cls.neighbor(element[:-1], cls.NW) + (cls.E,)

    @classmethod
    def manhattan(cls, element_a, element_b):
        if cls.element_type(element_b) != cls.element_type(element_a):
            raise ValueError('elements are not the same type')
        if cls.element_type(element_a) == 'face':
            qa, ra = element_a
            qb, rb = element_b
            sa, sb = -qa-ra, -qb-rb
            return int((abs(qa - qb) + abs(ra - rb) + abs(sa - sb)) / 2)
        if cls.element_type(element_b) == 'vertex':
            qa, ra, da = element_a
            qb, rb, db = element_b
            sa, sb = -qa-ra, -qb-rb
            dq = qa - qb
            dr = ra - rb
            ds = sa - sb
            if da == db:
                return max(2*abs(dq), 2*abs(dr), 2*abs(ds))
            else:
                if da == cls.E:
                    dq *= -1
                    dr *= -1
                    ds *= -1
                dq -= 2
                dr += 1
                ds += 1
                return (max(2*abs(dq), 2*abs(dr), 2*abs(ds)) + 1 if dq >= 0
                        else max(2*abs(dq), 2*abs(dr), 2*abs(ds)) - 1)


    @classmethod
    def euclidian(cls, element_a, element_b):
        if cls.element_type(element_a) == 'face':
            a = np.array(cls.to_pixel_face_center(element_a))
        elif cls.element_type(element_a) == 'edge':
            a = np.array(cls.to_pixel_edge_mid(element_a))
        elif cls.element_type(element_a) == 'vertex':
            a = np.array(cls.to_pixel_vertex(element_a))
        if cls.element_type(element_b) == 'face':
            b = np.array(cls.to_pixel_face_center(element_b))
        elif cls.element_type(element_b) == 'edge':
            b = np.array(cls.to_pixel_edge_mid(element_b))
        elif cls.element_type(element_b) == 'vertex':
            b = np.array(cls.to_pixel_vertex(element_b))
        if a is None or b is None:
            raise ValueError('element not recognized')
        return sum((a-b)**2)**0.5

    @classmethod
    def ring(cls, center, radius):
        if radius == 0:
            return [center]
        current = center
        members = []
        for i in range(radius):
            current = cls.neighbor(current, cls.NE)
        for d in (cls.NW, cls.SW, cls.S, cls.SE, cls.NE, cls.N):
            for i in range(radius):
                members.append(current)
                current = cls.neighbor(current, d)
        return members

    @classmethod
    def circle(cls, center, radius):
        members = [center]
        for r in range(radius):
            members += cls.ring(center, r + 1)
        return members

    @classmethod
    def rotate(cls, coords, center, degrees):
        if degrees % 60 != 0:
            raise AttributeError("degrees must be divisible by 60")
        q, r = coords
        s = -q-r
        cq, cr= center
        cs = -cq-cr
        newq, newr, news = q-cq, r-cr, s-cs
        for i in range(degrees // 60):
            newq, newr, news = -news, -newq, -newr
        newq, newr, = newq + cq, newr + cr
        return newq, newr

    @classmethod
    def line_between(cls, face1, face2, rounding='half ceiling'):
        distance = cls.manhattan(face1, face2)
        coord1, coord2 = (cls.to_pixel_face_center(f) for f in (face1, face2))
        face1, face2, coord1, coord2 = (np.array(f) for f in (face1,
                                                              face2,
                                                              coord1,
                                                              coord2))
        step = (coord2 - coord1) / distance
        faces = [coord1 + step * i for i in range(distance)] + [coord2]
        return tuple((cls.to_face_pixel(tuple(f), rounding) for f in faces))

    def __init__(self, qmax, qmin=None, rmax=None,
                 rmin=None, smax=None, smin=None, holes=()):
        qmin = -qmax if qmin is None else qmin
        rmax = qmax if rmax is None else rmax
        rmin = -rmax if rmin is None else rmin
        smax = qmax if smax is None else smax
        smin = -smax if smin is None else smin
        self.qrange = range(qmin, qmax + 1)
        self.rrange = range(rmin, rmax + 1)
        self.srange = range(smin, smax + 1)
        super().__init__(holes)

    def calc_all_elements(self):
        f, e, v = set(), set(), set()
        for q in self.qrange:
            for r in self.rrange:
                if -q - r in self.srange and not (q, r) in self.holes:
                    f.add((q, r))
        for face in f:
            vertices = self.neighbors(face, 'vertex')
            edges = self.neighbors(face, 'edge')
            for vertex in vertices:
                v.add(self.proper(vertex))
            for edge in edges:
                e.add(self.proper(edge))
        self._all_faces, self._all_edges, self._all_vertices = f, e, v

    def calc_bounds(self):
        bounds = []
        to_check = set()
        for h in self.holes:
            for v in self.neighbors(h, 'vertex'):
                to_check.add(v)
        # Set q constant
        for face in [(max(self.qrange), r) for r in self.rrange]:
            to_check.add(self.proper(face + (self.E,)))
        for face in [(min(self.qrange), r) for r in self.rrange]:
            to_check.add(self.proper(face + (self.W,)))
        # Set r constant
        for face in [(q, max(self.rrange)) for q in self.qrange]:
            to_check.add(self.proper(face + (self.NNW,)))
        for face in [(q, min(self.rrange)) for q in self.qrange]:
            to_check.add(self.proper(face + (self.SSE,)))
        # Set s constant
        for face in [(q, max(self.srange) - q) for q in self.qrange]:
            to_check.add(self.proper(face + (self.SSW,)))
        for face in [(q, min(self.srange) - q) for q in self.qrange]:
            to_check.add(self.proper(face + (self.NNE,)))
        for vertex in to_check:
            if vertex in self.all_vertices:
                if not self.interior(vertex):
                    bounds.append(vertex)
        self.need_to_calc_bounds = False
        self.need_to_calc_bbox = True
        self._bounds = bounds
