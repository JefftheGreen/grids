# /usr/bin/env python
# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
import collections
import new_math as math
from PyQt4 import QtGui, QtCore
import warnings
from . import drawing


# Abstract base class for representing grids.
#
#   Attributes:
#       face_tagsets, edge_tagsets, vertex_tagsets:
#           Dictionaries.
#
#           Contain tagsets for face, edge, and vertex elements, respectively.
#
#           Have the form name: {elem: value}.
#               name is an arbitrary identifier for the tag.
#
#               elem is the tuple coordinates of an element.
#
#               value is the value of the tag for elem.
#
#           Tags store information about grid elements that is relatively
#           static and that is independent of any other information. For
#           example, the terrain of a grid square, the product of its
#           coordinates, or a unique identifier.
#
#           Tags for an element can be retrieved with self.tags()
#
#
#       property_maps:
#           Dictionary.
#
#           Contains property maps.
#
#           Has the form {category:
#                           {name:
#                               {property:
#                                   (value, priority)}}}
#               category is an arbitrary identifier. Different functions or
#               methods may use different categories. For example, drawing
#               functions may use 'tkinter' or 'qt', while pathfinding
#               algorithms may use 'passable' or 'cost'.
#
#               name is an arbitrary identifier for the property map.
#
#               property is the name of the property, e.g. 'color'.
#
#               value is a rule used to determine the property value for an
#               element.
#                   It may be a callable object or a sequence.
#
#                   If a callable object, it is called with the results of
#                   self.tags[element] as an argument and should return the
#                   value of property for element.
#
#                   If a sequence, it should have a callable object as its
#                   first element should be a callable object and following
#                   elements should be strings. The strings are names of
#                   instance variables or methods of the Grid or one of a few
#                   variables in the properties() method. The first element is
#                   called with the variables indicated by the following
#                   strings as arguments and should return the value of the
#                   property for element.
#
#                   If the callable object returns None, the value-priority
#                   pair is ignored.
#
#                   In addition to Grid instance variables, the following
#                   variables may be indicated:
#                       element: the coordinates for the element
#
#               priority is a rule used to determine the priority value of the
#               associated property value. Priority values are used to
#               determine which property value is used if multiple property
#               maps return different values for the same property; the
#               property value with the greatest priority value is used.
#                   It may use the same formats as value, but must return
#                   an object comparable with 0.
#
#           Property maps are used to determine information about grid
#           elements that is dynamic and dependent on other information. For
#           example, movement cost, passability, or drawing options.
#
#           Properties for an element can be accessed using self.properties()


class Grid(metaclass=ABCMeta):
    # Values are the length of coordinates used for internal representation.
    COORD_LEN = {'face': 2, 'edge': 3, 'vertex': 3}

    # Values should be lists of identifiers acceptable for that element type.
    # Two types should not share an identifier.
    COORD_DIRS = {'face': [], 'edge': [], 'vertex': []}

    # Keys should be direction identifiers, values should be angles clockwise
    # from the vector (0, 0) -> (1,0), in radians
    DIR_ANGLES = {}

    # Get the "proper" name of an element, as used internally.
    #   element:
    #       the coordinates of the element. sequence
    # Returns the coordinates of the element as a tuple.
    @classmethod
    @abstractmethod
    def proper(cls, element):
        pass

    # Get the pixel coordinates of the center of a face.
    #   face:
    #       the coordinates of the face. sequence.
    # Returns the pixel coordinates as a tuple of floats.
    @staticmethod
    @abstractmethod
    def to_pixel_face_center(face):
        pass

    # Get the angle of an edge clockwise from horizontal, in radians.
    @staticmethod
    @abstractmethod
    def to_pixel_face_vertices(face):
        pass

    # Get the pixel coordinates of the midpoint of an edge.
    #   edge:
    #       the coordinates of the edge. sequence.
    # Returns the pixel coordinates as a tuple of floats
    @staticmethod
    @abstractmethod
    def to_pixel_edge_mid(edge):
        pass

    # Get the pixel coordinates of a vertex.
    #   vertex:
    #       the coordinates of the vertex. sequence.
    # Returns the pixel coordinates as a tuple of floats
    @staticmethod
    @abstractmethod
    def to_pixel_vertex(vertex):
        pass

    # Get the face a pixel is in.
    #   vertex:
    #       the coordinates of the pixel. two-element sequence of floats or
    #       integers.
    #   round:
    #       the round method. acceptable values determined by subclass;
    #       usually a decimal rounding mode.
    # Returns the face coordinates as a tuple.
    @staticmethod
    @abstractmethod
    def to_face_pixel(vertex, round):
        pass


    @classmethod
    @abstractmethod
    def get_adjacency_map(cls, element_type, neighbor_type):
        pass

    # Get the equivalent coordinates of an element (usually edge or vertex).
    #   element:
    #       the coordinates of the element. sequence
    # Returns a tuple of element coordinates (tuples).
    # e.g. ((0,0,NE), (0,0,NW)).
    @classmethod
    @abstractmethod
    def aliases(cls, element):
        pass

    @classmethod
    @abstractmethod
    def manhattan(cls, element1, element2):
        pass

    @classmethod
    @abstractmethod
    def euclidian(cls, element1, element2):
        pass


    @classmethod
    def element_type(cls, element):
        if len(element) == cls.COORD_LEN['face']:
            if element[-1] in cls.COORD_DIRS['face'] or \
                            len(cls.COORD_DIRS['face']) == 0:
                return 'face'
        if len(element) == cls.COORD_LEN['edge']:
            if element[-1] in cls.COORD_DIRS['edge'] or \
                            len(cls.COORD_DIRS['edge']) == 0:
                return 'edge'
        if len(element) == cls.COORD_LEN['vertex']:
            if element[-1] in cls.COORD_DIRS['vertex'] or \
                            len(cls.COORD_DIRS['vertex']) == 0:
                return 'vertex'

    @classmethod
    def neighbors(cls, element, neighbor_type, include_dir=False):
        adj_map = cls.get_adjacency_map(element, neighbor_type)
        if adj_map is not None:
            if element[-1] in sum(cls.COORD_DIRS.values(), []):
                elem_coord = element[:-1]
            else:
                elem_coord = element
            adjacent = {}
            for adj in adj_map:
                acr, d = adj_map[adj]
                if adj[-1] in sum(cls.COORD_DIRS.values(), []):
                    adj = math.psum(elem_coord, adj[:-1]) + (adj[-1],)
                else:
                    adj = math.psum(elem_coord, adj)
                if acr[-1] in sum(cls.COORD_DIRS.values(), []):
                    acr = math.psum(elem_coord, acr[:-1]) + (acr[-1],)
                else:
                    acr = math.psum(elem_coord, acr[:-1])
                adjacent[adj] = acr, d if include_dir else acr
            return adjacent

    @classmethod
    def neighbor(cls, element, d):
        element = cls.proper(element)
        element_type = cls.element_type(element)
        try:
            if element_type == 'face':
                return math.psum(cls.NEIGHBORS[element_type][d], element)
            elif element_type in ['edge', 'vertex']:
                element_coords, element_dir = element[:-1], element[-1]
                diff = cls.NEIGHBORS[element_type][element_dir][d]
                num_diff, new_dir = diff[:-1], diff[-1]
                return tuple(math.psum(num_diff, element_coords) + (new_dir,))
        except KeyError:
            return

    # Initialize the grid.
    @abstractmethod
    def __init__(self, holes = ()):
        # Empty containers for tagsets.
        self.face_tagsets = {}
        self.edge_tagsets = {}
        self.vertex_tagsets = {}
        # Helper object to retrieve tags. Usage: self.tags[element_coords]
        self.tags = TagFetcher(self)
        # Empty container for property maps
        self.property_maps = {'tkinter': {},
                              'GTK': {},
                              'Qt': {},
                              'cost':
                                  {'default':
                                       {'cost':
                                            (lambda x: 1, -float('inf'))}},
                              'passable':
                                  {'default':
                                       {'passable':
                                            (lambda x: True, -float('inf'))}}}
        # Hidden variables behind all_faces/edges/vertices properties.
        self._all_faces = []
        self._all_edges = []
        self._all_vertices = []
        # Hidden variable behind bounds property
        self._bounds = []
        # Set to True when size or shape of Grid may have changed.
        self.need_to_calc_elements = True
        self.need_to_calc_bounds = True
        self._a_bbox = None
        self._o_bbox = None
        self.need_to_calc_bbox = True
        self.holes = list(holes)
        self.tagsets_frozen = False
        self.property_maps_frozen = False

    @property
    def all_faces(self):
        if self.need_to_calc_elements:
            self.calc_all_elements()
        return self._all_faces

    @property
    def all_edges(self):
        if self.need_to_calc_elements:
            self.calc_all_elements()
        return self._all_edges

    @property
    def all_vertices(self):
        if self.need_to_calc_elements:
            self.calc_all_elements()
        return self._all_vertices

    @property
    def bounds(self):
        if self.need_to_calc_bounds:
            self.calc_bounds()
        return self._bounds

    @property
    def o_bbox(self):
        if self.need_to_calc_bbox:
            self.calc_bbox()
        return self._o_bbox

    @property
    def a_bbox(self):
        if self.need_to_calc_bbox:
            self.calc_bbox()
        return self._a_bbox

    # Calculate all elements that exist and puts them in the lists
    # self._all_faces, self._all_edges, and self._all_vertices as tuples.
    @abstractmethod
    def calc_all_elements(self):
        pass

    def interior(self, element):
        element_type = self.element_type(element)
        adj = self.neighbors(element, 'face')
        for face in adj:
            if not face in self.all_faces:
                return False
        return True

    def calc_bounds(self):
        bounds = []
        for vertex in self.all_vertices:
            if not self.interior(vertex):
                bounds.append(vertex)
        self.need_to_calc_bounds = False
        self.need_to_calc_bbox = True
        self._bounds = bounds

    def calc_bbox(self):
        boundary_points = [self.to_pixel_vertex(p) for p in self.bounds]
        self._o_bbox = math.bounding_box(boundary_points, False)
        self._a_bbox = math.bounding_box(boundary_points, True)
        self.need_to_calc_bbox = False

    # Add a new tagset or update an existing one.
    #   name:
    #       name of the tagset. if a tagset of the same name already exists,
    #       it will be updated with the given values. usually string.
    #   element_type:
    #       the type of element the tagset contains. 'face', 'vertex',
    #       or 'edge'.
    #   tags:
    #       the tagset to be added. dictionary of the form:
    #           {element1: value1, element2: value2}
    #       elements should be tuples. values are arbitrary.
    # Raises ValueError if tags contains elements of the wrong type.
    # Raises ValueError if tags contains different names for the same element
    # with different values.
    def add_tagset(self, name, element_type, tags):
        if self.tagsets_frozen:
            warnings.warn('tried to add to frozen tagsets', RuntimeWarning)
            return
        size = self.COORD_LEN[element_type]
        dirs = self.COORD_DIRS[element_type]
        target = {'vertex': self.vertex_tagsets,
                  'edge': self.edge_tagsets,
                  'face': self.face_tagsets}[element_type]
        tagset = {}
        for element in tags:
            # Check that element is of element_type
            if not len(element) == size:
                raise ValueError('%s is not the correct length for a \
                                 %s' % (element, element_type))
            if not (element[-1] in dirs
                    or len(self.COORD_DIRS[element_type]) == 0):
                raise ValueError('%s is not a permissible direction for a \
                                 %s' % (element[-1], element_type))
            if self.proper(element) in tagset:
                # Raise Value Error if a different value was given to the
                # element under a different name.
                if tagset[self.proper(element)] != tags[element]:
                    raise ValueError('%s given two different values under \
                                     different names.' % self.proper(element))
            else:
                # Use only "official" name for element
                tagset[self.proper(element)] = tags[element]
        if name in target:
            target[name].update(tagset)
        else:
            target[name] = collections.defaultdict(lambda: None, tagset)

    # Add a new property map or update an existing one.
    #   category:
    #       the category of the property map. different categories are used
    #       by different processes. examples include 'tkinter', 'passable',
    #       and 'cost'.
    #   name:
    #       the name of the property map. arbitrary.
    #   prop_map:
    #       the property map. dictionary of form:
    #           {prop1: (value1, priority1), prop2: (value2, priority2)}
    #       see the description of the Grid class for more information.
    def add_property_map(self, category, name, prop_map):
        if self.property_maps_frozen:
            warnings.warn('tried to add to frozen property map',
                          RuntimeWarning)
            return
        for prop in prop_map:
            # Check that prop is a 2-length tuple
            if not isinstance(prop_map[prop], tuple):
                raise TypeError('Value for property %s is not a tuple, \
                                instead %s.' % (prop, prop_map[prop]))
            elif len(prop_map[prop]) != 2:
                raise ValueError('Value for property %s does not have two \
                                 elements, instead %s.'
                                 % (prop, prop_map[prop]))
            value, priority = prop_map[prop]
            # Check that value is callable or has a callable first element.
            while True:
                if callable(value):
                    break
                elif isinstance(value, collections.Sequence):
                    if callable(value[0]):
                        break
                raise ValueError('Value for property %s is not callable and \
                                 does not have a callable first element.')
            # Check that priority is callable or has a callable first element.
            while True:
                if callable(priority):
                    break
                elif isinstance(priority, collections.Sequence):
                    if callable(priority[0]):
                        break
                else:
                    try:
                        priority < 0
                        priority > 0
                        priority == 0
                        break
                    except TypeError:
                        pass
                raise ValueError('Priority for property %s is not callable \
                                 and does not have a callable first element.')
        # Add property map or update existing with new values.
        if category in self.property_maps:
            if name in self.property_maps:
                self.property_maps[category][name].update(prop_map)
            else:
                self.property_maps[category][name] = collections.defaultdict(
                    lambda: None, prop_map)
        else:
            self.property_maps[category] = {name: collections.defaultdict(
                lambda: None, prop_map)}

    # Get properties of an element.
    #   element:
    #       the coordinates of the element. tuple.
    #   category:
    #       the category of property maps to search in.
    def properties(self, element, category):
        # Can append to non-specified entries of properties.
        properties = collections.defaultdict(lambda: [])
        # Names of variables other than instance/class available to functions.
        variables = {'element': element, 'grid': self}
        for name in self.property_maps[category]:
            property_map = self.property_maps[category][name]
            for prop in property_map:
                value_func, priority_func = property_map[prop]
                # Call value_func to get value.
                if callable(value_func):
                    value = value_func(self.tags[element])
                else:
                    args = list(value_func[1:])
                    func = value_func[0]
                    # Search instance variables/methods, then variables set
                    # at beginning of method, for value_func arguments.
                    for i in range(len(args)):
                        if args[i] in dir(self):
                            args[i] = getattr(self, args[i])
                        elif args[i] in variables:
                            args[i] = variables[args[i]]
                        else:
                            raise ValueError('%s is not a valid variable name \
                                             to pass to a property value.')
                    value = func(*args)
                # Call value_func to get value.
                if callable(priority_func):
                    priority = priority_func(self.tags[element])
                elif isinstance(priority_func, collections.Sequence):
                    args = list(priority_func[1:])
                    func = priority_func[0]
                    # Search instance variables/methods, then variables set
                    # at beginning of method, for priority_func arguments.
                    for i in range(len(args)):
                        if args[i] in dir(self):
                            args[i] = getattr(self, args[i])
                        elif args[i] in variables:
                            args[i] = variables[args[i]]
                        else:
                            raise ValueError('%s is not a valid variable name \
                                             to pass to a property priority.')
                    priority = func(*args)
                else:
                    priority = priority_func
                if value is not None:
                    properties[prop].append((value, priority))
        highest_properties = {}
        # Get only the value with the highest priority.
        for prop in properties:
            highest_properties[prop] = sorted(properties[prop],
                                              key=lambda x: -x[1])[0][0]
        return highest_properties

    def draw_qt(self, scale=0.7, rotate=False):
        import sys
        app = QtGui.QApplication(sys.argv)
        viewer = QtGui.QGraphicsView()
        scene = drawing.GridScene(viewer, self, rotate)
        viewer.setScene(scene)
        viewer.setRenderHints(QtGui.QPainter.Antialiasing)
        viewer.scale(scale, scale)
        viewer.show()
        sys.exit(app.exec_())


# Helper class for retrieving tags of grid elements.
# Usage is tagfetcher[index], where index is the coordinates (as a tuple) of
# the element. For example (with conventional naming), tagfetcher[0,0] would
# retrieve the tags of face 0,0, while tagfetcher[0,0,Grid.N] would retrieve
# tags for the north edge of face 0,0.
class TagFetcher():
    # Initialize the tag fetcher
    def __init__(self, parent):
        self.parent = parent

    def __getitem__(self, index):
        type = self.parent.element_type(index)
        if type == 'face':
            return self.tags(index, self.parent.face_tagsets)
        elif type == 'edge':
            return self.tags(index, self.parent.edge_tagsets)
        elif type == 'vertex':
            return self.tags(index, self.parent.vertex_tagsets)
        raise IndexError('TagFetcher index out of range')

    def tags(self, index, tagsets):
        tags = {}
        for name in tagsets:
            for elem in tagsets[name]:
                if elem == index:
                    tags[name] = tagsets[name][elem]
                    break
        return tags