# /usr/bin/env python
# -*- coding: utf-8 -*-

from . import grid
import new_math.packer as packer
import new_math as math
import collections
import numbers
from abc import ABCMeta


class LinkedGrid(grid.Grid):
    # Get the "proper" name of an element, as used internally.
    #   element:
    #       the coordinates of the element. sequence
    # Returns the coordinates of the element as a tuple.
    def proper(self, element):
        return (element[0],) + self[element[0]].proper(element[1:])

    # Get the pixel coordinates of the center of a face.
    #   face:
    #       the coordinates of the face. sequence.
    # Returns the pixel coordinates as a tuple of floats.
    def to_pixel_face_center(self, face):
        # Use the subgrid's method
        return self[face[0]].to_pixel_face_center(face[1:])

    def to_pixel_face_vertices(self, face):
        # Use the subgrid's method
        return self[face[0]].to_pixel_face_vertices(face[1:])

    def to_pixel_edge_ends(self, edge):
        # Use the subgrid's method
        return self[edge[0]].to_pixel_edge_ends(edge[1:])

    # Get the pixel coordinates of the midpoint of an edge.
    #   edge:
    #       the coordinates of the edge. sequence.
    # Returns the pixel coordinates as a tuple of floats
    def to_pixel_edge_mid(self, edge):
        # Use the subgrid's method
        return self[edge[0]].to_pixel_edge_mid(edge[1:])

    def edge_angle(self, edge):
        # Use the subgrid's method
        return self[edge[0]].edge_angle(edge[1:])

    # Get the pixel coordinates of a vertex.
    #   vertex:
    #       the coordinates of the vertex. sequence.
    # Returns the pixel coordinates as a tuple of floats
    def to_pixel_vertex(self, vertex):
        # Use the subgrid's method
        return self[vertex[0]].to_pixel_vertex(vertex[1:])

    # Get the face a pixel is in.
    #   vertex:
    #       the coordinates of the pixel. two-element sequence of floats or
    #       integers.
    #   round:
    #       the round method. acceptable values determined by subclass;
    #       usually a decimal rounding mode.
    # Returns the face coordinates as a tuple.
    def to_face_pixel(self, component, point, rounding):
        # Use the subgrid's method
        return self[component].to_face_pixel(point, rounding)

    def get_adjacency_map(self, component, element_type, neighbor_type):
        # Use the subgrid's method
        return self[component].get_adjacency_map(element_type, neighbor_type)

    # Get the equivalent coordinates of an element (usually edge or vertex).
    #   element:
    #       the coordinates of the element. sequence
    # Returns a tuple of element coordinates (tuples).
    # e.g. ((0,0,NE), (0,0,NW)).
    def aliases(self, element):
        # Use the subgrid's method
        base_aliases = self[element[0]].aliases(element[1:])
        # Rename the elements
        return tuple([(element[0],) + a for a in base_aliases])

    def manhattan(self, element1, element2):
        # If the elements are from the same subgrid, use the subgrid's.
        if element1[0] == element2[0]:
            return self[element1[0]].manhattan(element1[1:],
                                               element2[1:])
        else:
            # Infinity if on different subgrids
            return float('inf')

    def euclidian(self, element1, element2):
        # If the elements are from the same subgrid, use the subgrid's.
        if element1[0] == element2[0]:
            return self[element1[0]].euclidian(element1[1:],
                                               element2[1:])
        else:
            # Infinity if on different subgrids
            return float('inf')

    def element_type(self, element):
        # Use the subgrid's method
        return self[element[0]].element_type(element[1:])

    def neighbors(self, element, neighbor_type):
        # Get neighbors on the same subgrid
        neighbors = self[element[0]].neighbors(element[1:], neighbor_type)
        # Get neighbors via links
        if self.element_type(element) == 'face':
            for faces in self._links:
                if element in faces:
                    other_face = (set(faces) - {element}).pop()
                    # Figure out the order of the elements so we can get
                    # the correct edge
                    element_index = faces.index(element)
                    edge = self._links[faces][element_index]
                    if neighbor_type == 'face':
                        neighbors[other_face] = edge
                    elif neighbor_type == 'edge':
                        if edge in neighbors:
                            neighbors[edge] = neighbors[edge] + (other_face,)
                        else:
                            neighbors[edge] = other_face
        return neighbors

    def neighbor(self, element, d):
        # Use the subgrid's method
        return (element[0],) + self[element[0]].neighbor(element[1:], d)

    def element_parent_subgrid(self, element):
        return element[0]

    def __init__(self, subgrids, links):
        super().__init__()
        # Contains the subgrids
        self._subgrids = list(subgrids)
        # Contains the links.
        # Format (elem1, elem2): (connecting_elem1, connecting_elem2)
        self.links = Linkage(self)
        for l in links:
            source, sink = l
            source_dir, sink_dir = links[l]
            self.links.add_link(source, source_dir, sink, sink_dir)
        self.layout = Layout(self)
        self.layout.pack('2d pack')
        # Get valid coordinate directions from subgrids
        self.COORD_DIRS = {}
        for c in self:
            for elem_type in c.COORD_DIRS:
                for d in c.COORD_DIRS[elem_type]:
                    if not hasattr(self, str(d)):
                        setattr(self, str(d), d)
        # Merge tagsets and property maps, then freeze subgrids'
        self.consolidate_tagsets()
        self.consolidate_property_maps()

    # Can be accessed by integer indexes, returning subgrids, or slices,
    # returning links between the subgrids start and stop
    def __getitem__(self, index):
        if isinstance(index, int):
            try:
                return self._subgrids[index]
            except IndexError:
                raise IndexError('LinkedGrid index out of range')
        elif isinstance(index, slice):
            try:
                return self._links[(index.start, index.stop)]
            except IndexError:
                raise IndexError('LinkedGrid index out of range')

    # Use subgrids' Identifier objects (e.g. directions) if own aren't
    # initialized
    def __getattr__(self, name):
        for c in self:
            if hasattr(type(c), name):
                return getattr(type(c), name)
        raise AttributeError("'{0}' object has no attribute \
                             '{1}'".format(type(self), name))

    def __len__(self):
        return len(self._subgrids)

    # All elements from all subgrids
    def calc_all_elements(self):
        f, e, v = set(), set(), set()
        for c in range(len(self._subgrids)):
            for face in self[c].all_faces:
                f.add((c,) + face)
            for edge in self[c].all_edges:
                e.add((c,) + edge)
            for vertex in self[c].all_vertices:
                v.add((c,) + vertex)
        self._all_faces, self._all_edges, self._all_vertices = f, e, v

    def calc_bounds(self):
        self._bounds = [c.bounds for c in self]

    def calc_bbox(self):
        flat_bounds = []
        for c in range(len(self.bounds)):
            for b in self.bounds[c]:
                flat_bounds.append((c,) + b)
        boundary_points = [self.to_pixel_vertex(p) for p in flat_bounds]
        self._o_bbox = math.bounding_box(boundary_points, False)
        self._a_bbox = math.bounding_box(boundary_points, True)
        self.need_to_calc_bbox = False

    # Called at initialization to combine tagsets from subgrids, then freeze
    # the subgrids' tagsets.
    def consolidate_tagsets(self):
        for i in range(len(self)):
            c = self[i]
            for name in c.face_tagsets:
                tagset = c.face_tagsets[name]
                if name not in self.face_tagsets:
                    self.face_tagsets[name] = {}
                for element in tagset:
                    self.face_tagsets[name][(i,) + element] = tagset[element]
            for name in c.edge_tagsets:
                tagset = c.edge_tagsets[name]
                if name not in self.edge_tagsets:
                    self.edge_tagsets[name] = {}
                for element in tagset:
                    self.edge_tagsets[name][(i,) + element] = tagset[element]
            for name in c.vertex_tagsets:
                tagset = c.vertex_tagsets[name]
                if name not in self.vertex_tagsets:
                    self.vertex_tagsets[name] = {}
                for element in tagset:
                    self.vertex_tagsets[name][(i,) + element] = tagset[element]
            c.tagsets_frozen = True

    # Called at initialization to combine property maps from subgrids, then
    # freeze the subgrids' tagsets.
    def consolidate_property_maps(self):
        for c in self:
            for category in c.property_maps:
                if category not in self.property_maps:
                    self.property_maps[category] = {}
                for name in c.property_maps[category]:
                    map = c.property_maps[category][name]
                    new_name = str(name) + str(id(c))
                    self.property_maps[category][new_name] = map
            c.property_maps_frozen = True


class Linkage():
    def __init__(self, parent):
        self.parent = parent
        self._links = {}

    def add_link(self, source, source_dir, sink, sink_dir, directed=True):
        source_grid, source_elem = source[0], source[1:]
        sink_grid, sink_elem = sink[0], source[1:]
        if self.parent[source_grid].neighbor(source_elem, source_dir):
            if self.parent[sink_grid].neighbor(sink_elem, sink_dir):
                self._links[(source, sink)] = (source_dir, sink_dir)
                if not directed:
                    self._links[(sink, source)] = (sink_dir, source_dir)
            else:
                raise ValueError('invalid sink direction')
        else:
            raise ValueError('invalid source direction')

    def __getitem__(self, item):
        if not isinstance(item, slice):
            raise TypeError('Linkage index must be slice, not {0}'.format(
                type(item)))
        if isinstance(item.start, tuple) and isinstance(item.stop, tuple):
            if (item.start, item.stop) in self._links:
                return self._links[(item.start, item.stop)]
        else:
            results = {}
            for link in self._links:
                if ((link[0] == item.start or
                     link[0][0] == item.start or
                     item.start is None) and
                    (link[1] == item.stop or
                     link[1][0] == item.stop or
                     item.stop is None)):
                    results[link] = self._links[link]
            if len(results) > 0:
                return results

    def __iter__(self):
        for ss_pair in self._links:
            yield tuple(zip(ss_pair, self._links[ss_pair]))


# Object used to create and store layouts for LinkedGrids.
class Layout():
    # Initialize Layout object
    #   parent:
    #       LinkedGrid the Layout belongs to and describes.
    def __init__(self, parent):
        self.parent = parent
        self._positions = {}
        self.rotated = None
        self.current_method = None

    # Automatically produce a layout.
    #   method:
    #       the method to use. One of:
    #           '2d pack' - packs around the point 0,0 to minimize area
    #           '1d sorted width' - packs in order from narrow to wide
    #           '1d sorted height' - packs in order from short to tall
    #   rotate:
    #       if True, use subgrids' oriented bounding boxes. if False, use
    #       subgrid's aligned bounding boxes
    def pack(self, method, rotate=False):
        self.current_method = method
        if method == '2d pack':
            self.rotated = rotate
            p = packer.Packer()
            for i in range(len(self.parent)):
                bounds = [self.parent.to_pixel_vertex((i,) + v)
                          for v in self.parent.bounds[i]]
                bbox = math.bounding_box(bounds, not rotate)
                p.add_rect(bbox.width, bbox.height, data=i)
            packed = p.pack(padding=0.5)
            min_x = min([bbox.center().x for bbox in packed])
            min_y = min([bbox.center().y for bbox in packed])
            for bbox in packed:
                x = bbox.center().x - min_x
                y = bbox.center().y - min_y
                self._positions[bbox.data] = x, y
        elif method == '1d sorted width':
            bounds = [[self.parent.to_pixel_vertex((i,) + v)
                       for v in self.parent.bounds[i]]
                      for i in range(len(self.parent))]
            bboxes = [math.bounding_box(b) for b in bounds]
            order = sorted(list(range(len(self.parent))),
                           key=lambda i: bboxes[i].width)
        elif method == '1d sorted height':
            bounds = [[self.parent.to_pixel_vertex((i,) + v)
                       for v in self.parent.bounds[i]]
                      for i in range(len(self.parent))]
            bboxes = [math.bounding_box(b) for b in bounds]
            order = sorted(list(range(len(self.parent))),
                           key=lambda i: bboxes[i].width)
        else:
            raise ValueError('unrecognized pack method')

    # Move a component manually
    #   component:
    #       the component to move. Tuple.
    #   new_position:
    #       the new position. 2-element tuple if 2d, 1-element tuple or Real
    #       for 1d
    def move_component(self, component, new_position):
        if self.current_method[0:2] == '2d':
            # Check that position is valid for 2d
            if not isinstance(new_position, collections.Sequence):
                raise TypeError('for 2d packing positions must be sequences')
            elif len(new_position) != 2:
                raise ValueError('for 2d packing positions must be 2d')
            elif component not in self._positions:
                raise ValueError('unrecognized component \
                                 {0}'.format(component))
            elif not all([isinstance(c, numbers.Real) for c in new_position]):
                raise ValueError('coordinates must be numeric')
            else:
                # Reassign position
                self._positions[component] = new_position
        elif self.current_method[0:2] == '1d':
            # Check that position is valid for 1d
            if isinstance(new_position, collections.Sequence):
                if len(new_position) != 1:
                    raise ValueError('for 1d packing positions must be 1d')
                elif not isinstance(new_position[0], numbers.Real):
                    raise ValueError('coordinates must be numeric')
            elif not isinstance(new_position, numbers.Real):
                raise ValueError('coordinates must be numeric')
            else:
                # Reassign position
                self._positions[component] = new_position

    def __getitem__(self, key):
        return self._positions[key]


class TwoAndAHalfDGrid(grid.Grid, metaclass=ABCMeta):
    def __init__(self, altitude_map, *args, **kwargs):
        self._altitude_map = altitude_map
        super().__init__(*args, **kwargs)

    def manhattan(self, element_a, element_b, alt_a=None, alt_b=None):
        if self.element_type(element_b) != self.element_type(element_a):
            raise ValueError('elements are not the same type')
        if self.element_type(element_a) == 'face':
            alt_a = self.altitude(element_a) if alt_a is None else alt_a
            alt_b = self.altitude(element_b) if alt_b is None else alt_b
            manhattan_2d = super().manhattan(element_a, element_b)
            return manhattan_2d + abs(alt_a - alt_b)
        else:
            return super().manhattan(element_a, element_b)

    def euclidian(self, element_a, element_b, alt_a=None, alt_b=None):
        euclidian_2d = super().euclidian(element_a, element_b)
        return (euclidian_2d ** 2 + (alt_a - alt_b) ** 2) ** 0.5

    # Check whether an element is continuous—whether the adjacent quadrants
    # have the same altitude. This is always true for faces.
    #   element:
    #       the element to check for continuity. Tuple.
    # Returns True if continuous, False if not.
    def continuous(self, element):
        if self.element_type(element) == 'face':
            return True
        else:
            adj_faces = self.neighbors(element, 'face', True)
            altitudes = [self.altitude(face, d) for face, d in adj_faces]
            return set(altitudes) == 1

    # Get the altitude of an element or quadrant
    #   element:
    #       the element to get the altitude of. Tuple.
    #   face_quadrant:
    #       the quadrant of a face to get the altitude of. Only relevant for
    #       faces. None returns the altitude of the face center, identifier
    #       of an edge direction the altitude at the midpoint of the outer
    #       edge of that quadrant, and identifer of a vertex direction the
    #       altitude at the intersection of two quadrants.
    # Returns a numeric altitude or None, if not continuous.
    def altitude(self, element, face_quadrant=None):
        element_type = self.element_type(element)
        if element_type == 'face':
            # Altitude at the center of the face
            if face_quadrant is None:
                # Get the average of the quadrant altitudes
                return (sum(self._altitude_map[element]) /
                        len(self._altitude_map[element]))
            # Altitude at the midpoint of the outer edge of a quadrant
            elif face_quadrant in self.COORD_DIRS['edge']:
                # Get the quadrant altitude
                return self._altitude_map[element][face_quadrant]
            # Altitude at the intersection between two quadrants
            elif face_quadrant in self.COORD_DIRS['vertex']:
                # Find the two adjacent quadrants
                index = self.COORD_DIRS['vertex'].index(face_quadrant)
                adjacent = (self.COORD_DIRS['edge'][index],
                            self.COORD_DIRS['edge'][(index + 1) %
                                                    len(self.COORD_DIRS[
                                                            'edge'])])
                # Get the average of the adjacent quadrants' altitudes
                return sum([self.altitude(element, a) for a in adjacent]) / 2
        else:
            # Altitude undefined if not continuous
            if self.continuous(element):
                # Get adjacent quadrants' altitudes
                adj_faces = self.neighbors(element, 'face', True)
                altitudes = [self.altitude(face, d) for face, d in adj_faces]
                # Get the average of adjacent quadrants
                return sum(altitudes) / len(altitudes)
