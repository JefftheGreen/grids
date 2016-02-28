# /usr/bin/env python
# -*- coding: utf-8 -*-


import copy

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt

import new_math as math
from graphics.qt_classes import DraggableScene


class GridScene(DraggableScene):
    SCALE = 100

    DEFAULT_LABEL_PEN = QtGui.QPen(Qt.black, 0.01 * SCALE)
    DEFAULT_LABEL_PEN.setCosmetic(True)
    DEFAULT_LABEL_BRUSH = QtGui.QBrush(Qt.white)

    DEFAULT_FACE_PEN = QtGui.QPen()
    DEFAULT_FACE_PEN.setCosmetic(True)
    DEFAULT_FACE_BRUSH = QtGui.QBrush(Qt.white)

    DEFAULT_EDGE_PEN = QtGui.QPen(Qt.black, 0.03 * SCALE, cap=Qt.RoundCap)
    DEFAULT_EDGE_PEN.setCosmetic(True)

    DEFAULT_VERTEX_PEN = QtGui.QPen()
    DEFAULT_VERTEX_PEN.setCosmetic(True)
    DEFAULT_VERTEX_BRUSH = QtGui.QBrush(Qt.black)

    DEFAULT_OPTS = {'face': {'pen': DEFAULT_FACE_PEN,
                             'brush': DEFAULT_FACE_BRUSH,
                             'label text': True,
                             'label font': QtGui.QFont('Arial',
                                                       weight=100),
                             'label size': 33,
                             'label pen': DEFAULT_LABEL_PEN,
                             'label brush': DEFAULT_LABEL_BRUSH},
                    'edge': {'pen': DEFAULT_EDGE_PEN,
                             'label text': False,
                             'label font': QtGui.QFont('Arial',
                                                       weight=87),
                             'label size': 20,
                             'label pen': DEFAULT_LABEL_PEN,
                             'label brush': DEFAULT_LABEL_BRUSH},
                    'vertex': {'size': 25,
                               'pen': DEFAULT_VERTEX_PEN,
                               'brush': DEFAULT_VERTEX_BRUSH,
                               'label text': False,
                               'label font': QtGui.QFont('Arial',
                                                         weight=87),
                               'label size': 33,
                               'label pen': DEFAULT_LABEL_PEN,
                               'label brush': DEFAULT_LABEL_BRUSH}}


    LINK_PEN = None #QtGui.QPen(Qt.black, 10, Qt.SolidLine, Qt.FlatCap)

    ARROW_HEAD = QtGui.QPolygonF([QtCore.QPointF(0,0),
                                  QtCore.QPointF(40,15),
                                  QtCore.QPointF(0,30)])

    def __init__(self, parent, grid, rotate=False):
        super().__init__(parent)
        self.grid = grid
        self.draggable_items = []
        self.rotate = rotate
        self.grid_groups = []
        self.link_arcs = {}
        self.link_pens = {}
        self.grabbed_group = None
        self.hovering_group = None
        self.last_move_event_pos = None
        self.setBackgroundBrush(QtGui.QBrush(Qt.gray))
        for g in self.subgrids:
            grid_group = GridGroup(g,
                                   self.DEFAULT_OPTS,
                                   self.SCALE,
                                   self.grid,
                                   (self.subgrids.index(g)
                                    if self.subgrids[0] is not self.grid
                                    else None))
            grid_group.setZValue(len(self.grid_groups))
            self.grid_groups.append(grid_group)
            self.draggable_items.append(grid_group.outline)
            self.addItem(grid_group)
        if hasattr(self.grid, 'layout'):
            for g in range(len(self.subgrids)):
                x, y = self.grid.layout[g]
                x *= self.SCALE
                y *= -self.SCALE
                self.grid_groups[g].setPos(x + 600 * g, y - 100 * g)
        if hasattr(self.grid, 'links'):
            for link in self.grid.links:
                self.make_link_arc(link)

    @property
    def subgrids(self):
        if hasattr(self.grid, '_subgrids'):
            return self.grid._subgrids
        else:
            return [self.grid]

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Control:
            hovering_grid = None
            for gg in self.grid_groups:
                if gg.outline.contains(gg.mapFromScene(self.parent().mapToScene(self.parent().mapFromGlobal(QtGui.QCursor.pos())))):
                    gg.outline.setPen(QtGui.QPen(Qt.yellow, 30))
                    gg.setOpacity(1)
                    hovering_grid = gg
                    print(len(gg.faces) + len(gg.edges) + len(gg.vertices))
            for gg in self.grid_groups:
                if hovering_grid:
                    if gg != hovering_grid:
                        gg.setOpacity(0.25)
                        gg.outline.setPen(QtGui.QPen(Qt.transparent))
                else:
                    gg.setOpacity(1)
                    gg.outline.setPen(QtGui.QPen(Qt.transparent))

    def keyReleaseEvent(self, event):
        if event.key() == Qt.Key_Control:
            for gg in self.grid_groups:
                gg.setOpacity(1)
                gg.outline.setPen(QtGui.QPen(Qt.transparent))

    def mouseMoveEvent(self, event):
        super().mouseMoveEvent(event)
        key_modifiers = QtGui.QApplication.keyboardModifiers()
        if key_modifiers == QtCore.Qt.ControlModifier:
            hovering_grid = None
            for gg in self.grid_groups:
                if gg.outline.contains(gg.mapFromScene(event.scenePos())):
                    gg.outline.setPen(QtGui.QPen(Qt.yellow, 30))
                    gg.setOpacity(1)
                    hovering_grid = gg
            for gg in self.grid_groups:
                if hovering_grid:
                    if gg != hovering_grid:
                        gg.setOpacity(0.25)
                        gg.outline.setPen(QtGui.QPen(Qt.transparent))
                else:
                    gg.setOpacity(1)
                    gg.outline.setPen(QtGui.QPen(Qt.transparent))


    def on_drag_start(self, click_event, move_event):
        self.grabbed_group = None
        items = self.items(click_event.scenePos())
        for gg in self.grid_groups:
            if gg.outline in items:
                self.grabbed_group = gg
                break
        self.last_move_event_pos = (click_event.scenePos().x(),
                                    click_event.scenePos().y())
        self.on_drag(click_event, move_event)

    def on_drag(self, click_event, move_event):
        if self.grabbed_group:
            dx = self.last_move_event_pos[0] - move_event.scenePos().x()
            dy = self.last_move_event_pos[1] - move_event.scenePos().y()
            self.last_move_event_pos = (move_event.scenePos().x(),
                                        move_event.scenePos().y())
            self.move_subgrid(self.grid_groups.index(self.grabbed_group),
                              dx, dy)

    def on_drag_end(self, click_event, release_event):
        self.grabbed_group = None
        self.last_move_event_pos = None

    def move_subgrid(self, subgrid, dx, dy):
        x = self.grid_groups[subgrid].pos().x()
        y = self.grid_groups[subgrid].pos().y()
        self.grid_groups[subgrid].setPos(QtCore.QPointF(x - dx, y - dy))
        if len(self.link_arcs) > 0:
            links = self.grid.links[subgrid:]
            if links is not None:
                if self.grid.links[:subgrid] is not None:
                    links.update(self.grid.links[:subgrid])
                for link in links:
                    link_dirs = self.grid.links[link[0]:link[1]]
                    self.make_link_arc(tuple(zip(link, link_dirs)))
            else:
                if self.grid.links[:subgrid] is not None:
                    links = self.grid.links[:subgrid]
                    for link in links:
                        link_dirs = self.grid.links[link[0]:link[1]]
                        self.make_link_arc(tuple(zip(link, link_dirs)))

    def make_link_arc(self, link):
        (source, source_dir), (sink, sink_dir) = link
        source_subgrid = self.subgrids[source[0]]
        sink_subgrid = self.subgrids[sink[0]]
        if link in self.link_arcs:
            self.removeItem(self.link_arcs[link])
        source_neighbor = self.grid.neighbor(
            self.grid.neighbor(source, source_dir), source_dir)
        source_x, source_y = self.get_element_location(source)
        source_neighbor_x, source_neighbor_y = self.get_element_location(
            source_neighbor)
        sink_neighbor = self.grid.neighbor(self.grid.neighbor(sink, sink_dir),
                                           sink_dir)
        sink_x, sink_y = self.get_element_location(sink)
        sink_neighbor_x, sink_neighbor_y = self.get_element_location(
            sink_neighbor)
        start_point = QtCore.QLineF(source_x, source_y, source_neighbor_x, source_neighbor_y).pointAt(0.125)
        end_point = QtCore.QLineF(sink_x, sink_y, sink_neighbor_x, sink_neighbor_y).pointAt(0.125)
        end_x, end_y = end_point.x(), end_point.y()
        bezier = QtGui.QPainterPath(start_point)
        bezier.cubicTo(source_neighbor_x, source_neighbor_y,
                       sink_neighbor_x, sink_neighbor_y,
                       end_x, end_y)
        bezier_item = QtGui.QGraphicsPathItem(bezier)
        if link in self.link_pens:
            bezier_pen = self.link_pens[link]
        else:
            if self.LINK_PEN is None:
                hue = math.binary_divide_range(range(360),
                                               len(self.link_pens) + 1)
                color = QtGui.QColor.fromHsv(hue, 255, 255)
                bezier_pen = QtGui.QPen(color, 10, Qt.SolidLine, Qt.FlatCap)
            else:
                bezier_pen = self.LINK_PEN
            self.link_pens[link] = bezier_pen
        bezier_item.setPen(bezier_pen)
        head_item = QtGui.QGraphicsPolygonItem(self.ARROW_HEAD)
        rotation = sink_subgrid.DIR_ANGLES[sink_dir] + math.pi
        h_width, h_height = (head_item.boundingRect().width(),
                             head_item.boundingRect().height())
        head_item.setPos(end_x - 6*h_width/8, end_y - h_height/2)
        head_item.setTransformOriginPoint(6 * h_width/8, h_height/2)
        head_item.setRotation(rotation * 180 / math.pi)
        head_item.setBrush(QtGui.QBrush(bezier_pen.color()))
        head_item.setPen(QtGui.QPen(Qt.NoPen))
        link_item = QtGui.QGraphicsItemGroup()
        link_item.addToGroup(head_item)
        link_item.addToGroup(bezier_item)
        link_item.setZValue(len(self.grid_groups))
        self.link_arcs[link] = link_item
        self.addItem(link_item)

    def get_element_location(self, element):
        subgrid_id, subgrid_element = element[0], element[1:]
        subgrid = self.subgrids[subgrid_id]
        group = self.grid_groups[subgrid_id]
        x, y = None, None
        if subgrid.element_type(subgrid_element) == 'face':
            x, y = self.grid.to_pixel_face_center(element)
            x *= self.SCALE
            y *= -self.SCALE
            unadj_pos = QtCore.QPointF(x, y)
            pos = group.mapToScene(unadj_pos)
        elif subgrid.element_type(subgrid_element) == 'edge':
            # TODO: Write edge code
            pass
        elif subgrid.element_type(subgrid_element) == 'vertex':
            # TODO: Write vertex code
            pass
        return pos.x(), pos.y()


class GridGroup(QtGui.QGraphicsItemGroup):
    scale = 100

    def __init__(self, grid, default_opts, scale, parent_grid=None, index=None):
        super().__init__()
        self.grid = grid
        self.scale = scale
        self.default_opts = default_opts
        self.index = index
        self.parent_grid = grid if parent_grid is None else parent_grid
        self.label_queue = []
        self.faces, self.edges, self.vertices, self.labels = {}, {}, {}, {}
        self.element_group = QtGui.QGraphicsItemGroup()
        self.label_group = QtGui.QGraphicsItemGroup()
        self.addToGroup(self.element_group)
        self.addToGroup(self.label_group)
        self.setAcceptHoverEvents(True)
        self.setAcceptsHoverEvents(True)
        self.draw()
        self.outline = QtGui.QGraphicsPathItem(self.shape())
        outline_pen = QtGui.QPen(Qt.yellow, 20)
        outline_pen.setCosmetic(True)
        self.outline.setPen(QtGui.QPen(outline_pen))
        self.outline.setZValue(-1)
        blur = QtGui.QGraphicsBlurEffect()
        blur.setBlurHints(QtGui.QGraphicsBlurEffect.PerformanceHint)
        blur.setBlurRadius(15)
        self.outline.setGraphicsEffect(blur)
        self.addToGroup(self.outline)

    def draw(self):
        import time
        self.label_queue = []
        times = []
        for face in self.grid.all_faces:
            t = time.clock()
            self.draw_face(face)
            times.append(t-time.clock())
        for edge in self.grid.all_edges:
            t = time.clock()
            self.draw_edge(edge)
            times.append(t-time.clock())
        for vertex in self.grid.all_vertices:
            t = time.clock()
            self.draw_vertex(vertex)
            times.append(t-time.clock())
        for label in self.label_queue:
            t = time.clock()
            self.draw_label(*label)
            times.append(t-time.clock())

    def draw_face(self, face):
        points = [QtCore.QPointF(p[0] * self.scale, -p[1] * self.scale)
                  for p in self.grid.to_pixel_face_vertices(face)]
        poly = QtGui.QPolygonF(points)
        item = QtGui.QGraphicsPolygonItem(poly)
        proper_face = (face if self.index is None
                       else (self.index,) + face)
        face_properties = self.parent_grid.properties(proper_face, 'Qt')
        properties = copy.copy(self.default_opts['face'])
        print(self.index, face, properties['brush'].color().saturation())
        properties.update(face_properties)
        item.setPen(properties['pen'])
        item.setBrush(properties['brush'])
        if properties['label text'] is not False:
            lx, ly = self.grid.to_pixel_face_center(face)
            label_pos = lx * self.scale, -ly * self.scale
            self.label_queue.append((face, label_pos))
        self.faces[face] = item
        self.element_group.addToGroup(item)

    def draw_edge(self, edge):
        points = [QtCore.QPointF(p[0] * self.scale, -p[1] * self.scale)
                  for p in self.grid.to_pixel_edge_ends(edge)]
        line = QtCore.QLineF(*points)
        item = QtGui.QGraphicsLineItem(line)
        proper_edge = (edge if self.index is None
                       else (self.index,) + edge)
        edge_properties = self.parent_grid.properties(proper_edge, 'Qt')
        properties = copy.copy(self.default_opts['edge'])
        properties.update(edge_properties)
        item.setPen(properties['pen'])
        if properties['label text'] is not False:
            lx, ly = self.grid.to_pixel_edge_mid(edge)
            label_pos = lx * self.scale, -ly * self.scale
            self.label_queue.append((edge, label_pos,
                                     self.grid.edge_angle(
                                         edge) * 180 / math.pi))
        self.edges[edge] = item
        self.element_group.addToGroup(item)

    def draw_vertex(self, vertex):
        cx, cy = self.grid.to_pixel_vertex(vertex)
        cx *= self.scale
        cy *= -self.scale
        proper_vertex = (vertex if self.index is None
                         else (self.index,) + vertex)
        vertex_properties = self.parent_grid.properties(proper_vertex, 'Qt')
        properties = copy.copy(self.default_opts['vertex'])
        properties.update(vertex_properties)
        size = properties['size']
        item = QtGui.QGraphicsEllipseItem(cx - size / 2, cy - size / 2,
                                          size, size)
        item.setPen(properties['pen'])
        item.setBrush(properties['brush'])
        if properties['label text'] is not False:
            self.label_queue.append((vertex, (cx, cy)))
        self.vertices[vertex] = item
        self.element_group.addToGroup(item)

    def draw_label(self, element, location, angle=0):
        cx, cy = location
        element_type = self.grid.element_type(element)
        proper_element = (element if self.index is None
                          else (self.index,) + element)
        element_properties = self.parent_grid.properties(proper_element, 'Qt')
        properties = copy.copy(self.default_opts[element_type])
        properties.update(element_properties)
        if properties['label text'] in [None, True]:
            label = str(element).replace(' ', '').replace("'", '')
        else:
            label = str(properties['label text'])
        item = QtGui.QGraphicsSimpleTextItem(label)
        font = properties['label font']
        if properties['label size'] is not None:
            font.setPixelSize(properties['label size'])
            height = properties['label size']
        else:
            height = font.pointSize()
        item.setFont(font)
        item.setPen(properties['label pen'])
        item.setBrush(properties['label brush'])
        width = item.boundingRect().width()
        item.setTransformOriginPoint(width / 2, height / 2)
        item.setPos(cx - 0.5 * width, cy - 0.5 * height)
        item.setRotation(angle)
        item.translate(math.sin(-angle / 180 * math.pi) * height * 0.6,
                       math.cos(-angle / 180 * math.pi) * height * 0.6)
        self.labels[element] = item
        self.label_group.addToGroup(item)

    def shape(self):
        shape = QtGui.QPainterPath()
        for f in self.faces.values():
            shape.addPolygon(f.polygon())
        for v in self.vertices.values():
            shape += v.shape()
        shape = shape.simplified()
        return shape


def test():
    import grids.hexgrid as h
    import grids.special as s
    g0 = h.HexGrid(2)
    g0.holes = [(1, 0), (-1, 0)]
    g1 = h.HexGrid(0)
    g2 = h.HexGrid(0)
    tags = {}
    for face in g0.line_between((1,-2), (-1, 2), 'half ceiling'):
        tags[face] = True
    g0.add_tagset('line', 'face', tags)

    def brush(tags):
        if 'line' in tags:
            if tags['line'] == True:
                return QtGui.QBrush(Qt.red)
        return None

    map = {'brush': (brush, 1)}
    g0.add_property_map('Qt', 'line', map)
    links = {((0, -2, 0), (1, 0, 0)): (g0.NE, g1.SW),
                                ((0, -1,-1), (1,0,0)): (g0.N, g1.S),
                                ((0, -2, 1), (1,0, 0)): (g0.SE, g1.NW),
                                ((0,2,0), (1,0,0)): (g0.SW, g1.NE),
                                ((0,2,-1), (1,0,0)): (g0.NW, g1.SE),
                                ((0,1,1), (1,0,0)): (g0.S, g1.N)}
    g = s.LinkedGrid((g0, g1), links)
    g0.draw_qt(0.5)


if __name__ == '__main__':
    test()
