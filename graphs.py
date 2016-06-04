# /usr/bin/env python
# -*- coding: utf-8 -*-

import igraph


def make_igraph(grid, cost_property, passable_property,
                vertex_element='face', directed=True):
    if vertex_element == 'face':
        all_vertex_elements = grid.all_faces
    elif vertex_element == 'edge':
        all_vertex_elements = grid.all_edges
    elif vertex_element == 'vertex':
        all_vertex_elements = grid.all_vertices
    else:
        raise ValueError('vertex_element not recognized')
    graph = igraph.Graph(len(all_vertex_elements), directed=directed)
    graph.vs['coordinates'] = list(all_vertex_elements)
    for vertex in all_vertex_elements:
        adjacent = grid.neighbors(vertex, vertex_element)
        for adj_vertex in [val for val in adjacent
                           if val in all_vertex_elements]:
            source = graph.vs.select(coordinates_eq=vertex)[0]
            sink = graph.vs.select(coordinates_eq=adj_vertex)[0]
            link = adjacent[adj_vertex]
            cost = grid.properties(link, 'cost')[cost_property]
            passable = grid.properties(link, 'passable')[passable_property]
            if passable:
                graph.add_edge(source, sink, coordinates=link, cost=cost)
    if not directed:
        graph.simplify()
    return graph


def test():
    from grids import hexgrid
    h = hexgrid.HexGrid(5)
    graph = make_igraph(h, 'cost', 'passable', 'face', False)
    print(*list(h.neighbors((0,0,h.E), 'vertex').keys()))
    import cairo
    igraph.drawing.plot(graph, 'test.png')

if __name__ == '__main__':
    test()