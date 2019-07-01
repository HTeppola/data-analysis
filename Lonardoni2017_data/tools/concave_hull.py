from scipy.spatial import Delaunay
import numpy as np
import networkx as nx
import pylab as plt
import math

alpha = .1


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def hull(alpha, points, plot=False):
    edges = set()
    edge_points = []
    edge_2_remove = set()
    out_2_remove = set()

    points = np.array(points)

    tri = Delaunay(points)

    def add_edge(i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            edge_2_remove.add((i, j))
            return
        edges.add((i, j))
        edge_points.append(points[[i, j]])

    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        # Semiperimeter of triangle
        s = (a + b + c) / 2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(ia, ib)
            add_edge(ib, ic)
            add_edge(ic, ia)
        else:
            out_2_remove.add((ia, ib))
            out_2_remove.add((ib, ic))
            out_2_remove.add((ic, ia))

    K = nx.Graph()

    contour = [edge for edge in edges if edge not in edge_2_remove]
    contour = [edge for edge in contour if (edge[1], edge[0]) not in edge_2_remove]
    out = [edge for edge in out_2_remove if edge not in contour]
    out = [edge for edge in out if (edge[1], edge[0]) not in contour]

    K.add_edges_from(contour)
    for n in K.nodes(): K.node[n]['pos'] = points[n]
    if plot:
        plt.figure(500)
        ax = plt.gca()
        nx.draw_networkx_edges(K, pos=nx.get_node_attributes(K, 'pos'), ax=ax)

        x, y = list(zip(*points))
        plt.scatter(x, y)
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        simpleaxis(plt.gca())

    # plt.show(block=True)
    rpos = []
    ipos = []

    nodelist = K.nodes()

    if len(nodelist) <= 3:
        tri_trlist = []
        tri_pos = tri.points
        tri_list2 = tri.vertices
        for i in range(len(tri_list2)):
            tri_trlist.append([tri_pos[tri_list2[i][0]], tri_pos[tri_list2[i][1]], tri_pos[tri_list2[i][2]]])
            ipos.append(i)
            print("TRIANGOLO!")
        K.add_edges_from(edges)
        for n in K.nodes(): K.node[n]['pos'] = points[n]
        if plot:
            plt.figure(500)
            ax = plt.gca()
            nx.draw_networkx_edges(K, pos=nx.get_node_attributes(K, 'pos'), ax=ax)
            x, y = list(zip(*points))
            plt.scatter(x, y)
            plt.xlim([0, 1])
            plt.ylim([0, 1])
            simpleaxis(plt.gca())

    findex = [nodes for nodes in K.nodes() if len(K.edges(nodes)) > 2 for _ in range(len(K.edges(nodes)) / 2)]
    if len(findex) == 0:
        if len(nodelist):
            findex.append(nodelist[0])

    elen = 0
    pivot = findex[elen]
    rpos.append(points[pivot])
    ipos.append(pivot)
    for i in range(len(contour) - 1):

        for edge in K.edges(pivot):
            if K.degree(edge[1]) % 2 == 0:  # edge[1] not in ipos and
                pivot = edge[1]
                rpos.append(points[pivot])
                ipos.append(pivot)
                if i > 0:
                    K.remove_edge(edge[0], edge[1])
                break
                #        K.remove_edge(edge[1],edge[0])
                #        if not flag:
                #            print "something is working bad: MAYBE DOUGHNUT",pivot,elen,len(findex)
                #            try:
                #                nodelist.append(findex[elen])
                #                elen+=1
                #                pivot=findex[elen]
                #                rpos.append(points[pivot])
                #                ipos.append(pivot)
                #            except:
                #                break

    tri2 = Delaunay(rpos)
    tri_pos = tri2.points
    tri_list2 = tri2.vertices

    tri_list = []
    for i in range(len(tri_list2)):
        ele = [0, 0, 0]
        for j in range(3):
            ele[j] = ipos[tri_list2[i][j]]
        # print ele
        tri_list.append(ele)
        # print tri_list

    tri_trlist = []
    for i in range(len(tri_list)):
        ele = tri_list[i]
        flag = 0
        #        for edge in edge_2_remove:
        #            if edge[0] in ele and edge[1] in ele:
        #                flag=1
        # print ele,out
        for edge in out:

            if edge[0] in ele and edge[1] in ele:
                flag = 1
                # print "destroyed",edge[0],edge[1],"--->", ele

        if flag == 0:
            tri_trlist.append([tri_pos[tri_list2[i][0]], tri_pos[tri_list2[i][1]], tri_pos[tri_list2[i][2]]])
            #    for ele in tri_trlist:
            #        plt.plot([ele[i%3][0] for i in xrange(len(ele)+1)],[ele[i%3][1] for i in xrange(len(ele)+1)],'g',lw=5)
            #    if len(ipos)==0:
            #            print ipos,points,tri.vertices,nodelist
            #            time.sleep(10000)
    return tri_trlist, ipos