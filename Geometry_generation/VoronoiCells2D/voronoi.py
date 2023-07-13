#lines are stored in (point, slope) form, where slope==None for vertical
#polygons are stored as a clockwise list of vertices
#all polygons are assumed to be convex

from math import sqrt

def bounded_voronoi(bounding_polygon, point_set):
    """Returns a [point]->[polygon] hash, where:
       -the point is from point_set,
       -polygon is the points of point's voronoi polygon
    All points in point_set should be contained in bounding_polygon.
    """
    polygons = {}
    for point in point_set:
        bound = bounding_polygon
        for other_point in point_set:
            if other_point == point:
                continue
            hp = halfplane(point, other_point)
            bound = polygon_intersect_halfplane(bound, hp)
        polygons[point] = bound
    return polygons

def update_diagram(diagram, new_point, bounding_polygon):
    point_set = diagram.keys()
    changed = []
    for point in point_set:
        bound = diagram[point]
        hp = halfplane(point, new_point)
        new_bound = polygon_intersect_halfplane(bound, hp)
        if bound != new_bound:
            diagram[point] = new_bound
            changed.append(point)
    bound = bounding_polygon
    for other_point in point_set:
        hp = halfplane(new_point, other_point)
        bound = polygon_intersect_halfplane(bound, hp)
    diagram[new_point] = bound

def halfplane(inner_point, outer_point):
    """Returns a pair (line, {'lt'|'gt'}), where 'lt' means that inner_point is 
    below the line (or to the left for a vertical line), and the line is the 
    perpendicular bisector of the segment defined by the two given points.
    """
    mp = midpoint(inner_point, outer_point)
    s = slope(inner_point, outer_point)
    if s == None:
        new_slope = 0
    elif s == 0:
        new_slope = None
    else:
        new_slope = -(1/s)
    perpendicular_bisector = (mp, new_slope) 
    if new_slope == None:
        if inner_point[0] < mp[0]:
            half = 'lt'
        else:
            half = 'gt'
        return (perpendicular_bisector, half)
    inner_b = inner_point[1]-inner_point[0]*new_slope
    real_b = mp[1]-mp[0]*new_slope
    if inner_b < real_b:
        half = 'lt'
    else:
        half = 'gt'
    return (perpendicular_bisector, half)
    
def polygon_intersect_halfplane(polygon, halfplane):
    """Returns a new polygon which is the intersection of the two arguments.
    This is not general use: it assumes that there is a nonempty intersection.
    """
    n = len(polygon)
    point_locations = []
    #print "polygon_intersect_halfplane: polygon", polygon
    #print "polygon_intersect_halfplane: halfplane", halfplane
    for point in polygon:
        point_locations.append(point_in_halfplane(point, halfplane))
    if len(filter(bool, point_locations)) == n:
        return polygon
    while not (point_locations[0] and not point_locations[-1]):
        #print point_locations
        #print polygon
        point_locations = lshift_list(point_locations)
        polygon = lshift_list(polygon)
    #print "polygon_intersect_halfplane: polygon", polygon
    #print "polygon_intersect_halfplane: point_locations", point_locations
    new_polygon = []
    intersects = 0
    for p in range(n):
        #print "P=", p
        if point_locations[p]:
            new_polygon.append(polygon[p])
        elif intersects == 0:
            intersects = 1
            if point_locations[p] == None:
                new_polygon.append(polygon[p])
            else:
                segment = (polygon[p-1], polygon[p])
                line = halfplane[0]
                new_point = segment_intersect_line(segment, line)
                new_polygon.append(new_point)
        if p == n-1:
            #print "On the last point of the polygon"
            if point_locations[p] == None:
                #print " And it's a None"
                new_point = polygon[p]
                if new_point != new_polygon[-1]:
                    new_polygon.append(new_point)
                    #print "  And it was appended"
            else:
                #print " And it isn't a None"
                segment = (polygon[p], polygon[0])
                line = halfplane[0]
                new_point = segment_intersect_line(segment, line)
                if new_point == None:
                    print "ERROR"
                    print "Found no intersection between segment and line"
                    print segment, line
                new_polygon.append(new_point)
    #print "----returning from polygon_intersect halfplane---"
    return new_polygon

def lshift_list(l):
    end = l[0]
    l = l[1:]
    l.append(end)
    return l

def point_in_halfplane(point, halfplane):
    #print "point_in_halfplane: point, halfplane", point, halfplane
    #special case: vertical line
    if halfplane[0][1] == None:
        x = halfplane[0][0][0]
        if halfplane[1] == 'lt':
            if point[0] < x:
                return True
            elif point[0] == x:
                return None
            else:
                return False
        else:
            if point[0] > x:
                return True
            elif point[0] == x:
                return None
            else:
                return False
    #otherwise: not vertical...
    real_b = halfplane[0][0][1]-halfplane[0][0][0]*halfplane[0][1]
    point_b = point[1]-point[0]*halfplane[0][1]
    if halfplane[1] == 'lt':
        if point_b < real_b:
            return True
    else:
        if point_b > real_b:
            return True
    if real_b == point_b:
        return None
    return False

def segment_intersect_line(segment, line):
    #returns either a point (x, y) or None; endpoints not included
    seg_line = (segment[0], slope(segment[0], segment[1]))
    intersection = line_intersect(seg_line, line)
    if intersection == None:
        return None
    if point_on_segment(intersection, segment):
        return intersection
    else:
        return None

def line_intersect(l1, l2):
    # test le coefficient directeur
    if l1[1] == l2[1]:  
        return None
    ((x1, y1), m1) = l1
    ((x2, y2), m2) = l2
    if m1 == None:
        x = x1
        y = m2*(x-x2)+y2
        return (x, y)
    if m2 == None:
        x = x2
        y = m1*(x-x1)+y1
        return (x, y)
    x = (m1*x1-y1-m2*x2+y2)/(m1-m2)
    y = m1*(x-x1)+y1
    return (x, y)


def point_on_segment(point, segment):
    #IMPORTANT: the point is assumed to be on the line defined by the segment
    d1 = dist(segment[0], segment[1])
    d2 = dist(point, segment[0])
    d3 = dist(point, segment[1])
    return ((d3 <= d1) and (d2 <= d1))

def dist(p1, p2):
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def midpoint(p1, p2):
    return ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)

def normalize(v):
    length = sqrt(v[0]**2 + v[1]**2)
    return (v[0]/length, v[1]/length)

def slope(p1, p2):
    if (p2[0] == p1[0]):
        #special case: vertical line
        return None
    return (p2[1]-p1[1])/(p2[0]-p1[0])
