import sys
import os
from collections import defaultdict
from bisect import bisect_left
import subprocess
import toolshed as ts

from interlap import InterLap, Interval as IntervalSet, reduce as ireduce
import numpy as np


def overlaps(s1, e1, s2, e2):
    """
    >>> overlaps(2, 4, 3, 5)
    True
    >>> overlaps(2, 4, 4, 5)
    False
    >>> overlaps(2, 200, 3, 5)
    True
    >>> overlaps(3, 5, 2, 200)
    True
    >>> overlaps (3, 5, 5, 6)
    False
    >>> overlaps (5, 6, 3, 5)
    False
    """
    return not (e1 <= s2 or s1 >= e2)

def split_ranges(ranges, splitters, varflags): # if range is in splitters, it is removed from potential constraint regions; import my version of interlap
    """
    >>> split_ranges([(1018, 1034)], [(1022, 1034)], ['VARFALSE'])
    ([[(1018, 1022)]], [['VARFALSE']])

    >>> split_ranges([(1018, 1034)], [(1030, 1032)], ['VARFALSE'])
    ([[(1018, 1030)], [(1032, 1034)]], [['VARFALSE'], ['VARFALSE']])

    >>> split_ranges([(1018, 1034)], None, ['VARFALSE'])
    ([[(1018, 1034)]], [['VARFALSE']])

    >>> split_ranges([(1018, 1034), (1045, 1069)], None, ['VARFALSE', 'VARTRUE'])
    ([[(1018, 1034), (1045, 1069)]], [['VARFALSE', 'VARTRUE']])

    >>> split_ranges([(1018, 1034), (1045, 1069)], [(1030, 1032), (1047, 1050)], ['VARFALSE', 'VARTRUE'])
    ([[(1018, 1030)], [(1032, 1034)], [(1045, 1047)], [(1050, 1069)]], [['VARFALSE'], ['VARFALSE'], ['VARTRUE'], ['VARTRUE']])

    >>> split_ranges([(1018, 1034), (1045, 1069)], [(1030, 1050)], ['VARFALSE', 'VARTRUE'])
    ([[(1018, 1030)], [(1050, 1069)]], [['VARFALSE'], ['VARTRUE']])

    >>> split_ranges([(1018, 1034)], [(1022, 1024), (1028, 1034)], ['VARFALSE'])
    ([[(1018, 1022)], [(1024, 1028)]], [['VARFALSE'], ['VARFALSE']])

    >>> split_ranges([(18, 24), (28, 35), (55, 60)], [(28, 35), (55, 57)], ['VARFALSE', 'VARTRUE', 'VARTRUE'])
    ([[(18, 24)], [(57, 60)]], [['VARFALSE'], ['VARTRUE']])

    >>> split_ranges([(12, 18), (22, 28), (32, 39), (42, 48)],
    ...                 [(12, 18), (22, 26),           (44, 48)], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])
    ([[(26, 28)], [(32, 39)], [(42, 44)]], [['VARFALSE'], ['VARFALSE'], ['VARTRUE']])

    >>> split_ranges([(11, 18), (22, 28), (32, 39), (42, 48)],
    ...                 [(12, 18), (22, 26),           (44, 48)], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])
    ([[(11, 12)], [(26, 28)], [(32, 39)], [(42, 44)]], [['VARFALSE'], ['VARFALSE'], ['VARFALSE'], ['VARTRUE']])

    >>> split_ranges([(38826782, 38826890), (38827874, 38828144), (38828232, 38828286), (38834219, 38834405), (38834632, 38834759), (38834935, 38835008), (38837089, 38837266), (38842900, 38843135), (38845349, 38845487), (38847123, 38847183), (38847381, 38847501), (38848917, 38848968), (38849071, 38849201), (38850109, 38850221), (38851128, 38851287), (38851370, 38851483), (38852287, 38852501), (38852849, 38852912), (38853015, 38853214), (38853353, 38853472), (38855530, 38855611), (38855700, 38855757), (38857795, 38857952), (38858148, 38858212), (38858320, 38858416), (38858687, 38858777), (38860611, 38860612)],
    ...             [(38826782, 38826890), (38827874, 38828144), (38828232, 38828286), (38834219, 38834405), (38834632, 38834759), (38834935, 38835008), (38837089, 38837266), (38842900, 38843135)],
    ...             [ 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])
    ([[(38845349, 38845487), (38847123, 38847183), (38847381, 38847501), (38848917, 38848968), (38849071, 38849201), (38850109, 38850221), (38851128, 38851287), (38851370, 38851483), (38852287, 38852501), (38852849, 38852912), (38853015, 38853214), (38853353, 38853472), (38855530, 38855611), (38855700, 38855757), (38857795, 38857952), (38858148, 38858212), (38858320, 38858416), (38858687, 38858777), (38860611, 38860612)]], [['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE']])

    >>> split_ranges([(26782, 26890), (27874, 28144), (28232, 28286), (34219, 34405), (34632, 34759), (34935, 35008), (37089, 37266), (42900, 43135), (45349, 45487), (47123, 47183), (47381, 47501), (48917, 48968), (49071, 49201), (50109, 50221), (51128, 51287), (51370, 51483), (52287, 52501), (52849, 52912), (53015, 53214), (53353, 53472), (55530, 55611), (55700, 55757), (57795, 57952), (58148, 58212), (58320, 58416), (58687, 58777), (60611, 60612)],
    ...             [(27874, 28144), (28232, 28286), (34219, 34405), (34632, 34759), (34935, 35008), (37089, 37266), (42900, 43135)],
    ...             [ 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])
    ([[(26782, 26890)], [(45349, 45487), (47123, 47183), (47381, 47501), (48917, 48968), (49071, 49201), (50109, 50221), (51128, 51287), (51370, 51483), (52287, 52501), (52849, 52912), (53015, 53214), (53353, 53472), (55530, 55611), (55700, 55757), (57795, 57952), (58148, 58212), (58320, 58416), (58687, 58777), (60611, 60612)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE']])

    >>> split_ranges([(3012265, 3012435), (3012435, 3012437)], [(3012265, 3012442), (3012442, 3012445)], ['VARFALSE', 'VARTRUE'])
    ([], [])
    """
    results=[]
    if splitters is None:
        return [ranges], [varflags]
    results=[x._vals for x in IntervalSet(ranges).split(splitters)]
    vf=[]; res=[]
    for i, ivs in enumerate(results): # res
        v=[]
        for j, iv, in enumerate(ivs):
            for k, r in enumerate(ranges):
                if overlaps(iv[0], iv[1], r[0], r[1]):
                    v.append(varflags[k])
                    break
        vf.append(v)
    assert len(results) == len(vf)
    res = [res]; vf = vf
    return results, vf

def get_ranges(last, vstart, vend, exon_starts, exon_ends, chrom=1): # NOTE: new model version
    """
    >>> get_ranges(874772, 874778, 874827, [874655], [874827])
    ([(874772, 874777), (874777, 874827)], 874827, ['VARFALSE', 'VARTRUE'])
    
    >>> get_ranges(61018, 62029, 62029, (
    ... 60174, 60370, 60665, 60925, 62029, 62216, 62453,
    ... 62675, 63052, 63398, 63652, 63868, 64512, 64764,
    ... 65018, 65671), (60281,
    ... 60565, 60808, 61033, 62134, 62379, 62587, 62824,
    ... 63209, 63559, 63779, 64102, 64691, 64946, 65084,
    ... 65985))
    ([(61018, 61033)], 61018, ['VARFALSE'])

    >>> get_ranges(350, 346, 369, (350, 350), (400, 400))
    ([], 369, [])
    
    >>> get_ranges(0, 1, 1, (0, 20), (10, 30)) # needs varflag
    ([(0, 1)], 0, ['VARTRUE'])

    >>> get_ranges(0, 3, 3, (0, 20), (10, 30)) # second region needs varflag
    ([(0, 2), (2, 3)], 0, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(3, 6, 6, (0, 20), (10, 30)) # second region needs varflag
    ([(3, 5), (5, 6)], 3, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(6, 7, 7, (0, 20), (10, 30)) # needs varflag
    ([(6, 7)], 6, ['VARTRUE'])

    >>> get_ranges(7, 9, 9, (0, 20), (10, 30)) # needs varflag for second region only
    ([(7, 8), (8, 9)], 7, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(65400, 65404, 65425, (65242, 65242), (65601, 65601))
    ([(65400, 65403), (65403, 65425)], 65425, ['VARFALSE', 'VARTRUE'])
    
    >>> get_ranges(65405, 65404, 65425, (65242, 65242), (65601, 65601))
    ([], 65425, [])

    >>> get_ranges(61018, 61990, 62001, (60925, 62000), (61033, 62040))
    ([(61018, 61033), (62000, 62001)], 62001, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(61018, 62023, 62030, (60925, 62000), (61033, 62040)) # this situation in the future could possibly become (62022, 62030) with the vartrue flag
    ([(61018, 61033), (62000, 62022), (62022, 62030)], 62030, ['VARFALSE', 'VARFALSE', 'VARTRUE'])
    
    >>> get_ranges(62000, 62023, 62050, (60925, 62000, 62045), (61033, 62040, 62060)) # this situation in the future could possibly become (62022, 62030) with the vartrue flag
    ([(62000, 62022), (62022, 62040), (62045, 62050)], 62050, ['VARFALSE', 'VARTRUE', 'VARTRUE'])

    >>> get_ranges(62000, 62023, 62070, (60925, 62000, 62045), (61033, 62040, 62060)) # this situation in the future could possibly become (62022, 62030) with the vartrue flag
    ([(62000, 62022), (62022, 62040), (62045, 62060)], 62070, ['VARFALSE', 'VARTRUE', 'VARTRUE'])

    >>> get_ranges(56, 95, 95, range(0, 1000, 10), range(5, 1000, 10))
    ([(60, 65), (70, 75), (80, 85), (90, 94), (94, 95)], 60, ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])

    >>> get_ranges(0, 10, 10, range(0, 100, 10), range(5, 100, 10))
    ([(0, 5)], 0, ['VARFALSE'])

    >>> get_ranges(98320, 98324, 98324, [14537, 14540, 15243, 15812, 16718, 17156, 17398, 18462, 18739, 19230, 19885, 24207, 26058, 27040, 28168, 28434, 28900, 29938, 32062, 34212, 34940, 36058, 36850, 38476, 48595, 49413, 49835, 50111, 52360, 53750, 64266, 64724, 70062, 74665, 82034, 83129, 84984, 91647, 96432, 98175], [14540, 14675, 15306, 15946, 16822, 17265, 17477, 18605, 18860, 19344, 19988, 24300, 26216, 27110, 28335, 28676, 29000, 30064, 32200, 34293, 35025, 36157, 36886, 38512, 48701, 49489, 49952, 50217, 52506, 53791, 64429, 64859, 70218, 74830, 82223, 83372, 85151, 91769, 96574, 98323])
    ([(98320, 98323)], 98323, ['VARFALSE'])

    >>> get_ranges(74667, 98324, 98324, [14537, 14540, 15243, 15812, 16718, 17156, 17398, 18462, 18739, 19230, 19885, 24207, 26058, 27040, 28168, 28434, 28900, 29938, 32062, 34212, 34940, 36058, 36850, 38476, 48595, 49413, 49835, 50111, 52360, 53750, 64266, 64724, 70062, 74665, 82034, 83129, 84984, 91647, 96432, 98175], [14540, 14675, 15306, 15946, 16822, 17265, 17477, 18605, 18860, 19344, 19988, 24300, 26216, 27110, 28335, 28676, 29000, 30064, 32200, 34293, 35025, 36157, 36886, 38512, 48701, 49489, 49952, 50217, 52506, 53791, 64429, 64859, 70218, 74830, 82223, 83372, 85151, 91769, 96574, 98323])
    ([(74667, 74830), (82034, 82223), (83129, 83372), (84984, 85151), (91647, 91769), (96432, 96574), (98175, 98323)], 98323, ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE'])

    >>> get_ranges(4980013, 4980017, 4980019, [4977206, 4977209, 4978695, 4979932, 4980230, 4982744, 4983938, 4985117, 4987882, 4988407, 4989652, 4992118, 4993005, 4994330, 4995278, 4998358, 5001052, 5017596, 5062631], [4977209, 4977318, 4978747, 4980018, 4980239, 4982768, 4984022, 4985255, 4987955, 4988492, 4989788, 4992186, 4993065, 4994527, 4995323, 4998419, 5001083, 5017601, 5062651])
    ([(4980013, 4980016), (4980016, 4980018)], 4980019, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(933226, 933256, 933260, [896196, 896921, 899238, 899515, 904127, 907996, 909967, 914189, 915861, 916586, 921957, 930626, 931761, 933129], [896475, 897023, 899363, 899669, 904254, 908184, 910121, 914257, 915954, 916694, 922073, 930789, 931898, 933259])
    ([(933226, 933255), (933255, 933259)], 933260, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(62000, 62023, 62070, (60925, 62000, 62045, 62080), (61033, 62040, 62060, 62100))
    ([(62000, 62022), (62022, 62040), (62045, 62060)], 62070, ['VARFALSE', 'VARTRUE', 'VARTRUE'])

    >>> get_ranges(74300, 75689, 75703, [69160, 69613, 70817, 71170, 71437, 72279, 72937, 73246, 73698, 74143, 74298, 75693], [69407, 70063, 70928, 71312, 71630, 72420, 73063, 73483, 73804, 74298, 74301, 75699])
    ([(74300, 74301), (75693, 75699)], 75703, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(56832, 57012, 57080, [56792, 57019, 57447, 58103, 59303, 59430, 59682, 60936, 61007], [56846, 57073, 57501, 58175, 59339, 59497, 59850, 61007, 61010])
    ([(56832, 56846), (57019, 57073)], 57080, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(20756, 20888, 20942, [20569, 20896, 20899, 21009, 21573, 22326, 22633, 22851, 23403, 26831, 40040], [20771, 20899, 20905, 21230, 21688, 22548, 22768, 23039, 23835, 26921, 40131])
    ([(20756, 20771), (20896, 20899), (20899, 20905)], 20942, ['VARFALSE', 'VARTRUE', 'VARTRUE'])
    """

    varflag=[]
 
    assert last >= exon_starts[0]
    #assert vstart <= exon_ends[-1]
    assert all(s < e for s, e in zip(exon_starts, exon_ends))

    istart = bisect_left(exon_starts, last) - 1
    if istart == -1: istart = 0
    [(0, 6), (10, 11)]

    assert exon_starts[istart] <= last, (exon_starts[istart], last, istart)
    if exon_ends[istart] <= last:
        istart += 1
        if istart < len(exon_starts): # in case last and exon ends are equal, it will loop through again, but I want that last loop
            last = exon_starts[istart]
    start = last

    if vstart < vend: # moved here because there are variants in UTRs that do not exist in coding exon space
        last = vend
    ranges = []
    while start < vstart and istart < len(exon_starts): #<= lets it capture 0 length regions, so I removed it and the +1 allows it to make 1 bp regions when two variants are right next to one another
        ranges.append((start, exon_ends[istart])) #removed +1 from exon_ends[istart] + 1, because IntervalSet is already in 0-based half-open format
        istart += 1
        varflag.append("VARFALSE") # unless using vstart, there is no variant
        try: 
            if exon_starts[istart] > vstart and exon_starts[istart] < vend and ranges[-1][1] < vstart:
                varflag.append("VARTRUE")
                if vend < exon_ends[-1] and vend > exon_starts[-1]:
                    ranges.append((exon_starts[-1], vend))
                    break
                elif vend > exon_ends[-1] or (vend < exon_starts[istart+1] and vend > exon_ends[istart]):
                    ranges.append((exon_starts[istart], exon_ends[istart]))  
                    break
                elif vend < exon_ends[istart]:
                    ranges.append((exon_starts[istart], vend))
                    break
                else:
                    while vend > exon_starts[istart]:
                        ranges.append((exon_starts[istart], exon_ends[istart]))
                        istart+=1
                    varflag.append("VARTRUE")
                    break
        except IndexError:
            pass
        if ranges[-1][1] >= vstart: # equal to is now possible, since we are including variant start+1 and ranges are in 0-based half-open
            #print vstart, vend, exon_ends[istart], ranges[-1][1]
            if ranges[-1][0]-(vstart-1)==0:
                ranges[-1] = (ranges[-1][0], vstart)
                varflag[-1]="VARTRUE"
                break
            ranges[-1] = (ranges[-1][0], vstart-1) #removed +1 from vstart + 1, because IntervalSet is already in 0-based half-open format
            varflag.append("VARTRUE") #variant contained at end coordinate = TRUE; this indicates the region contains the variant, therefore should be considered 0 bp, get a 0 coverage and a 0 cpg score
            if exon_ends[istart-1] < vend:
                if exon_ends[-1] < vend:
                    #varflag.append("VARTRUE")
                    ranges.append((vstart-1, exon_ends[istart-1]))
                    if exon_starts[-1] < vend:
                        if ranges[-1][1] < exon_starts[-1]:
                            varflag.append("VARTRUE")
                            ranges.append((exon_starts[-1], exon_ends[-1]))
                    break
                ranges.append((vstart-1, exon_ends[istart-1]))
                if vend < exon_starts[istart]:
                    ranges[-1]=(ranges[-1][0], exon_ends[istart-1])
                    break
                if vend < exon_ends[istart]:
                    ranges.append((exon_starts[istart], vend))
                    varflag.append("VARTRUE")
                    break
                else:
                    ranges.append((exon_starts[istart], exon_ends[istart]))
                    varflag.append("VARTRUE")
                    break
            ranges.append((vstart-1, vend))
            break
        if exon_ends[-1]==exon_ends[istart-1] and vstart > exon_ends[-1]:
            last=exon_ends[-1]
            break
        start = exon_starts[istart]

    return ranges, last, varflag

def get_ranges_w_variant(last, vstart, vend, exon_starts, exon_ends, chrom=1): # NOTE: the previous version of the model where end coordinate contains variant
    """
    >>> get_ranges_w_variant(61018, 62029, 62029, (
    ... 60174, 60370, 60665, 60925, 62029, 62216, 62453,
    ... 62675, 63052, 63398, 63652, 63868, 64512, 64764,
    ... 65018, 65671), (60281,
    ... 60565, 60808, 61033, 62134, 62379, 62587, 62824,
    ... 63209, 63559, 63779, 64102, 64691, 64946, 65084,
    ... 65985))
    ([(61018, 61033)], 61018, ['VARFALSE'])

    >>> get_ranges_w_variant(38865400, 38865404, 38865425, (38865242, 38865242), (38865601, 38865601))
    ([(38865400, 38865404)], 38865425, ['VARTRUE'])
    
    >>> get_ranges_w_variant(38865405, 38865404, 38865425, (38865242, 38865242), (38865601, 38865601))
    ([], 38865425, [])

    >>> get_ranges_w_variant(617350, 617346, 617369, (617350, 617350), (617400, 617400))
    ([], 617369, [])

    >>> get_ranges_w_variant(61018, 61990, 62001, (60925, 62000), (61033, 62040))
    ([(61018, 61033)], 62001, ['VARFALSE'])

    >>> get_ranges_w_variant(61018, 62023, 62030, (60925, 62000), (61033, 62040)) #varflag for last one only
    ([(61018, 61033), (62000, 62023)], 62030, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges_w_variant(56, 95, 95, range(0, 1000, 10), range(5, 1000, 10)) #varflag for last one only
    ([(60, 65), (70, 75), (80, 85), (90, 95)], 60, ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])

    >>> get_ranges_w_variant(1, 10, 10, range(0, 100, 10), range(5, 100, 10))
    ([(1, 5)], 1, ['VARFALSE'])

    >>> get_ranges_w_variant(0, 10, 10, range(0, 100, 10), range(5, 100, 10))
    ([(0, 5)], 0, ['VARFALSE'])

    >>> get_ranges_w_variant(50, 59, 59, (50, 61), (60, 70))
    ([(50, 59)], 50, ['VARTRUE'])

    >>> get_ranges_w_variant(1562576, 1562675, 1562675,
    ... (1560665, 1560925, 1562029, 1562216, 1562453, 1562675, 1563052, 1563398, 1563652, 1563868, 1564512, 1564764, 1565018, 1565671),
    ... (1560808, 1561033, 1562134, 1562379, 1562587, 1562824, 1563209, 1563559, 1563779, 1564102, 1564691, 1564946, 1565084, 1565985))
    ([(1562576, 1562587)], 1562576, ['VARFALSE'])
    """

    #f4=open('deletioncut.txt','a') #code removed by deletions

    varflag=[]
 
    assert last >= exon_starts[0]
    assert vstart <= exon_ends[-1]
    #assert vstart >= last, (vstart, last, exon_starts) # deletion can overlap UTR/intron, but this is controlled for in exac-regions.py
    assert all(s < e for s, e in zip(exon_starts, exon_ends))

    istart = bisect_left(exon_starts, last) - 1
    if istart == -1: istart = 0
    [(0, 6), (10, 11)]

    assert exon_starts[istart] <= last, (exon_starts[istart], last, istart)
    if exon_ends[istart] <= last:
        istart += 1
        if istart < len(exon_starts): # in case last and exon ends are equal, it will loop through again, but I want that last loop
            last = exon_starts[istart]
    start = last

    if vstart < vend: # moved here because there are variants in UTRs that do not exist in coding exon space
        last = vend
    ranges = []; writedel=True
    while start < vstart and istart < len(exon_starts): #<= lets it capture 0 length regions, so I removed it and the +1 allows it to make 1 bp regions when two variants are right next to one another
        ranges.append((start, exon_ends[istart])) #removed +1 from exon_ends[istart] + 1, because IntervalSet is already in 0-based half-open format
        #if writedel and vstart < vend and vend > exon_starts[istart]:
        #    f4.write("\t".join(map(str,[chrom,vstart,last]))+"\n") # removed by deletion, but only if it is within an exon
        #    writedel=False
        istart += 1
        varflag.append("VARFALSE") # unless using vstart, there is no variant
        if ranges[-1][1] >= vstart: # equal to is now possible, since we are including variant start+1 and ranges are in 0-based half-open
            ranges[-1] = (ranges[-1][0], vstart) #removed +1 from vstart + 1, because IntervalSet is already in 0-based half-open format
            varflag[-1]="VARTRUE" #variant contained at end coordinate = TRUE
            break
        start = exon_starts[istart]

    return ranges, last, varflag

import doctest
res = doctest.testmod()
if res.failed != 0:
    sys.stderr.write("FAILING TESTS")
    sys.exit(1)


def path(p):
    return os.path.expanduser(os.path.expandvars(p))

def floatfmt(v, prec="%.2f"):
    return (prec % v).rstrip('0').rstrip('.')

def read_coverage(chrom, cov=10, length=249250621, path="data/exacv2.chr{chrom}.cov.txt.gz"): #length may need to be fixed in the future, if new chromosome lengths are established in GRCh38 #or data/Panel.chr{chrom}.coverage.txt.gz for exacv1
    """
    read ExAC coverage from a single chrom into a numpy array. If no length is
    given, just use the one length from chrom 1.
    path is expected to contain Panel.chr*
    cov is the column to pull
    """

    cols = "chrom	pos	mean	median	1	5	10	15	20	25	30	50 100".split()
    coli = cols.index(str(cov)) + 1
    path=path.format(**locals()) 
    # just extract the position (2) and the requested column
    p = subprocess.Popen("tabix {path} {chrom} | cut -f 2,{coli} ".format(**locals()), # {path}/exacv2.chr{chrom}.cov.txt.gz #{path}/Panel.chr{chrom}.coverage.txt.gz
            stdout=subprocess.PIPE, stderr=sys.stderr,
            shell=True,
            executable=os.environ.get("SHELL"))

    cov = np.zeros(length, dtype=np.float32)
    j = 0
    for line in p.stdout:
        pos, val = line.split()
        cov[int(pos)-1] = float(val)
        j += 1
        #if j > 100000: break
    assert j > 0, ("no values found for", chrom, path)
    p.wait()
    if p.returncode != 0:
        raise Exception("bad: %d", p.returncode)
    return cov


def read_exons(gtf, chrom, cutoff, coverage_array, exclude):
    genes = defaultdict(IntervalSet)
    splitters = defaultdict(IntervalSet)


    interlaps = []
    split_iv = InterLap()
    # preempt any bugs by checking that we are getting a particular chrom
    assert gtf[0] == "|", ("expecting a tabix query so we can handle chroms correctly")
    #f1 = open("selfchaincut.txt","a")
    #f2 = open("segdupscut.txt","a")
    #f3 = open("coveragecut.txt","a")
    for bed in exclude:
        # expecting a tabix query so we can handle chroms correctly
        a = "|tabix {bed} {chrom}".format(chrom=chrom, bed=bed)
    
        # any file that gets sent in will be used to split regions (just like
        # low-coverage). For example, we split on self-chains as well.
#TODO: comment this block if you don't want any filtering by self-chains or segdups
        for toks in (x.strip().split("\t") for x in ts.nopen(a)): # adds self chains and segdups to splitters list, so that exons can be split, and they are removed from CCRs
            s, e = int(toks[1]), int(toks[2])
            split_iv.add((s, e))
            #if len(toks) > 3:
            #    f1.write("\t".join(toks)+"\n") # self chain
            #else:
            #    f2.write("\t".join(toks)+"\n") # segdups
                

    for toks in (x.rstrip('\r\n').split("\t") for x in ts.nopen(gtf) if x[0] != "#"):
        if toks[2] not in("CDS", "stop_codon") or toks[1] not in("protein_coding"): continue
        #if toks[0] != "1": break
        start, end = map(int, toks[3:5])
        gene = toks[8].split('gene_name "')[1].split('"', 1)[0]
        assert start <= end, toks
        key = toks[0], gene

        #cutoff = 0.3

        # find sections of exon under certain coverage.
#TODO: comment this if we don't want coverage cutoff filtering
        if coverage_array[start-1:end].min() < cutoff: # doesn't bother to run these operations if there is not one bp below the cutoff
            #splitters[key].add([(start - 1, end)]) #this takes out the whole exon for one section of poor coverage
            a = coverage_array[start - 1: end]
            #print str(start-1),end,a
            is_under, locs = False, [] # generates "locs" for each exon"
            if a[0] < cutoff:
                locs.append([start - 1])
                is_under = True # so you can initialize is_under
            for pos, v in enumerate(a[1:], start=start): #enumerates positions in the coverage array starting at the beginning of the exon
                if v < cutoff:
                    if not is_under:
                        is_under = True
                        locs.append([pos - 1]) #start, coverage is in bed format, so pos-1 is necessary, since splitters are open left and right side
                else:
                    if is_under:
                        is_under = False
                        locs[-1].append(pos) #end
            if is_under:
                locs[-1].append(end) # in this case would end splitter at the end of the exon
            splitters[key].add(map(tuple, locs))
            #for i in locs:
            #    f3.write(chrom+"\t"+"\t".join(map(str,i))+"\n")

        for s, e in split_iv.find((start - 1, end)):
            splitters[key].add([(s, e)])

        genes[key].add([(start-1, end)]) # converts GTF exon coordinates to BED format (subtracts 1 from exon start)
    # sort by start so we can do binary search.
    genes = dict((k, sorted(v._vals)) for k, v in genes.iteritems())
    #ends = dict((k, sorted(v)) for k, v in ends.iteritems())
    splits, starts, ends = {}, {}, {}
    splitters = dict(splitters)
    for chrom_gene, sends in genes.iteritems():
        starts[chrom_gene] = [s[0] for s in sends]
        ends[chrom_gene] = [s[1] for s in sends]
        if chrom_gene in splitters:
            splits[chrom_gene] = splitters[chrom_gene]._vals

    return starts, ends, splits


def get_cdna_start_end(cdna_position, v):
    cdna_start, cdna_end = cdna_position.split("/")
    #cdna_end = cdna_end.rstrip("-?")
    if cdna_start[0] == "?": # deletion
        _, cdna_end = cdna_start.split("-")
        cdna_end = int(cdna_end)
        cdna_start = cdna_end - len(v.REF)
    elif "-" == cdna_start:
        cdna_start = "na"
        cdna_end = "na"
    elif "-" in cdna_start: 
        if "?" in cdna_start:
            cdna=cdna_start.split("-")
            cdna_start = int(cdna[0])
            cdna_end = "na"
        else:
            try:
                cdna_start, cdna_end = map(int, cdna_start.split("-"))
            except:
                print(v.REF, v.ALT, cdna_start, v.INFO['CSQ'])
                raise
    else:
        cdna_start = int(cdna_start)
        cdna_end = cdna_start + len(v.REF)
    return cdna_start, cdna_end

def isfunctional(csq):
    return any(c in csq['Consequence'] for c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant', 'missense_variant', 'protein_altering_variant', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion')) \
    or (('splice_donor_variant' in csq['Consequence'] or 'splice_acceptor_variant' in csq['Consequence']) and 'coding_sequence_variant' in csq['Consequence'])

def ismissense(csq):
    return any(c in csq['Consequence'] for c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant', 'missense_variant'))

def issynonymous(csq):
    return any(c in csq['Consequence'] for c in ('synonymous_variant', 'stop_retained_variant', 'start_retained_variant'))

def cg_content(seq):
    if len(seq) == 0: return 0.0
    return 2.0 * seq.count('CG') / len(seq)
