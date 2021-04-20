from collections import defaultdict
from operator import itemgetter
from copy import copy, deepcopy
import time
import heapq
import math
import sys
from functools import reduce, cmp_to_key

class Components(object):
    def __init__(self,proteins,peptides,edges,maxpepdeg=1e+20,progress=None):
        self.prot = proteins
        self.pep = peptides
        self.edges = edges
        self.invedges = defaultdict(set)
        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)
        badpep = set()
        for pep in list(self.invedges.keys()):
            if len(self.invedges[pep]) > maxpepdeg:
                badpep.add(pep)
        self.pr2comp = {}
        self.component = []
        if progress:
            progress.stage("Finding components...",len(self.prot))
        for pr in self.prot:
            if pr in self.pr2comp:
                continue
            thecomp = set()
            queue = set([pr])
            while len(queue) > 0:
                pr0 = queue.pop()
                thecomp.add(pr0)
                peps = ((self.edges[pr0] & self.pep) - badpep)
                pr1 = reduce(set.union,list(map(self.invedges.get,peps)),set())
                pr1 &= self.prot
                pr1 -= thecomp
                queue.update(pr1)
            for pr1 in thecomp:
                self.pr2comp[pr1] = len(self.component)-1
            self.component.append((thecomp,
                                   self.pep & reduce(set.union,
                                              list(map(self.edges.get,thecomp)))))
            if progress:
                progress.update()
        if progress:
            progress.done()
            progress.stage("Components: %d"%len(self.component))

    def __iter__(self):
        return next(self)

    def __next__(self):
        for pr,pep in self.component:
            yield pr,pep

def cmp(a,b):
    if a < b:
        return -1
    elif a > b:
        return 1
    return 0

class Dominator:
    def __init__(self,peptides,edges,proteins):
        self.peptides = set(peptides)
        self.proteins = set(proteins)
        self.edges = copy(edges)
        self.eliminated = set()
        self.equivalentto = defaultdict(set)
        self.containedby = defaultdict(set)
        self.invedges = defaultdict(set)
        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)

    def dominate(self,prsortkey=None,tiebreakers=None,eliminate=0,promiscuouslimit=1e+20,relationships=True,progress=None):
        pr2tbpep = defaultdict(set)
        tbpep2fdr = defaultdict(lambda :1e+20)
        tbcmp = cmp
        if tiebreakers:
            pr2tbpep = tiebreakers[0]
            tbpep2fdr = tiebreakers[1]
            tbcmp = tiebreakers[2]
        if progress:
            progress.message("Proteins: %d peptides: %d"%(len(self.proteins),len(self.peptides)))
        if prsortkey != None:
            lprots = sorted(self.proteins,key=prsortkey)
        else:
            lprots = list(self.proteins)
        # print self.proteins
        peps = set()
        for pep in self.peptides:
            if len(self.invedges[pep]) <= promiscuouslimit:
                peps.add(pep)
        pepedges = {}
        for pr in lprots:
            pepedges[pr] = (self.edges[pr]&peps)
        invedges = defaultdict(set)
        for i,pr in enumerate(lprots):
            for v in pepedges[pr]:
                invedges[v].add(i)
        contained = set()
        if progress:
            progress.stage("Checking for dominated proteins",len(lprots))
        for i in range(0,len(lprots)):
            # print i,pepedges[lprots[i]]
            if progress:
                progress.update()
            if len(pepedges[lprots[i]]) <= eliminate:
                contained.add(lprots[i])
                if relationships:
                    self.eliminated.add(lprots[i])
                continue
            # print contained
            checked = set()
            for pep in pepedges[lprots[i]]:
                # print invedges[pep], checked, contained
                for j in invedges[pep]:
                    # print j, j in checked, j in contained
                    if j == i:
                        continue
                    if j in checked:
                        continue
                    checked.add(j)
                    # if lprots[j] in contained:
                    #     continue
                    # print " ",j,pepedges[lprots[j]]
                    if pepedges[lprots[i]] >= pepedges[lprots[j]]:
                        if pepedges[lprots[i]] > pepedges[lprots[j]]:
                            if relationships:
                                self.containedby[lprots[i]].add(lprots[j])
                            contained.add(lprots[j])
                        elif prsortkey != None and i < j:
                            if relationships:
                                self.equivalentto[lprots[i]].add(lprots[j])
                                self.containedby[lprots[i]].add(lprots[j])
                            contained.add(lprots[j])
                        else:
                            tbi = sorted(map(tbpep2fdr.get,pr2tbpep[lprots[i]]-pr2tbpep[lprots[j]]),key=cmp_to_key(tbcmp))
                            tbj = sorted(map(tbpep2fdr.get,pr2tbpep[lprots[j]]-pr2tbpep[lprots[i]]),key=cmp_to_key(tbcmp))
                            tbi.extend([1e+20]*max(len(tbj)-len(tbi),0))
                            tbj.extend([1e+20]*max(len(tbi)-len(tbj),0))
                            # if len(tbi) > 0:
                            #     print >>sys.stderr, (tbi,i),(tbj,j)
                            if (tbi,i) < (tbj,j):
                                # print >>sys.stderr, (tbi,i),(tbj,j)
                                if relationships:
                                    self.equivalentto[lprots[i]].add(lprots[j])
                                    self.containedby[lprots[i]].add(lprots[j])
                                contained.add(lprots[j])

        if progress:
            progress.done()
        # print contained
        self.proteins = set(lprots) - contained
        self.edges = pepedges
        try:
            self.peptides = reduce(set.union,list(map(self.edges.get,self.proteins)))
        except TypeError:
            self.peptides = set()
        self.invedges = defaultdict(set)

        # Fix equivalence and contained relationships...
        if relationships:
            for p in self.proteins:
                self.equivalentto[p] = self.collapse(p,self.equivalentto)
                self.containedby[p] = (self.collapse(p,self.containedby) - self.equivalentto[p])
            for k in list(self.equivalentto):
                if k not in self.proteins:
                    del self.equivalentto[k]
            for k in list(self.containedby):
                if k not in self.proteins:
                    del self.containedby[k]
            # print self.proteins
            # print self.equivalentto
            # print self.containedby

        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)
        if progress:
            progress.message("Proteins: %d peptides: %d"%(len(self.proteins),len(self.peptides)))

    @staticmethod
    def collapse(p0,edges):
        collapsed = set()
        collapsed.add(p0)
        # print collapsed
        tocheck = set()
        # print p0,edges[p0]
        if p0 in edges:
            tocheck.update(edges[p0])
        while len(tocheck) > 0:
            # print collapsed,tocheck
            p1 = tocheck.pop()
            collapsed.add(p1)
            if p1 in edges:
                tocheck.update(edges[p1]-collapsed)
        return collapsed

    def force(self,forced,progress=None):
        nunique = defaultdict(int)
        if progress:
            progress.stage("Counting unique peptides per protein",len(self.peptides))
        for pep in self.peptides:
            if progress:
                progress.update()
            hitset = self.invedges[pep]
            if len(hitset) == 1:
                nunique[next(iter(hitset))] += 1
        if progress:
            progress.done()
        self.forced = set()
        if progress:
            progress.stage("Finding forced proteins",len(self.proteins))
        for pr in self.proteins:
            if progress:
                progress.update()
            if nunique[pr] >= forced:
                self.forced.add(pr)
        if progress:
            progress.done()
            progress.message("Forced proteins: %d"%len(self.forced))
        self.proteins -= self.forced
        if progress:
            progress.stage("Determining uncovered peptides",len(self.forced))
        self.covered = set()
        for pr in self.forced:
            if progress:
                progress.update()
            self.covered.update(self.edges[pr])
        if progress:
            progress.done()
            progress.message("Covered peptides: %d"%len(self.covered))
        self.peptides -= self.covered
        self.dominate(relationships=False,progress=progress)
