#!/usr/bin/python

import random
import sys

histogram = {}
for i in xrange(10000):
	i = int(random.gauss(5, 0.001) + 0.5)
	histogram[i] = histogram.get(i, 0) + 1

m = max(histogram.values())
for k, v in sorted(histogram.items()):
	print k, '*' * (v * 50 / m)


