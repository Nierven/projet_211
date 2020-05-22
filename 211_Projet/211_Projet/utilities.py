from math import atan2, sqrt
from sys import platform
from os import environ, popen
from copy import deepcopy

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __sub__(self, p):
        return Point(self.x - p.x, self.y - p.y)

    def distance(self, p):
        dist_x = p.x - self.x
        dist_y = p.y - self.y
        return sqrt(dist_x**2 + dist_y**2)

    def angle(self, p):
        dist_x = p.x - self.x
        dist_y = p.y - self.y
        return arctan2(dist_y, dist_x)

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

    def __str__(self):
        return '({};{})'.format(self.x, self.y)

    def __repr__(self):
        return '<point>: ' + str(self)

linewidth = 45
def print_header(str):
    l = len(str)
    side_l = int((linewidth - l - 2) / 2)

    offset = 0
    if (l % 2 != 0 and linewidth % 2 == 0 or 
        l % 2 == 0 and linewidth % 2 != 0):
        offset = 1

    print()
    print('=' * linewidth)
    print('=' * (side_l + offset) + ' ' + str + ' ' + '=' * side_l)
    print('=' * linewidth)
    print()

# Print iterations progress
def print_progress_bar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def get_threads_count():
    """ Returns the number of available threads on a posix/win based system"""
    if platform == 'win32':
        return (int) (environ['NUMBER_OF_PROCESSORS'])
    else:
        return (int) (popen('grep -c cores /proc/cpuinfo').read())