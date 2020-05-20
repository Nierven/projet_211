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