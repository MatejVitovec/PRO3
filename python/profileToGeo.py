import numpy as np

data = np.loadtxt('python/nozzle.dat')

file = open('python/nozzle.geo', 'w')

Line = 'lc = 1e-2;\n'
file.write(Line)

for i in range(len(data)):

    point = 'Point(' + str(i + 1) + ') = '
    coords = '{' + str(data[i][0]) + ', ' + str(data[i][1]) + ', ' + '0, lc};' + '\n'

    file.write(point + coords)

data = np.flip(data, axis=0)

for i in range(len(data)):

    point = 'Point(' + str(i + 1 + len(data)) + ') = '
    coords = '{' + str(data[i][0]) + ', ' + str(-data[i][1]) + ', ' + '0, lc};' + '\n'

    file.write(point + coords)

file.write('\n')



Line = 'Line(1) = {'
for i in range(len(data)):
    if not(i == len(data) - 1):
        Line = Line + str(i+1) + ', '
    else:
        Line = Line + str(i+1) + '};\n'

file.write(Line)

Line = 'Line(2) = {' + str(len(data)) + ', '  + str(len(data) + 1) + '};\n'
file.write(Line)

Line = 'Line(3) = {'
for i in range(len(data)):
    if not(i == len(data) - 1):
        Line = Line + str(i + 1 + len(data)) + ', '
    else:
        Line = Line + str(i + 1 + len(data)) + '};\n'

file.write(Line)

Line = 'Line(4) = {' + str(len(data)*2) + ', 1};\n'
file.write(Line)

Line = 'Curve Loop(1) = {1, 2, 3, 4};\n'
file.write(Line)

Line = 'Plane Surface(1) = {1};\n'
file.write(Line)

Line = 'h = 0.01;\n'
Line = Line + 'Extrude {0,0,h} { Surface{1}; Layers{1}; Recombine; }\n'
file.write(Line)

file.close()