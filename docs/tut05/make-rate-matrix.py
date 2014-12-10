from __future__ import print_function, division

n = 20

print('row states')
lines = []
for i in range(n):
    line = ','.join(('[%d]' % i) for j in range(n) if i != j)
    lines.append('[' + line + ']')
print(',\n'.join(lines))

print('column_states')
lines = []
for i in range(n):
    line = ','.join(('[%d]' % j) for j in range(n) if i != j)
    lines.append('[' + line + ']')
print(',\n'.join(lines))

