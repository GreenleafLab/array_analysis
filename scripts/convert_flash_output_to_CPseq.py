import sys

addrs, seqs, phreds = [],[],[]
with open(sys.argv[1],'r') as f:
	for i, lin in enumerate(f.readlines()):
		if i % 4 == 0:
			addrs.append(lin.replace('\n',''))
		elif i % 4 == 1:
			seqs.append(lin.replace('\n',''))
		elif i % 4 == 3:
			phreds.append(lin.replace('\n',''))


assert len(addrs) == len(seqs)
assert len(addrs) == len(phreds)

with open(sys.argv[2], 'w') as f:
	for i in range(len(addrs)):
		f.write('%s\t%s\t%s\n' % (addrs[i],seqs[i],phreds[i]))
