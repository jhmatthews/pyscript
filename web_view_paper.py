import os, sys

path = "/Data/matthewsj/mypapers/antithesis/bibs"


x = sys.argv[1]

url = None

entry = False

poss = []

fnames = os.listdir(path)

for ff in fnames:

	fname = "%s/%s" % (path, ff)

	f = open(fname, "r")

	for line in f:


		common = 0

		if x in line: 
			entry = True

		if entry:
			data = line.split()

			if data[0] == "adsurl":

				if url != None: print "Error: multiple matches. Going to open the last one."
				url = data[2][1:-2]

				entry = False




if url!=None:
	os.system("open -a 'Google Chrome' %s" % url)

else:
	print "None found"





