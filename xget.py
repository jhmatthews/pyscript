import os, sys

host = sys.argv[1]
file = sys.argv[2]
mode = sys.argv[3]

if host == "allx10":
	hostname = "allx10.physics.ox.ac.uk"
	home = "/home/matthewsj/"
	user = "matthewsj"

elif host == "arc":
	hostname = "arcus-b.arc.ox.ac.uk"
	home = "/home/phys-mh3d/phys1504/data/"
	user = "phys1504"

if mode == "r":
	command = "scp -r {}@{}:{}{} .".format(user, hostname, home, file)
else:
	command = "scp {}@{}:{}{} .".format(user, hostname, home, file)

print(command)
os.system(command)

