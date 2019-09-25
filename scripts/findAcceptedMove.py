import sys

fname = sys.argv[1]

fid = open(fname);
lines = fid.readlines();
fid.close()

fid_output = open("iter-accp.txt", 'w')
for i in range(0, len(lines)):
  if "ACCEPTED" in lines[i]:
    idx= i;
    while not ( "Iter" in lines[idx] ):
        idx-=1
    #print( lines[idx])
    splitted = lines[idx].split()
    fid_output.write( "%s\n" %splitted[2][0:-1] )

fid_output.close()
