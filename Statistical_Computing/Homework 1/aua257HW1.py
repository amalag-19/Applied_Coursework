## Q2
for x in range(1,11):
  print "Iteration ", x

# Q3
def sqdiff(inpfile,whichcol1,whichcol2,outfile):
    col1=[] # initialize the vector that will hold the extracted column
    col2=[]
    if isinstance(inpfile, basestring):
      f=open(inpfile,'r') # this is a pointer to the file
    else:
      f=inpfile
    inp=f.readline() # read the line from the file
    while not(inp==""): # repeat until end of file
        # string.split: this splits the input (inp) into multiple items
        # map function: apply the function string.strip to each element returned by string.split
        # string.strip removes whitespace characters 
        outp = map(str.strip,str.split(inp)) # parse lines by blank space separators
        col1.append(outp[whichcol1-1]) # add the new column element to the list called "col"
        col2.append(outp[whichcol2-1])
        inp=f.readline() # read next line
        
    f.close() # close the file
    # col is the vector that holds the column of interest; write this out to file
    sum=0
    if (not(outfile=="")):
        g=open(outfile,"w")
        for k in range(0,(len(col1))):
            sqdiff = (float(col1[k])-float(col2[k]))**2
            sum=sum+sqdiff
            g.write(str(sqdiff)+'\n') # to write to file, convert value to string + add newlines
        g.close()

def total(inpfile):
  f=open(inpfile,'r')
  inp=f.readline()
  tot=0
  while not(inp==""):
    tot=tot+float(inp)
    inp=f.readline() # read next line
  f.close()
  print tot
    
sqdiff("pyEg1.dat",1,3,"out1.dat")

import urllib2
data = urllib2.urlopen("http://www.stat.psu.edu/~mharan/540/hwdir/pyEg2.dat")
sqdiff(data,300,650,"out2.dat")
total("out2.dat")
