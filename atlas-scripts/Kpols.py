#!/usr/bin/python3

#do facet computation, broken up by dimension and by facet in the FPP
#compute(G,d,[n_0,...,n_m]):
#n_i is the number of facets of dimension i (computed by number_fundamental_facets(G))
#this runs n_d atlas jobs: each one computes the KTypePols, for the facets
#of dimension d, one for each fundamental facet of dimension d
#the output of job k is to a file (for example:)
#E8_s/Kpols_E8_s_dim_1_ff_3
#the directory

import sys, getopt,os, multiprocessing, time
from multiprocessing import Process
cpu_count=multiprocessing.cpu_count()
print("Number of cores: ", cpu_count)
print("Using at most half the cores")
#argv:  F4_s 2 10,5,5,10   (just dimension 2)  OR
#       F4_s 10,5,5,10     (all dimensions)
#       group dimension #fundamental facets of all dimensions

def atlas(arg):
    group=arg[0]
    dim=arg[1]
    i=arg[2]
    print("doing atlas, dim=", dim , "  ff_number=",i)
    if not os.path.exists(group):
          os.makedirs(group)
    atlas_output_file="\"" + group + "/Kpols_" + "group_" + group + "_dim_" + str(dim) + "_ff_" + str(i) + "\""
    run_file_name="run_"+str(i)+".at"
    print("input file: ", run_file_name)
    print("output file: ", atlas_output_file)
    run_file = open(run_file_name,"w")
    run_file.write(">" + atlas_output_file + " TEST(" + group +","+ str(dim) + "," + str(i)  + ")\n")
    run_file.close()
#    print("Here is run_file_name for i=:",i)
#    with open(run_file_name, 'r') as f:
#          print(f.read())
    print(os.popen("../atlas polsSMALLEST.at < "+ run_file_name).read())
    print("done")

def main(argv):
   args=sys.argv
   print("args: ",args)
   if len(args)==4:
       group=args[1]
       dim=args[2]
       print("Only dimension: ", dim)
       ff_numbers_string=args[3]
       ff_numbers=list(map(int, ff_numbers_string.split(',')))  #this is an array of integers
       ff_number=ff_numbers[int(dim)]
       print("group: ", group)
       print("dim=",dim)
       print("ff_numbers: ", ff_numbers)
       print("ff_number: ", ff_number)
       max_cores= max(ff_number, cpu_count//2)
       print("max number of cores: ",max_cores)
       pool = multiprocessing.Pool()
       pool = multiprocessing.Pool(processes=max_cores)
       mylist=[]
       for i in range(ff_number):
             mylist.append((group,dim,i))
             print("mylist: ", mylist)
             outputs=pool.map(atlas,mylist)
   else:
       group=args[1]
       ff_numbers_string=args[2]
       ff_numbers=list(map(int, ff_numbers_string.split(',')))  #this is an array of integers
       print("all dimension: 0 to",len(ff_numbers))
       for dim in range(len(ff_numbers)):
             ff_number=ff_numbers[int(dim)]
             print("doing dimension ",dim)
             print("group: ", group)
             print("dim=",dim)
             print("ff_numbers: ", ff_numbers)
             print("ff_number: ", ff_number)
             pool = multiprocessing.Pool()
             pool = multiprocessing.Pool(processes=ff_number)
             mylist=[]
             for i in range(ff_number):
                   mylist.append((group,dim,i))
                   print("mylist: ", mylist)
                   outputs=pool.map(atlas,mylist)
             


if __name__ == "__main__":
   main(sys.argv[1:])



