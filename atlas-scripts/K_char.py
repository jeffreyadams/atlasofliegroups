#!/bin/python3

#python script to compute K_characters in parallel
#output to a valid atlas file, of the form
#K_chars=[...]
#sample usage: Kpols.py -g G2_s -d 2 -f "facets_G2/facets_G2_s_dim
#compute number of fundamental facets of dimension d, then look for files
#facets_G2/facets_G2_s_dim2_i (0\le i\le d-1)
#command line arguments:
#-g: group
#-d: dimension (for naming output files)
#-f: facet_file_stem
#-o: output_file_stem (optional: override default naming convention)
#default output files: directory "group"/"group_dim" + dim + "_" + j + ".at"
#example: G2_s/G2_s_dim1_3.at

#to create a facet file (example):
#in atlas:
#>facetsE6dim4 facets(E6_s,4)   {4 dim facets}
#>facetsE6dim4 facets(E6_s)     {all facets}
import sys, time, os, getopt, subprocess, gc,re 
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp   #only for cpu count
progress_step_size=10 #how often to report progress

#simple reporting function
def report(q,i,proc):
   return("(i)[size of q]<status_i>: (" + str(i) + ")[:" + str(q.qsize()) +  "]<" + str(proc.poll()) + ">")

#call the atlas process proc, with arguments:
#q: queue of facet data, obtained from facet_file
#procs: array of processes
#i: number of process
def atlas_compute(i):
   print("starting atlas_compute process #", i)
   proc=procs[i]
   facet_file=facet_file_stem + "_" + str(i)
   print("reading data from file: ", facet_file)
   file=open(facet_file,"r")
   facets_computed_local=0  #facets computed by this process
   outputfile=group + "/" + output_file_stem + "_" + str(i) + ".at"  #atlas requires quotes if filename includes a directory
   f=open(outputfile,"wb", buffering=4000000)
   firstline=True
   for facet in file:
#      atlas_arg='{}'.format("\n prints(\"(\"," + facet +  "[0]" + ",\",\",K_data(K_char(" + group + "," + facet + "))" + ",\")\")\n" ).encode('utf-8'
      atlas_arg=b"".join([b"\n prints(\"(\",", bytes(facet,'utf-8') ,  b"[0]",b",\",\",K_data(K_char(",bytes(group,'utf-8'),b",",bytes(facet,'utf-8'),b"))",b",\")\")\n"])
#      print("atlas_arg: ", atlas_arg)
      proc.stdin.write(atlas_arg)
      proc.stdin.flush()
      line = proc.stdout.readline().decode('ascii')
      line=re.sub(r'\s+', '', line) + "\n"   #strip spaces, add newline back
      if firstline:
         print("first line of process ", i, " written at ", time.ctime())
         line="set K_chars_dim_" + str(dim) + "_" + str(i) +"=[" + line
         f.write(line.encode())
         firstline=False
      else:
         line=","+line
         f.write(line.encode())
      facets_computed_local+=1
   f.write("]".encode())
   f.close()
   return(i,facets_computed_local)

def main(argv):
   print("Starting at ", time.ctime())
   global group
   global dim
   global facet_file_stem
   global output_file_stem
   cpu_count=mp.cpu_count()
   print("----------------------------")
   print("Number of cores: ", cpu_count)
   opts, args = getopt.getopt(argv, "g:d:f:o:")
   output_file_stem_in=""
   if len(opts)==0:
      print("Usage: \n-g: group\n-d: dimension\n-f: facet_file_stem\n-o: output_file_stem (optional)")
      exit()
   for opt,arg in opts:
      print(opt, " ", arg)
   for opt, arg in opts:
       if opt in ('-g'):
          group=arg
       elif opt in ('-d'):
          dim=arg
       elif opt in ('-f'):
          facet_file_stem=arg
       elif opt in ('-o'):
          output_file_stem_in=arg
   print("group: ", group)
   print("facet dimension: ", dim)
   if output_file_stem_in !="":
      output_file_stem=output_file_stem_in
   else:
      output_file_stem=group + "_dim" + str(dim)
   if not os.path.exists(group):
      os.makedirs(group)
   #run atlas once to determine number of fundamental facets
   proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE,stdout=PIPE)
   atlas_arg='{}'.format("\n prints(#facets_fundamental(" + group + "," + str(dim) + "))\n").encode('utf-8')
#   print("atlas_arg to count fundamental facets: ", atlas_arg)
   proc.stdin.write(atlas_arg)
   proc.stdin.flush()
   number_fundamental_facets =int(proc.stdout.readline().splitlines()[0])
   print("number of fundamental facets: ", number_fundamental_facets)
   proc.kill()
   global queues
   queues=[]
   number_facets=[]
   with mp.Manager() as manager:
      global procs
      procs=[]
      T=[]   #array of results from the atlas processes
   for i in range(number_fundamental_facets):
      proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE,stdout=PIPE)
      procs.append(proc)
   print("number of processes started: ", len(procs))
   with concurrent.futures.ProcessPoolExecutor(number_fundamental_facets) as P:
      for i in range(number_fundamental_facets):
         proc=procs[i]
         Q=P.submit(atlas_compute,i)
         T.append(Q)
   total_facets=0
   print("fundamental facet #/#facets computed")
   for t in T:
      try:
         data = t.result()
#          data = (3,4)
      except Exception as exc:
         print("exception: ", exc)
#         print("data: ", data)
#         print('%r generated an exception: %s' % ("x",exc))
#         print("got data: ", data)
      else:
         (p,n)=data
         print(data)
         total_facets+=n
   print("number of facets in input files: ", number_facets)
   print("total number of facets: ", sum(number_facets))
   print("number of K_polynomials computed (sum of numbers from each process):", total_facets)

if __name__ == "__main__":
   main(sys.argv[1:])



