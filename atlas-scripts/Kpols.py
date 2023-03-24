#!/bin/python3

#python script to compute K_characters in parallel
#sample usage: Kpols.py -g G2_s -f facetsG2.txt -o out -c 10
#command line arguments:
#-g: group
#-f: facet_file
#-o: output_file_root
#-c: max_cores
#this will launch max_cores processes
#projess #j write output to: group/output_file_root_j
#example above: G2_s/out_1,..., G2_s/out_10

import sys, time, os, getopt, subprocess, queue
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing   #only for cpu count

cpu_count=multiprocessing.cpu_count()
print("Number of cores: ", cpu_count)


#call the atlas process running on proc, with arguments:
#output_file_root (see command line options)
#q: queue of facet data, obtained from facet_file
#procs: array of processes
#i: number of process
def report(q,i,proc):
   return("(i)[size of q]<status_i>: (" + str(i) + ")[:" + str(q.qsize()) +  "]<" + str(proc.poll()) + ">")
   
def atlas_compute(group,output_file_root,q,procs,i):
   count=0
   print("in atlas compute")
   print("group=", group)
   print("output_file_root=", output_file_root)
   proc=procs[i]
   file_name=  group + "/" + output_file_root + "_" + str(i)
   print("before while: ",report(q,i,proc))
   while not q.empty():
      count+=1
      print("inside while: ", report(q,i,proc))

      try:
         facet = q.get(False)
         print("got facet=", facet, " i=", i)
      except Queue.Empty:
         print("q is empty, i=", i)
         proc.stdin.write('{}'.format("\n quit").encode('utf-8'))
#         time.sleep(.2)
         report(q,i,proc)
      else:
         atlas_arg='{}'.format("\n  prints(\"(\",(" + str(i) + ")," + facet +  "[0]" + ",K_data(K_char(" + group + "," + facet + ")))").encode('utf-8')
         proc.stdin.write(atlas_arg)

      
#      facet=q.get()
#      atlas_arg='{}'.format("\n >> " + "\"" + group + "/" + output_file_root + "_" + str(i) + "\""   + " prints(\"(\"," + facet +  "[0]" + ",K_data(K_char(" + group + "," + facet + ")))").encode('utf-8')

#      atlas_arg='{}'.format("\n >> " + "\"" + file_name +  "\""   + " prints(\"(\"," + facet +  "[0]" + ",K_data(K_char(" + group + "," + facet + ")))").encode('utf-8')
#      atlas_arg='{}'.format("\n  prints(\"(\",(" + str(i) + ")," + facet +  "[0]" + ",K_data(K_char(" + group + "," + facet + ")))").encode('utf-8')
#      print("after stdin.write (facet): (", facet,")", report(q,i,proc),"\n")
#      time.sleep(.2)
#      print("done process:",i)

#   print("done running process #",i)
   print("outside while", report(q,i,proc))

#   proc.kill()
#   print("killed: ", i," :",  proc.poll())
#
#   print("before kill", report(q,i,proc))
#   killresult=proc.kill()
#   print("killresult: ", killresult)
#   print("after kill", report(q,i,proc))
#   time.sleep(.4)
#   proc.stdin.write('{}'.format("\n quit").encode('utf-8'))
#   killresult=proc.kill()
   time.sleep(.5)
#  print("after sleep", report(q,i,proc))
   return(i,count)


#usage:
#Kpols.py E8_s facetsE8_dim_4 KcharsE8_dim_4.at 
def main(argv):
   opts, args = getopt.getopt(argv, "g:c:f:o:")
   if len(opts)==0:
      print("Usage: \n-g: group\n-f: facet_file\n-o: output_file_root\n-c: number_of_cores\n")
      exit()
   for opt,arg in opts:
      print(opt, " ", arg)
   for opt, arg in opts:
       if opt in ('-g'):
          group=arg
       elif opt in ('-f'):
          facet_file=arg
       elif opt in ('-o'):
          output_file_root=arg
       elif opt in ('-c'):
          max_cores=int(arg)
   print("group: ", group)
   print("facet_file: ", facet_file)
   print("output_file_root: ", output_file_root)
   if not os.path.exists(group):
      os.makedirs(group)
   q=queue.Queue()
   file=open(facet_file,"r")
   #read facet_file and put each entry on q
   data=file.read().splitlines()
   for d in data:
      q.put(d)
   facets_input=q.qsize()
   print("number of facets in input file: ", facets_input)
   #initialize array of max_cores atlas processes
   Q=concurrent.futures.ThreadPoolExecutor()
   procs=[]
   T=[]   #array of results from the atlas processes
   for i in range(max_cores):
      filename=  group + "/" + output_file_root + "_" + str(i)
      print("filename: ", filename)
      with open(filename, "w") as outfile:
        proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE, stdout=outfile)
        procs.append(proc)
   print("number of procs: ", len(procs))
#   for P in procs:
#        print(P.poll())
   P=concurrent.futures.ThreadPoolExecutor()
   print("q=",q.qsize())
#   T.append(P.submit(atlas_compute,"G2_s","out",q,procs[0],1))

#   print("T=",T)
   for i in range(max_cores):
         print("submitting job i=",i)
         proc=procs[i]
         T.append(P.submit(atlas_compute,group,output_file_root,q,procs,i))
         print("appended proc")
   total_facets=0
   print("list of process numbers/number of facets computed")
   print("T: ", len(T))

   for t in T:
      try:
         data = t.result()
      except Exception as exc:
         print('%r generated an exception: %s' % ("x",exc))
      else:
         (p,n)=data
         total_facets+=n
         print(data)
   print("number of facets in input file: ", facets_input)
   print("number of K_polynomials computed: ", total_facets)
   print("status of processes: ")
   for i in range(len(procs)):
        if procs[i].poll() == None:
             print("still running: ", i)
   exit()
   
if __name__ == "__main__":
   main(sys.argv[1:])



