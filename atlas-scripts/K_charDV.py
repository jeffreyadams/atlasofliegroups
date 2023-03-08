#!/bin/python3

#python script to compute K_characters in parallel
#sample usage: Kpols.py -g G2_s -f facetsG2.txt -o out -c 10
#command line arguments:
#-g: group
#-f: facet_file
#-o: output_file_root
#-c: max_cores
#this will launch max_cores processes
#process #j write output to: group/output_file_root_j
#example above: G2_s/out_1,..., G2_s/out_10

#to create a facet file (example):
#in atlas:
#>facetsE6dim4 facets(E6_s,4)   {4 dim facets}
#>facetsE6dim4 facets(E6_s)     {all facets}
import sys, time, os, getopt, subprocess, queue
import concurrent.futures
from subprocess import Popen, PIPE, STDOUT
import multiprocessing   #only for cpu count

progress_step_size=100 #how often to report progress

#simple reporting function
def report(q,i,proc):
   return("(i)[size of q]<status_i>: (" + str(i) + ")[:" + str(q.qsize()) +  "]<" + str(proc.poll()) + ">")

#call the atlas process proc, with arguments:
#q: queue of facet data, obtained from facet_file
#procs: array of processes
#i: number of process
def atlas_compute(group,output_file_root,q,procs,i):
   facets_computed_local=0  #facets computed by this process
   global facets_computed #facets computed by all processes
   proc=procs[i]
#   outputfile="\"" + group + "/" + output_file_root + "_" + str(i)  + "\""   #atlas requires quotes around the filename if it includes a directory
   outputfile=group + "/" + output_file_root + "_" + str(i)
   f=open(outputfile,"w")
   while not q.empty():
      if q.qsize()%report_interval==0:
         report_string=str(q.qsize())
         print("size of q: "+ " "*(facets_digits - len(report_string)) + report_string,end="\r")
      report(q,i,proc)
      try:
         facet = q.get(False)
      except queue.Empty:
         quit_arg='{}'.format("\n quit").encode('utf-8')
         proc.stdin.write(quit_arg)
         print("quitting", i)
         report(q,i,proc)
      else:
#         atlas_arg='{}'.format("\n prints(\"(\"," + facet +  "[0]" + ",\",\",K_data(K_char(" + group + "," + facet + "))" + ",\"),\")\n" ).encode('utf-8')
         atlas_arg='{}'.format("\n >>\""+ outputfile + "\" prints(" + facet + ",K_data(K_char(" + group + "," + facet + "))" + ",\"),\")\n prints(\"HA\")\n").encode('utf-8')#         print(atlas_arg)
         z=proc.stdin.write(atlas_arg)
#         print("z=",z)
         proc.stdin.flush()
         line = proc.stdout.readline().decode('ascii')
         if not line == "HA\n":
            print("output line:", line)
#         f=open(outputfile,"a")
#         f.write(line.decode('ascii'))
#      print("DONE")

      facets_computed+=1  #total facets computed by all processed
      facets_computed_local+=1
   return(i,facets_computed_local)

def main(argv):
   global facets_computed
   global report_interval
   global facets_digits
   report_interval=100
   facets_computed=0
   cpu_count=multiprocessing.cpu_count()
   print("----------------------------")
   print("Number of cores: ", cpu_count)
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
   global facets_input
   facets_input=q.qsize()
   facets_digits=len(str(q.qsize()))
   print("number of facets in input file: ", facets_input)
#   outputfile=group + "/" + output_file_root + "_"
#   f=open(outputfile,"w")
   procs=[]
   T=[]   #array of results from the atlas processes
   print("starting ", max_cores, "atlas processes")
   for i in range(max_cores):
      filename=  group + "/" + output_file_root + "_" + str(i)
      with open(filename, "w") as outfile:
         proc=subprocess.Popen(["../atlas","polsParallel.at"], stdin=PIPE,stdout=PIPE)
         procs.append(proc)
   print("number of processes started: ", len(procs))
   P=concurrent.futures.ThreadPoolExecutor()
   for i in range(max_cores):

         proc=procs[i]
         T.append(P.submit(atlas_compute,group,output_file_root,q,procs,i))
   total_facets=0
   print("list of process numbers/number of facets computed")
   for t in T:
      try:
         data = t.result()
      except Exception as exc:
         print('%r generated an exception: %s' % ("x",exc))
      else:
         (p,n)=data
         total_facets+=n
         if p is 0:
            print()
         print(data)
   print("number of facets in input file: ", facets_input)
   print("number of K_polynomials computed: ", total_facets, "(", facets_computed,")")

#   print("status of processes: ")
#   for i in range(len(procs)):
#      print("i: ", i, procs[i].poll())
#   for p in procs:
#      print(proc.poll())
#   exit()

if __name__ == "__main__":
   main(sys.argv[1:])



