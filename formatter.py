import os
import sys

def write_file(filename, folder, run_name, run_type, iso_limit):

  path = os.getcwd();
  curr_folder = path.split("/")[-1]

  print("#PBS -N " + run_name + "_" + filename)
  print("#PBS -l nodes=1:ppn=8")
  print("#PBS -l mem=10GB")
  print("#PBS -q una")
  print("#PBS -M lmthomps@hamilton.edu")
  print("#PBS -m bea")
  print("#PBS -j oe")
  print("#PBS -l walltime=12:00:00")
  print("#PBS -r y")
  print()
  
  print("echo \"\"")
  print("echo \"Beginning " + run_name + " on " +  filename + " job\"")
  print("echo \"\"")
  print()
  
  print("cd $PBS_O_WORKDIR")
  print()
  
  print("cd ../vcc-reductions/code")
  print("date")

  print("./optimized/vcc --preconfiguration=fsocial --k=2 ", end="")
  print("--solver_time_limit=36000", end="")
  print(" --decompose_limit=" + iso_limit, end="")
  print(" --run_type=\"" + run_type + "\"", end="")
  
  print(" ../../../data_sets/graphs/" + folder + "/" + filename + ".graph", end="")
  print(" >> ../../" + curr_folder + "/" + run_name + "/outputs/" + filename)

  print("echo -n \"Finished at: \"")
  print("date")
  print("echo \"\"")
  print("time")


def main():

        filename = sys.argv[1]
        folder = sys.argv[2]
        run_name = sys.argv[3]
        run_type = sys.argv[4]
        iso_limit = sys.argv[5];


        write_file(filename, folder, run_name, run_type, iso_limit);

if __name__ == "__main__":
        main()
