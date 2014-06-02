import subprocess
import sys
from subprocess import Popen, PIPE

def main():
    # nbm_fa = object from normalize_by_median
    cmd = ['python', 'bar_func.py']

    #ex1 = call(cmd)
    ex1 = Popen(cmd, stdout=PIPE, stderr=PIPE)

    cmd2 =['mv', 'out.txt','myout.py']
    ex2 = Popen(cmd2, stdout=PIPE, stderr=PIPE)

    cmd3 =['python', 'myout.py']
    ex2 = Popen(cmd3, stdout=PIPE, stderr=PIPE)
    com = ex2.communicate()
    print com

if __name__ == '__main__':
    main()
